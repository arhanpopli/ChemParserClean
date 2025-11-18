# Options & Cache URL Implementation Summary

## Problem
The rendering options (aromatic circles, show carbons, show methyls, etc.) were not being sent to the backend servers. Both MoleculeViewer and mol2chemfig were ignoring user settings, always rendering with default options.

## Solution Overview
Implemented a **query-parameter-based options system** where:
1. Extension sends options as URL parameters to the server
2. Server renders SVG with the specified options
3. Server caches the SVG with **options encoded in the filename**
4. Each unique combination of (SMILES/nomenclature + options) gets its own cached file

## Changes Made

### 1. MoleculeViewer Backend (`app/api.py`)

#### `/img/smiles` endpoint (lines ~639-690)
**Before**: Always called `smiles_to_svg(smiles, width, height, {})`  
**After**: Parses options from query parameters and passes them to rendering

```python
# Parse rendering options from query parameters
options = {
    'show_carbons': request.args.get('show_carbons', 'false').lower() == 'true',
    'show_methyls': request.args.get('show_methyls', 'false').lower() == 'true',
    'aromatic_circles': request.args.get('aromatic_circles', 'false').lower() == 'true',
    'fancy_bonds': request.args.get('fancy_bonds', 'false').lower() == 'true',
    'atom_numbers': request.args.get('atom_numbers', 'false').lower() == 'true',
    'flip_horizontal': request.args.get('flip_horizontal', 'false').lower() == 'true',
    'flip_vertical': request.args.get('flip_vertical', 'false').lower() == 'true',
    'hydrogens_mode': request.args.get('hydrogens_mode', 'keep')
}

error, svg = smiles_to_svg(smiles, width, height, options)
```

**Cache Key with Options**:
```python
# Create cache key with options encoded in filename
options_str = '_'.join([
    k for k, v in options.items() 
    if v is True or (k == 'hydrogens_mode' and v != 'keep')
])
cache_identifier = f"smiles_{smiles[:10]}_{options_str}" if options_str else f"smiles_{smiles[:10]}"
```

**Example Cache URLs**:
- Default: `smiles_c1ccccc1_a1b2c3.svg`
- With aromatic circles: `smiles_c1ccccc1_aromatic_circles_a1b2c3.svg`
- With multiple options: `smiles_c1ccccc1_aromatic_circles_show_carbons_a1b2c3.svg`

#### `/img/nomenclature` endpoint (lines ~692-750)
Same changes applied - parses options, passes to rendering, encodes in cache filename.

### 2. Chrome Extension (`content.js`)

#### MoleculeViewer API Calls (lines ~800-845)
**Before**: Simple URL with only SMILES/nomenclature  
**After**: Uses `URLSearchParams` to build complete options query string

```javascript
// Build options query string for MoleculeViewer
const optionsParams = new URLSearchParams({
  smiles: moleculeData.smiles,
  width: '300',
  height: '200',
  json: 'true',
  show_carbons: settings.showCarbons.toString(),
  show_methyls: settings.showMethyls.toString(),
  aromatic_circles: settings.aromaticCircles.toString(),
  fancy_bonds: settings.fancyBonds.toString(),
  atom_numbers: settings.atomNumbers.toString(),
  flip_horizontal: settings.flipHorizontal.toString(),
  flip_vertical: settings.flipVertical.toString(),
  hydrogens_mode: settings.hydrogensMode
});

apiUrl = `${MOLECULE_VIEWER_API}/img/smiles?${optionsParams.toString()}`;
```

**Example API Request**:
```
http://localhost:5000/img/smiles?smiles=c1ccccc1&width=300&height=200&json=true&show_carbons=true&aromatic_circles=true&fancy_bonds=false&...
```

#### mol2chemfig API Calls (lines ~1008-1025)
Already updated in previous session - uses `selections` array with flags.

```javascript
const selections = [];
if (settings.m2cfAromaticCircles) selections.push('-o');
if (settings.m2cfShowCarbons) selections.push('-c');
if (settings.m2cfShowMethyls) selections.push('-m');
// ... etc

fetch(`${MOL2CHEMFIG_API}/m2cf/submit`, {
  method: 'POST',
  body: JSON.stringify({ 
    textAreaData: inputData,
    selections: selections,
    h2: settings.m2cfHydrogensMode
  })
})
```

## How It Works

### Workflow for MoleculeViewer:

1. **User enables options** in extension popup (e.g., "Aromatic Circles")
2. **Extension saves** to `chrome.storage.sync`
3. **Content script reads** settings when rendering formulas
4. **Content script builds URL** with all options as query parameters
5. **Flask backend receives** request, parses options
6. **Backend renders SVG** with specified options using RDKit
7. **Backend generates cache filename** encoding the options:
   - `smiles_CCO_aromatic_circles_show_carbons_abc123.svg`
8. **Backend saves SVG** to `svg-cache/` directory
9. **Backend returns JSON** with cache URL
10. **Extension displays SVG** from cache URL

### Cache Behavior:

**First Request** (benzene with aromatic circles):
```
GET /img/smiles?smiles=c1ccccc1&aromatic_circles=true
→ Renders SVG with aromatic circles
→ Saves to: smiles_c1ccccc1_aromatic_circles_xyz789.svg
→ Returns: { cache_url: "http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_xyz789.svg" }
```

**Second Request** (benzene with default options):
```
GET /img/smiles?smiles=c1ccccc1
→ Renders SVG WITHOUT aromatic circles
→ Saves to: smiles_c1ccccc1_abc123.svg (DIFFERENT FILE)
→ Returns: { cache_url: "http://localhost:5000/cache/smiles_c1ccccc1_abc123.svg" }
```

**Third Request** (benzene with aromatic circles again):
```
GET /img/smiles?smiles=c1ccccc1&aromatic_circles=true
→ Uses EXISTING cached file: smiles_c1ccccc1_aromatic_circles_xyz789.svg
→ Returns same cache URL
```

### Cache Expiration:
- Files expire after **24 hours** (`CACHE_EXPIRY_HOURS = 24`)
- Cleanup runs on backend startup
- Old files automatically deleted to save disk space

## Testing

### Manual Testing:
1. Open `test_options_cache.html` in browser
2. Click "Run Tests" for each molecule
3. Observe:
   - Different cache URLs for different options
   - Same cache URL when repeating same options
   - Visual differences in rendered structures

### Test Cases:
- ✅ Benzene default vs aromatic circles (should see circles)
- ✅ Propane with show_methyls (should see CH3 labels)
- ✅ Benzene with show_carbons (should see C labels)
- ✅ Combinations (aromatic_circles + show_carbons)
- ✅ Flip horizontal/vertical (should mirror structure)
- ✅ Adrenaline nomenclature with options

### Extension Testing:
1. Load extension in Chrome
2. Open webpage with chemical formulas (e.g., Wikipedia chemistry page)
3. Toggle options in extension popup
4. Reload page
5. Open DevTools Console → check API URLs being called
6. Verify console logs show: `Cache URL: http://localhost:5000/cache/smiles_...`

## Benefits

### 1. True Options Support
- Options actually affect rendering (previously ignored)
- Each option combination renders correctly
- Server-side rendering ensures consistent quality

### 2. Efficient Caching
- No redundant renders for same (molecule + options)
- Disk-based cache survives server restarts
- Automatic cleanup prevents cache bloat

### 3. Debugging & Transparency
- Cache URLs show exactly what options were used
- Filename encoding makes debugging easy
- Console logs track cache hits/misses

### 4. Scalability
- Can add new options without breaking existing cache
- Options encoded in URL = RESTful API design
- Cache shareable across multiple clients

## Files Modified

```
c:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\
├── MoleculeViewer\
│   ├── app\
│   │   └── api.py                          # ✅ Updated both img endpoints
│   └── test_options_cache.html             # ✅ New test page
│
└── chem-extension\
    ├── content.js                          # ✅ Updated MoleculeViewer API calls
    ├── popup.html                          # ✅ (From previous session)
    └── popup.js                            # ✅ (From previous session)
```

## API Documentation

### GET `/img/smiles`

**Query Parameters**:
- `smiles` (required): SMILES notation string
- `width` (optional): Image width (default: 300)
- `height` (optional): Image height (default: 200)
- `json` (optional): Return JSON (default: true)
- `show_carbons` (optional): Show carbon labels (true/false)
- `show_methyls` (optional): Show methyl labels (true/false)
- `aromatic_circles` (optional): Draw aromatic circles (true/false)
- `fancy_bonds` (optional): Enhanced bond rendering (true/false)
- `atom_numbers` (optional): Show atom indices (true/false)
- `flip_horizontal` (optional): Mirror horizontally (true/false)
- `flip_vertical` (optional): Flip vertically (true/false)
- `hydrogens_mode` (optional): Hydrogen display (keep/add/delete)

**Response JSON**:
```json
{
  "success": true,
  "smiles": "c1ccccc1",
  "cache_url": "http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_abc123.svg",
  "image_url": "http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_abc123.svg",
  "expires_in_hours": 24,
  "svg": "<svg>...</svg>"
}
```

### GET `/img/nomenclature`

**Query Parameters**:
- `nomenclature` (required): Chemical name (e.g., "benzene", "adrenaline")
- Same option parameters as `/img/smiles`

**Response JSON**:
```json
{
  "success": true,
  "nomenclature": "benzene",
  "smiles": "c1ccccc1",
  "cache_url": "http://localhost:5000/cache/nomenclature_benzene_aromatic_circles_xyz789.svg",
  "image_url": "...",
  "expires_in_hours": 24,
  "svg": "<svg>...</svg>"
}
```

## Common Issues & Solutions

### Issue: Options not applying
**Symptom**: All renders look the same regardless of options  
**Solution**: Check browser console for API URL - ensure options are in query string  
**Debug**: Look for `show_carbons=true` etc. in URL

### Issue: Cache not working
**Symptom**: Server renders every time, even for identical requests  
**Solution**: Check `svg-cache/` directory permissions  
**Debug**: Look for "Cached SVG: filename.svg" in server logs

### Issue: Old cached versions showing
**Symptom**: Options changed but still seeing old rendering  
**Solution**: Clear `svg-cache/` directory or wait 24 hours  
**Debug**: Check cache filenames - options should be in filename

### Issue: Extension options not syncing
**Symptom**: Options toggle in popup but don't affect rendering  
**Solution**: Reload webpage after changing options  
**Debug**: Check `chrome.storage.sync` in DevTools Application tab

## Future Enhancements

### 1. Cache Management UI
Add admin endpoint to:
- View all cached files
- Clear cache manually
- See cache hit/miss statistics

### 2. Cache Prewarming
Pre-render common molecules with popular option combinations

### 3. CDN Integration
Serve cached SVGs from CDN for faster worldwide access

### 4. Compression
Gzip cached SVG files to save disk space

### 5. Analytics
Track which options are most popular for UX improvements

## Verification Steps

### Backend Verification:
```bash
# 1. Check server is running with new code
curl "http://localhost:5000/img/smiles?smiles=CCO&aromatic_circles=true"

# 2. Check cache directory
ls MoleculeViewer/svg-cache/

# 3. Verify cache filename includes options
# Should see files like: smiles_CCO_aromatic_circles_abc123.svg
```

### Extension Verification:
1. Open extension popup → enable "Aromatic Circles"
2. Open ChatGPT, type "benzene"
3. Press F12 → Console tab
4. Look for: `API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true&...`
5. Look for: `Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_...`

### Visual Verification:
- **Aromatic Circles**: Should see circle inside benzene ring instead of alternating double bonds
- **Show Carbons**: Should see "C" labels on carbon atoms
- **Show Methyls**: Should see "CH3" labels on methyl groups
- **Flip**: Structure should be mirrored or upside-down

## Conclusion

The options system now works end-to-end:
- ✅ Options sent from extension to backend
- ✅ Backend renders with correct options
- ✅ Cache URLs encode options in filename
- ✅ Each option combination cached separately
- ✅ Visual differences visible in rendered structures

**Test**: Visit test page to see all options in action!  
**URL**: http://localhost:5000/test_options_cache.html
