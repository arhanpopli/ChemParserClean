# Chrome Extension mol2chemfig Options Implementation

## Overview
Added complete support for mol2chemfig-specific rendering options in the chemistry extension, allowing users to choose between MoleculeViewer (port 5000) and mol2chemfig (port 8000) renderers with distinct option sets for each.

## Changes Made

### 1. popup.html
**Added**: Complete mol2chemfig options section (lines ~350-450)

**Options included**:
- Show Carbons (-c flag)
- Aromatic Circles (-o flag)
- Show Methyls (-m flag)
- Fancy Bonds (-f flag)
- Atom Numbers (-n flag)
- Compact Mode (-z flag)
- Flip Horizontal (client-side CSS transform)
- Flip Vertical (client-side CSS transform)
- Hydrogen Treatment (keep/add/delete dropdown)

**ID naming convention**: All mol2chemfig options prefixed with `m2cf` to differentiate from MoleculeViewer options (e.g., `m2cfShowCarbonsToggle` vs `showCarbonsToggle`)

### 2. popup.js
**Added DOM references** (lines ~24-33):
```javascript
const mol2chemfigOptions = document.getElementById('mol2chemfigOptions');
const m2cfShowCarbonsToggle = document.getElementById('m2cfShowCarbonsToggle');
// ... 9 total option elements
```

**Added default settings** (lines ~57-65):
```javascript
m2cfShowCarbons: false,
m2cfAromaticCircles: false,
m2cfShowMethyls: false,
m2cfFancyBonds: false,
m2cfAtomNumbers: false,
m2cfCompact: false,
m2cfFlipHorizontal: false,
m2cfFlipVertical: false,
m2cfHydrogensMode: 'keep'
```

**Added settings loader** (lines ~85-93):
```javascript
if (m2cfShowCarbonsToggle) m2cfShowCarbonsToggle.checked = settings.m2cfShowCarbons;
// ... repeated for all 9 options
```

**Added event listeners** (lines ~215-295):
- Each option saves to chrome.storage.sync when changed
- Displays status message to user
- Pattern: `chrome.storage.sync.set({ m2cfShowCarbons: e.target.checked })`

**Updated `updateEngineInfo()` function**:
```javascript
if (engine === 'moleculeviewer') {
  // Show MoleculeViewer options, hide mol2chemfig options
} else if (engine === 'mol2chemfig') {
  // Hide MoleculeViewer options, show mol2chemfig options
}
```

### 3. content.js
**Added mol2chemfig settings to defaults** (lines ~431-440):
```javascript
// mol2chemfig rendering options
m2cfShowCarbons: false,
m2cfAromaticCircles: false,
m2cfShowMethyls: false,
m2cfFancyBonds: false,
m2cfAtomNumbers: false,
m2cfCompact: false,
m2cfFlipHorizontal: false,
m2cfFlipVertical: false,
m2cfHydrogensMode: 'keep'
```

**Updated API payload construction** (lines ~1008-1020):
Changed from incorrect `options` object to correct `selections` array:
```javascript
// Build selections array for mol2chemfig API
const selections = [];
if (settings.m2cfAromaticCircles) selections.push('-o');  // Aromatic circles
if (settings.m2cfShowCarbons) selections.push('-c');       // Show carbon labels
if (settings.m2cfShowMethyls) selections.push('-m');       // Show methyl labels
if (settings.m2cfFancyBonds) selections.push('-f');        // Fancy bonds
if (settings.m2cfAtomNumbers) selections.push('-n');       // Atom numbers
if (settings.m2cfCompact) selections.push('-z');           // Compact structure

fetch(`${MOL2CHEMFIG_API}/m2cf/submit`, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({ 
    textAreaData: inputData,
    selections: selections,
    h2: settings.m2cfHydrogensMode || 'keep'
  })
})
```

**Added client-side flip transforms** (3 locations in SVG rendering):
```javascript
// Apply client-side flip transforms for mol2chemfig
let transform = '';
if (settings.m2cfFlipHorizontal) transform += 'scaleX(-1) ';
if (settings.m2cfFlipVertical) transform += 'scaleY(-1) ';

svgImg.style.cssText = `
  display: inline-block;
  max-width: 300px;
  max-height: 200px;
  margin: 0 12px 8px 0;
  vertical-align: middle;
  cursor: pointer;
  ${transform ? `transform: ${transform.trim()};` : ''}
`;
```

Applied to:
1. Data URI svglink (line ~1118)
2. Raw SVG Blob URL (line ~1138)
3. Final svgImg creation (line ~1195)

## API Integration

### mol2chemfig API Endpoint: `/m2cf/submit`
**Correct payload format**:
```javascript
{
  textAreaData: "CCO",           // SMILES, MOL block, or ChemFig
  selections: ['-o', '-c'],      // Array of command-line flags
  h2: "keep"                     // Hydrogen treatment
}
```

**Command-line flags**:
- `-o` = Aromatic circles
- `-c` = Show carbon labels
- `-m` = Show methyl labels
- `-f` = Fancy bonds
- `-n` = Atom numbers
- `-z` = Compact structure

**Response format**: `{ chemfig, svglink, pdflink, error }`
- `svglink`: Inline SVG XML (NOT a URL)
- Priority: Display SVG (structural diagram), fallback to PDF

## User Workflow

1. **Open extension popup**: Click extension icon
2. **Choose renderer**: Select radio button
   - 🧪 MoleculeViewer (localhost:5000)
   - 📐 mol2chemfig (localhost:8000)
3. **Configure options**: Toggle options specific to selected renderer
4. **Apply**: Reload page to render formulas with new settings

### MoleculeViewer Options
- Aromatic Circles
- Show Carbons
- Show Methyls
- Fancy Bonds
- Atom Numbers
- Flip Horizontal/Vertical
- Hydrogen Mode (keep/add/delete)

### mol2chemfig Options
- Show Carbons (-c)
- Aromatic Circles (-o)
- Show Methyls (-m)
- Fancy Bonds (-f)
- Atom Numbers (-n)
- Compact Mode (-z)
- Flip Horizontal (client-side)
- Flip Vertical (client-side)
- Hydrogen Mode (keep/add/delete)

## Testing Checklist

### UI Tests
- [ ] Open extension popup, verify both option sections exist
- [ ] Select MoleculeViewer → verify moleculeViewerOptions visible, mol2chemfigOptions hidden
- [ ] Select mol2chemfig → verify mol2chemfigOptions visible, moleculeViewerOptions hidden
- [ ] Toggle each option → verify status message displays
- [ ] Close and reopen popup → verify settings persist

### MoleculeViewer Renderer Tests
- [ ] Select MoleculeViewer, enable aromatic circles → reload page with benzene SMILES
- [ ] Enable show carbons → verify C labels appear
- [ ] Enable flip horizontal → verify structure mirrors
- [ ] Change hydrogen mode → verify hydrogen display changes

### mol2chemfig Renderer Tests
- [ ] Select mol2chemfig, enable aromatic circles → reload with benzene
- [ ] Enable show carbons → verify C labels appear on structure
- [ ] Enable show methyls → verify CH3 groups labeled
- [ ] Enable compact mode → verify tighter structure
- [ ] Enable flip horizontal → verify immediate mirror effect (CSS transform)
- [ ] Enable flip vertical → verify upside-down effect
- [ ] Test combinations: aromatic circles + show carbons + fancy bonds

### API Integration Tests
- [ ] Verify `/m2cf/submit` receives correct `selections` array
- [ ] Check console logs for API payload
- [ ] Verify SVG renders as structural diagram (not text)
- [ ] Test with SMILES: `CCO`, `c1ccccc1`, `CC(=O)O`
- [ ] Test error handling (invalid SMILES)

### Edge Cases
- [ ] Switch renderers while page has formulas → verify reload prompt
- [ ] Rapidly toggle options → verify no race conditions
- [ ] Test on various websites (ChatGPT, Wikipedia chemistry pages)
- [ ] Test with malformed chemical notation

## Implementation Notes

### Why Separate Option Sets?
- MoleculeViewer and mol2chemfig use different backend APIs with different parameter formats
- MoleculeViewer: Uses query parameters on GET request
- mol2chemfig: Uses POST JSON payload with `selections` array of flags
- User-facing names similar, but backend integration differs

### Why Client-Side Flips?
User specified: "flip vertical and flip horizontal should be done in the extension itself"
- **Advantage**: Instant visual feedback (no server round-trip)
- **Advantage**: No re-rendering required (faster)
- **Implementation**: CSS `transform: scaleX(-1)` and `scaleY(-1)`
- **Applied**: During SVG image element creation in content.js

### Hydrogen Treatment
Both renderers support hydrogen display modes:
- **keep**: Show explicit hydrogens as specified in input
- **add**: Add all implicit hydrogens
- **delete**: Remove all hydrogens

MoleculeViewer parameter: `hydrogensMode`
mol2chemfig parameter: `h2`

## File Locations

```
c:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\chem-extension\
├── manifest.json          # Extension manifest
├── popup.html             # Settings UI (✅ UPDATED)
├── popup.js               # Settings logic (✅ UPDATED)
├── content.js             # Formula rendering (✅ UPDATED)
├── background.js          # Service worker (no changes)
└── IMPLEMENTATION_NOTES.md # This file
```

## Related Web Interface Files

```
c:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\templates\
├── index.html             # Main web UI with ChemFig textarea (✅ FIXED)
└── m2cf.html              # Standalone mol2chemfig component (✅ FIXED)
```

**Key fix**: Changed `applyM2CFOptions()` to handle inline SVG (not fetch as URL):
```javascript
// OLD (BROKEN): tried to fetch svglink as URL
const svgResponse = await fetch(result.svglink);

// NEW (FIXED): handles inline SVG XML
if (result.svglink && result.svglink.trim().startsWith('<')) {
  containerHTML = `<div>${result.svglink}</div>`;
}
```

## Next Steps

1. **Test Extension**: Load unpacked extension in Chrome, verify all options work
2. **Test Web Interface**: Open localhost:5000, test ChemFig Output textarea and options
3. **Cross-Test**: Verify web and extension use same backend correctly
4. **Performance Test**: Render 10+ formulas on page, check load times
5. **Edge Case Testing**: Invalid SMILES, server errors, network failures
6. **Documentation**: Update user-facing README with option descriptions

## Known Issues / Future Enhancements

- **Dark Mode**: SVG color inversion logic only in content.js (not web interface)
- **Zoom Controls**: Consider adding scale slider for SVG sizes
- **Batch Rendering**: Optimize when many formulas on single page
- **Keyboard Shortcuts**: Add hotkeys for common options
- **Export**: Add download SVG/PNG buttons in extension popup
- **Caching**: Consider caching rendered SVGs to reduce API calls

## Dependencies

- **Backend**: mol2chemfig Docker container (port 8000) must be running
- **Backend**: MoleculeViewer Flask server (port 5000) must be running
- **Browser**: Chrome/Edge with Manifest V3 support
- **API**: Both backends must be accessible via localhost HTTP (not HTTPS)
