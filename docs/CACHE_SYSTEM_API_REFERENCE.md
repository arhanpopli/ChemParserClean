# SVG Caching System - API Reference

## Overview
Complete API reference for the SVG caching system integrated with MoleculeViewer.

---

## Backend Endpoints

### 1. Cache SVG (Store a New SVG)
**Endpoint:** `POST /api/cache-svg`

**Purpose:** Send an SVG to the backend for caching with intelligent naming

**Request Body:**
```json
{
    "svg_content": "<svg xmlns='http://www.w3.org/2000/svg'...>...</svg>",
    "smiles_or_name": "caffeine",
    "options": {
        "aromatic_circles": true,
        "show_carbons": false,
        "show_methyls": true,
        "fancy_bonds": true,
        "atom_numbers": false,
        "hydrogens": "keep",
        "rotate": 0,
        "indentation": 4
    },
    "source": "mol2chemfig"
}
```

**Response (Success):**
```json
{
    "success": true,
    "cache_url": "/cache/caffeine_aromatic_methyls_fancy_m2cf_a1b2c3d4.svg",
    "filename": "caffeine_aromatic_methyls_fancy_m2cf_a1b2c3d4.svg"
}
```

**Response (Error):**
```json
{
    "success": false,
    "error": "SVG content and molecule identifier required"
}
```

**Status Codes:**
- `200` - SVG successfully cached
- `400` - Missing required fields
- `500` - Server error

**Notes:**
- `source` should be either `"mol2chemfig"` or `"moleculeviewer"`
- `options` object is flexible; only include options that are actually set
- Content is automatically hashed to ensure uniqueness

---

### 2. Serve Cached SVG (Retrieve an SVG)
**Endpoint:** `GET /cache/<filename>`

**Purpose:** Download or display a cached SVG file

**Example URL:**
```
http://localhost:5000/cache/caffeine_aromatic_m2cf_a1b2c3d4.svg
```

**Response:**
- SVG content with `Content-Type: image/svg+xml`
- Can be displayed directly in browser or downloaded

**Status Codes:**
- `200` - SVG found and served
- `400` - Invalid filename (contains .. or /)
- `404` - File not found

**Security:**
- Filenames are validated to prevent directory traversal attacks
- Only files in the cache directory can be accessed

---

### 3. Get Cache Info (View Statistics)
**Endpoint:** `GET /cache/info`

**Purpose:** Get information about cache usage

**Response:**
```json
{
    "success": true,
    "cache_info": {
        "cache_dir": "C:\\Users\\Kapil\\Personal\\PROJECTS\\Mol2chemfig\\MoleculeViewer\\svg-cache",
        "file_count": 42,
        "total_size_mb": 3.25,
        "expiry_hours": 24
    }
}
```

**Status Codes:**
- `200` - Statistics retrieved
- `500` - Error retrieving stats

---

### 4. Cleanup Cache (Remove Old Files)
**Endpoint:** `POST /cache/cleanup`

**Purpose:** Manually trigger removal of cache files older than expiry time

**Response:**
```json
{
    "status": "cleaned"
}
```

**What It Does:**
- Scans cache directory
- Removes any files older than `CACHE_EXPIRY_HOURS` (default: 24)
- Logs cleaned filenames to console
- Automatic cleanup also runs on server startup

---

### 5. Convert SMILES to SVG (Updated)
**Endpoint:** `POST /api/smiles-to-svg`

**Purpose:** Convert SMILES notation to SVG with options

**Request Body:**
```json
{
    "smiles": "C1=CC=CC=C1",
    "width": 600,
    "height": 500,
    "options": {
        "show_carbons": false,
        "show_methyls": false,
        "aromatic_circles": true,
        "fancy_bonds": true,
        "atom_numbers": false,
        "hydrogens": "keep",
        "flip_horizontal": false,
        "flip_vertical": false,
        "rotate": 0,
        "recalculate_coordinates": false
    }
}
```

**Response:**
```json
{
    "error": null,
    "svg": "<svg xmlns='http://www.w3.org/2000/svg'>...</svg>",
    "cache_url": "/cache/benzene_aromatic_mv_a1b2c3d4.svg",
    "smiles": "C1=CC=CC=C1",
    "info": {
        "formula": "C6H6",
        "molecular_weight": 78.11,
        "num_atoms": 6,
        "num_bonds": 6,
        "canonical_smiles": "c1ccccc1",
        "inchi": "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
    }
}
```

**New Field:** `cache_url` - URL to the cached SVG

---

## Frontend Functions

### JavaScript: cacheM2CFSVG()
**Purpose:** Send an SVG to the backend for caching

**Signature:**
```javascript
async function cacheM2CFSVG(svgContent, smiles, options = {}, source = 'mol2chemfig')
```

**Parameters:**
- `svgContent` (string) - SVG markup to cache
- `smiles` (string) - Molecule SMILES or name
- `options` (object) - Rendering options applied
- `source` (string) - "mol2chemfig" or "moleculeviewer"

**Returns:**
- Cache URL (string) if successful
- `null` if failed

**Example Usage:**
```javascript
const svgText = '<svg>...</svg>';
const cacheUrl = await cacheM2CFSVG(
    svgText,
    'caffeine',
    { aromatic_circles: true, show_carbons: true },
    'mol2chemfig'
);

if (cacheUrl) {
    console.log('SVG cached at:', cacheUrl);
    // Display cache URL to user
} else {
    console.warn('Failed to cache SVG');
}
```

---

## Complete Workflow Example

### Frontend Request Flow:

```
User enters "C1=CC=CC=C1" in search box
        ↓
Frontend calls submitM2CFMolecule()
        ↓
POST to http://localhost:8000/m2cf/submit
        ↓
Docker returns SVG content + chemfig code
        ↓
Frontend calls cacheM2CFSVG(svgContent, 'benzene', options, 'mol2chemfig')
        ↓
JavaScript POST to /api/cache-svg
        ↓
Backend creates: benzene_aromatic_m2cf_a1b2c3d4.svg
        ↓
Returns: /cache/benzene_aromatic_m2cf_a1b2c3d4.svg
        ↓
Frontend displays SVG + cache link
        ↓
User sees: 📍 Cache: /cache/benzene_aromatic_m2cf_a1b2c3d4.svg
```

---

## Cache Key Generation Logic

**Function:** `create_cache_key(smiles_or_name, options, source)`

### How Options Are Encoded:

1. **Aromatic Circles**
   - Option: `aromatic_circles: true`
   - Tag: `aromatic`

2. **Show Carbons**
   - Option: `show_carbons: true`
   - Tag: `carbons`

3. **Show Methyls**
   - Option: `show_methyls: true`
   - Tag: `methyls`

4. **Atom Numbers**
   - Option: `atom_numbers: true`
   - Tag: `atoms`

5. **Fancy Bonds**
   - Option: `fancy_bonds: true`
   - Tag: `fancy`

6. **Hydrogens Mode**
   - Option: `hydrogens: "add"`
   - Tag: `h_add` or `h_delete` (if not "keep")

7. **Rotation**
   - Option: `rotate: 90`
   - Tag: `rot90`

8. **Indentation** (Mol2ChemFig)
   - Option: `indentation: 4`
   - Tag: `indent4`

9. **Angle** (Mol2ChemFig)
   - Option: `angle: 30`
   - Tag: `angle30`

10. **H2 Setting** (Mol2ChemFig)
    - Option: `h2: "mcf"`
    - Tag: `h2_mcf`

### Examples:

```
Input:  smiles="C1=CC=CC=C1", options={}, source="mol2chemfig"
Output: benzene_default_m2cf

Input:  smiles="caffeine", options={"aromatic_circles": true}, source="mol2chemfig"
Output: caffeine_aromatic_m2cf

Input:  smiles="aspirin", options={"fancy_bonds": true, "rotate": 90}, source="moleculeviewer"
Output: aspirin_fancy_rot90_mv

Input:  smiles="ibuprofen", options={"show_carbons": true, "show_methyls": true, "atom_numbers": true}, source="mol2chemfig"
Output: ibuprofen_carbons_methyls_atoms_m2cf
```

---

## Error Handling

### Common Errors:

**Invalid molecule:**
```json
{
    "success": false,
    "error": "SVG content and molecule identifier required"
}
```

**File not found:**
```
GET /cache/nonexistent_file.svg
→ 404 Not Found
```

**Directory traversal attempt:**
```
GET /cache/../../../etc/passwd
→ 400 Invalid filename
```

---

## Configuration

### Edit Cache Settings:

File: `app/cache_manager.py`

```python
# Cache expiry time (in hours)
CACHE_EXPIRY_HOURS = 24

# Cache directory (relative to app directory)
CACHE_DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'svg-cache'
))
```

---

## Performance Notes

- **Cache Lookup:** O(n) where n = number of files in cache
- **Cache Saving:** ~10-50ms depending on SVG size
- **Cleanup:** Runs on server startup + on-demand via POST
- **Storage:** Typical molecule SVG = 5-50KB

---

## Integration Checklist

✅ Backend cache manager implemented
✅ API endpoints created
✅ Frontend caching function added
✅ Display cache URLs in UI
✅ Automatic cleanup on startup
✅ Content-based hashing for uniqueness
✅ Source tracking (mol2chemfig vs moleculeviewer)
✅ Option encoding in filenames

---

## Troubleshooting

### Cache URL not showing?
- Check browser console for errors
- Verify `/api/cache-svg` endpoint is accessible
- Ensure SVG content is valid

### Cache files not being saved?
- Check file permissions on svg-cache directory
- Verify disk space available
- Check backend logs for errors

### Old cache files not cleaning up?
- Call POST `/cache/cleanup` manually
- Check `CACHE_EXPIRY_HOURS` setting
- Verify cron/scheduler if auto-cleanup configured

### Getting 404 on cached SVG?
- Verify filename in cache directory
- Check filename doesn't contain special characters
- Ensure file extension is `.svg`
