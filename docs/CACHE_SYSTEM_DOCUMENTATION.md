# SVG Caching System Implementation

## Overview
Created an intelligent SVG caching system that stores generated molecular structures with meaningful, descriptive filenames based on the molecule and rendering options applied.

## How It Works

### 1. Cache Key Generation (`app/cache_manager.py`)
Generates descriptive cache filenames using molecule identifiers and options:

**Examples:**
- `benzene_default_m2cf_a1b2c3d4.svg` - Benzene with default options from Mol2ChemFig
- `caffeine_aromatic_carbons_methyls_mv_x7y8z9k0.svg` - Caffeine with aromatic circles, show carbons, show methyls from MoleculeViewer
- `aspirin_fancy_rotate90_m2cf_q1w2e3r4.svg` - Aspirin with fancy bonds, rotated 90°

### 2. File Structure
```
/cache/{name}_{options}_{source}_{hash}.svg
│
├── {name}          : Sanitized molecule name or SMILES prefix
│                    (e.g., "benzene", "caffeine", "acetic_acid")
│
├── {options}       : Comma-separated option tags
│                    - Aromatic circles: "aromatic"
│                    - Show carbons: "carbons"
│                    - Show methyls: "methyls"
│                    - Show atom numbers: "atoms"
│                    - Fancy bonds: "fancy"
│                    - Hydrogen mode: "h_add", "h_delete"
│                    - Rotation: "rot90", "rot180"
│                    - Indentation: "indent4", "indent8"
│                    - Default (no options): "default"
│
├── {source}        : Source application
│                    - "m2cf" = Mol2ChemFig
│                    - "mv" = MoleculeViewer
│
└── {hash}          : 8-character content hash (ensures uniqueness)
```

### 3. Backend Endpoints

#### `/api/cache-svg` (POST)
Save SVG to cache with intelligent naming

**Request:**
```json
{
    "svg_content": "<svg>...</svg>",
    "smiles_or_name": "caffeine",
    "options": {
        "aromatic_circles": true,
        "show_carbons": true
    },
    "source": "mol2chemfig"
}
```

**Response:**
```json
{
    "success": true,
    "cache_url": "/cache/caffeine_aromatic_carbons_m2cf_a1b2c3d4.svg",
    "filename": "caffeine_aromatic_carbons_m2cf_a1b2c3d4.svg"
}
```

#### `/cache/<filename>` (GET)
Serve cached SVG files

#### `/cache/info` (GET)
Get cache statistics (file count, total size, expiry)

#### `/cache/cleanup` (POST)
Manually trigger cleanup of expired cache files (>24 hours old)

### 4. MoleculeViewer Integration (`/api/smiles-to-svg`)
SVGs are automatically cached when converting SMILES to SVG:

```python
# Backend automatically caches the SVG
cache_url, cache_path = save_svg_to_cache(
    svg,
    smiles,
    options,
    source="moleculeviewer",
    base_url=request.host_url.rstrip('/')
)

# Returns cache_url in response
return jsonify({
    'error': None,
    'svg': svg,
    'cache_url': cache_url,  # NEW
    'smiles': smiles,
    'info': info
}), 200
```

### 5. Frontend Integration (`templates/index.html`)

#### JavaScript Function: `cacheM2CFSVG()`
```javascript
async function cacheM2CFSVG(svgContent, smiles, options = {}, source = 'mol2chemfig') {
    // Sends SVG to backend for caching
    // Returns the cache URL
}
```

#### Display Cache URL
When SVGs are rendered, the cache URL is displayed below the image:

```
📍 Cache: /cache/caffeine_aromatic_m2cf_a1b2c3d4.svg
```

Users can:
- Copy the cache link
- Share it directly (the SVG is permanent for 24 hours)
- Use it in documentation or reports

## Features

✅ **Descriptive Filenames** - Easy to identify what molecule and options were used
✅ **Content-Based Hashing** - Ensures uniqueness even if same options applied
✅ **Automatic Cleanup** - Removes old cache files after 24 hours
✅ **Cross-Source Tracking** - Distinguishes between Mol2ChemFig and MoleculeViewer
✅ **Option Tracking** - Preserves all rendering options in filename
✅ **Reusable Links** - Cache URLs can be shared and reused
✅ **Multiple Formats** - Handles various option combinations intelligently

## Example Workflow

### Scenario: Generate caffeine structure with different options

**Step 1:** Enter "caffeine" in search → SVG generated
- **Cache file:** `caffeine_default_m2cf_abc12345.svg`
- **Link:** `/cache/caffeine_default_m2cf_abc12345.svg`

**Step 2:** Add aromatic circles option → SVG regenerated
- **Cache file:** `caffeine_aromatic_m2cf_def67890.svg`
- **Link:** `/cache/caffeine_aromatic_m2cf_def67890.svg`

**Step 3:** Add show methyls → SVG regenerated
- **Cache file:** `caffeine_aromatic_methyls_m2cf_ghi13579.svg`
- **Link:** `/cache/caffeine_aromatic_methyls_m2cf_ghi13579.svg`

**Step 4:** Switch to MoleculeViewer tab, enter same SMILES
- **Cache file:** `caffeine_default_mv_jkl24680.svg`
- **Link:** `/cache/caffeine_default_mv_jkl24680.svg`

All four SVGs are stored separately and can be accessed via their unique cache URLs.

## Cache Manager API (`app/cache_manager.py`)

### Functions:

**`create_cache_key(smiles_or_name, options, source)`**
- Creates descriptive cache filename from molecule and options
- Returns: `"caffeine_aromatic_carbons_m2cf"`

**`save_svg_to_cache(svg_content, smiles_or_name, options, source, base_url)`**
- Saves SVG to cache directory
- Returns: `(cache_url, local_filepath)`

**`get_cached_svg(smiles_or_name, options, source)`**
- Checks if SVG already exists in cache
- Returns: filename if found, None otherwise

**`cleanup_old_cache()`**
- Removes files older than `CACHE_EXPIRY_HOURS` (default: 24)

**`get_cache_stats()`**
- Returns dictionary with cache statistics
- Includes: file count, total size, expiry hours, cache directory

## Storage Location
Cache files are stored in: `/MoleculeViewer/svg-cache/`

Each file is accessible via HTTP at: `http://localhost:5000/cache/{filename}`

## Configuration
Located in `app/cache_manager.py`:
- `CACHE_DIR` - Directory to store SVGs
- `CACHE_EXPIRY_HOURS` - Hours before cache cleanup (default: 24)

## Example Cache Filenames Breakdown

```
caffeine_aromatic_carbons_atoms_m2cf_d34e9021.svg
│         │         │       │     │        │
└─────────┴─────────┴───────┘     │        └─ Content hash (8 chars)
    Molecule + Options             │
                                    └─ Source (m2cf = Mol2ChemFig)

phenol_default_mv_c8f4b0e2.svg
│      │       │  │
└──┬───┘       │  └─ Hash
   │           └─ Source (mv = MoleculeViewer)
   └─ Molecule

CC(C)Cc1ccc_fancy_rotate90_m2cf_7a5d6f88.svg
│            │     │          │
└────┬────────┘     │          └─ Source
     │              └─ Options
     └─ SMILES prefix (first 10 chars)
```

## Future Enhancements

- [ ] Database logging of cache access
- [ ] Cache statistics dashboard
- [ ] Batch caching for multiple molecules
- [ ] Cache size limits with LRU eviction
- [ ] Compressed cache archive export
- [ ] Cache versioning for option changes
