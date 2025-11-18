# SVG Caching System - Quick Reference

## What Was Implemented

You now have an **intelligent SVG caching system** that:

1. ✅ **Generates unique cache links** for each molecule + option combination
2. ✅ **Stores SVGs with descriptive names** (not just random hashes)
3. ✅ **Tracks source** (Mol2ChemFig vs MoleculeViewer)
4. ✅ **Automatically cleans up** old cache files after 24 hours
5. ✅ **Reuses cached SVGs** if same molecule + options requested again

## Cache Filename Format

```
{molecule}_{options}_{source}_{hash}.svg
```

### Real Examples:

**Benzene (default):**
- `benzene_default_m2cf_a1b2c3d4.svg`

**Caffeine with aromatic circles:**
- `caffeine_aromatic_m2cf_def67890.svg`

**Aspirin with aromatic + carbons + methyls:**
- `aspirin_aromatic_carbons_methyls_m2cf_ghi13579.svg`

**Ibuprofen (MoleculeViewer) with fancy bonds:**
- `ibuprofen_fancy_mv_jkl24680.svg`

**Benzene rotated 90°:**
- `benzene_rot90_m2cf_mno35791.svg`

## How It Works

### 1. Generate a Molecule
```javascript
// User enters SMILES in Mol2ChemFig tab
"C1=CC=CC=C1"  (benzene)

// With options: aromatic_circles = true
// Generates SVG and automatically caches it
```

### 2. Backend Creates Cache File
```
File created: /svg-cache/benzene_aromatic_m2cf_abc12345.svg
```

### 3. User Sees Cache Link
```
📍 Cache: /cache/benzene_aromatic_m2cf_abc12345.svg
```

### 4. User Can:
- ✅ Copy and share the link
- ✅ Access the same SVG again without regenerating
- ✅ Use in documentation or presentations
- ✅ Reference the exact options used (in the filename!)

## Key Features

### 🎯 Smart Option Tracking
Options are encoded in filename so you can see what was used:
- `aromatic` = Aromatic circles enabled
- `carbons` = Show carbon atoms
- `methyls` = Show methyl groups
- `atoms` = Atom numbering
- `fancy` = Fancy bond rendering
- `rot90` = Rotated 90 degrees
- `h_add` = Hydrogens added
- `indent4` = Indentation set to 4

### 🚀 Reusability
Same molecule + options = Same cache file

```
First time: "caffeine" + "aromatic" 
→ Generates SVG
→ Saves to: caffeine_aromatic_m2cf_abc123.svg

Second time: "caffeine" + "aromatic"
→ REUSES same cache file!
→ No regeneration needed
```

### 🧹 Automatic Cleanup
- Cache files expire after 24 hours
- Old files are automatically deleted
- Prevents disk space issues

### 🔀 Source Separation
- `m2cf` = Mol2ChemFig (LaTeX format)
- `mv` = MoleculeViewer (standard format)

Same molecule can have different SVGs from each source!

```
caffeine_aromatic_m2cf_abc123.svg  ← Mol2ChemFig version
caffeine_aromatic_mv_def456.svg    ← MoleculeViewer version
```

## File Locations

**Cache Storage:** `C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\svg-cache\`

**Access URL:** `http://localhost:5000/cache/{filename}`

**API Endpoints:**
- POST `/api/cache-svg` - Cache an SVG
- GET `/cache/<filename>` - Serve cached SVG
- GET `/cache/info` - View cache statistics
- POST `/cache/cleanup` - Force cleanup

## Usage Examples

### Example 1: Benzene with Different Options

**Generate benzene (default):**
```
Input: C1=CC=CC=C1
Options: (none)
↓
Cache: benzene_default_m2cf_a1b2.svg
```

**Generate same benzene with aromatic circles:**
```
Input: C1=CC=CC=C1
Options: aromatic_circles = true
↓
Cache: benzene_aromatic_m2cf_c3d4.svg
```

**Different cache files, different URLs:**
- `/cache/benzene_default_m2cf_a1b2.svg`
- `/cache/benzene_aromatic_m2cf_c3d4.svg`

### Example 2: Caffeine from Search

**Search for "caffeine":**
```
Step 1: Nomenclature search → Converts to SMILES
Step 2: Generates SVG with options
Step 3: Caches as: caffeine_aromatic_carbons_m2cf_e5f6.svg
Step 4: Shows cache URL: /cache/caffeine_aromatic_carbons_m2cf_e5f6.svg
```

**Change to MoleculeViewer and search again:**
```
Step 1: Same "caffeine" SMILES
Step 2: Different rendering (MoleculeViewer format)
Step 3: Caches as: caffeine_aromatic_carbons_mv_g7h8.svg
```

**Result: Two different cache files for same molecule!**

## Testing the System

### 1. Open the web interface
```
http://localhost:5000/
```

### 2. Go to Mol2ChemFig tab

### 3. Try these examples:
```
Entry 1: C1=CC=CC=C1 (Benzene)
Options: Default
Cache file: benzene_default_m2cf_*.svg

Entry 2: CC(=O)O (Acetic Acid)
Options: aromatic + carbons + methyls
Cache file: acetic_acid_aromatic_carbons_methyls_m2cf_*.svg

Entry 3: caffeine (Chemical name search)
Options: aromatic + atoms
Cache file: caffeine_aromatic_atoms_m2cf_*.svg
```

### 4. View cache directory
```
C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\svg-cache\
```

All your generated SVGs are stored there with descriptive names!

## Backend Implementation Details

### Modified Files:

1. **`app/cache_manager.py`** (NEW)
   - Cache key generation
   - SVG saving with smart naming
   - Cleanup functions

2. **`app/api.py`** (UPDATED)
   - Imports cache manager
   - `/api/cache-svg` endpoint (NEW)
   - `/cache/info` endpoint (NEW)
   - `/api/smiles-to-svg` returns cache_url (UPDATED)

3. **`templates/index.html`** (UPDATED)
   - `cacheM2CFSVG()` function (NEW)
   - Displays cache URLs (UPDATED)
   - Integrates caching into workflows

## Configuration

Edit `app/cache_manager.py` to customize:

```python
# Cache expiry time (hours)
CACHE_EXPIRY_HOURS = 24

# Cache directory location
CACHE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'svg-cache'))
```

## Troubleshooting

### Q: Where are the cached SVGs stored?
**A:** In `/svg-cache/` folder. Check: `MoleculeViewer/svg-cache/`

### Q: Why does the same molecule have different cache files?
**A:** Because the options are different! Check the filename to see what options were used.

### Q: Can I delete cache files manually?
**A:** Yes, they're just regular SVG files in the svg-cache folder.

### Q: How long do cache files last?
**A:** 24 hours. Old files are automatically cleaned up.

### Q: What's the hash at the end of the filename?
**A:** Content hash to ensure uniqueness. Two SVGs with identical molecule+options still get unique files if content differs.

## Next Steps

You can now:
1. ✅ Generate molecules with various options
2. ✅ See cache links for each generation
3. ✅ Share cache URLs with others
4. ✅ Access cached SVGs directly via links
5. ✅ Monitor what options were used (in filename!)
