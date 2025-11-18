# Cache Structure Visual Diagram

## Before (Problems)

```
Chemparser/
├── MoleculeViewer/
│   ├── server.js                          (Port 5000)
│   └── svg-cache/                         ⚠️ MIXED FILES!
│       ├── 87970ea81928442541e33f0907da38bc.svg     ← MoleculeViewer file
│       ├── c_c_o_default_m2cf_4736b45e.svg         ← mol2chemfig file ❌
│       ├── c1ccccc1_default_m2cf_8ef9b140.svg      ← mol2chemfig file ❌
│       └── b23662e00ac02e4cc994bea7ba97448b.svg     ← MoleculeViewer file
│
├── mol2chemfig_server.py                  (Port 5001)
├── mol2chemfig_storage/                   (Empty - not being used!)
│
└── backend_data/                          (Docker internal use)
```

**Problems**:
- ❌ mol2chemfig files ending up in MoleculeViewer cache
- ❌ Confusion about which files belong to which system
- ❌ `mol2chemfig_storage/` was defined but not being used
- ❌ Hard to debug cache issues
- ❌ Hard to clear cache for one system without affecting the other

---

## After (Solution) ✅

```
Chemparser/
├── MoleculeViewer/
│   ├── server.js                          (Port 5000)
│   ├── cache/
│   │   └── moleculeviewer/                ✅ ONLY MoleculeViewer files
│   │       ├── 87970ea81928442541e33f0907da38bc.svg
│   │       └── b23662e00ac02e4cc994bea7ba97448b.svg
│   │
│   └── svg-cache/                         (Old - to be removed after testing)
│
├── mol2chemfig_server.py                  (Port 5001)
├── cache/
│   └── mol2chemfig/                       ✅ ONLY mol2chemfig files
│       ├── c_c_o_default_m2cf_4736b45e.svg
│       └── c1ccccc1_default_m2cf_8ef9b140.svg
│
├── mol2chemfig_storage/                   (Old - to be removed after testing)
│
└── backend_data/                          (Docker internal use - unchanged)
```

**Benefits**:
- ✅ Clear separation - each system has its own cache
- ✅ No file mixing between systems
- ✅ Easy to identify which cache belongs to which system
- ✅ Simple to clear cache for one system
- ✅ Better organization and maintenance

---

## Cache File Format Comparison

### MoleculeViewer Files
```
Location: MoleculeViewer/cache/moleculeviewer/
Format:   {md5_hash}.svg
Example:  87970ea81928442541e33f0907da38bc.svg
Hash of:  type:smiles:widthxheight
```

### mol2chemfig Files
```
Location: cache/mol2chemfig/
Format:   {sha256_hash}.{ext}
Example:  c_c_o_default_m2cf_4736b45e.svg
Hash of:  smiles:options
```

---

## Data Flow

### MoleculeViewer (Port 5000)
```
Request: GET /img/smiles?smiles=CCO
   ↓
Check: MoleculeViewer/cache/moleculeviewer/{hash}.svg
   ↓
   ├─ Found → Return cached SVG ✅
   │
   └─ Not found → Generate new SVG
                  ↓
                  Python: generate_svg.py (RDKit)
                  ↓
                  Save to: MoleculeViewer/cache/moleculeviewer/
                  ↓
                  Return new SVG ✅
```

### mol2chemfig (Port 5001)
```
Request: POST /api/generate {"smiles":"CCO"}
   ↓
Check: cache/mol2chemfig/{hash}.svg
   ↓
   ├─ Found → Return cached SVG ✅
   │
   └─ Not found → Call Docker backend (port 8000)
                  ↓
                  Backend: mol2chemfig LaTeX → SVG
                  ↓
                  Save to: cache/mol2chemfig/
                  ↓
                  Return new SVG ✅
```

---

## Server Configuration

### MoleculeViewer (server.js)
```javascript
const CACHE_DIR = path.join(__dirname, 'cache', 'moleculeviewer');
//                           ↑         ↑       ↑
//                           │         │       └─ Subfolder
//                           │         └─────────── Cache root
//                           └───────────────────── MoleculeViewer directory
```

### mol2chemfig (mol2chemfig_server.py)
```python
STORAGE_DIR = Path("cache") / "mol2chemfig"
//            ↑      ↑         ↑
//            │      │         └─ Subfolder
//            │      └─────────── Cache root (project level)
//            └────────────────── Path object
```

---

## Cache Management Endpoints

### MoleculeViewer
```
GET    /cache-info     → View cache statistics
DELETE /clear-cache    → Clear all cached SVGs
```

### mol2chemfig
```
GET  /api/cache/stats  → View cache statistics
POST /api/cache/clear  → Clear all cached SVGs
```

---

## Testing Commands

### Generate Test Files
```bash
# MoleculeViewer (creates file in MoleculeViewer/cache/moleculeviewer/)
curl "http://localhost:5000/img/smiles?smiles=CCO"

# mol2chemfig (creates file in cache/mol2chemfig/)
curl -X POST http://localhost:5001/api/generate \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO"}'
```

### Verify Separation
```bash
# Check MoleculeViewer cache
ls MoleculeViewer/cache/moleculeviewer/

# Check mol2chemfig cache
ls cache/mol2chemfig/

# Run automated tests
python test_cache_separation.py
```

---

## Migration Path

### Step 1: Current State
```
✓ Old code using old cache locations
✓ Files in old locations
```

### Step 2: Update Code (Done!)
```
✓ server.js updated
✓ mol2chemfig_server.py updated
✓ .gitignore updated
```

### Step 3: Test New System
```
⏳ Start servers with new code
⏳ Generate test files
⏳ Verify files go to new locations
⏳ Run test_cache_separation.py
```

### Step 4: Backup Old Caches (Optional)
```
⏳ Rename old directories to .backup
⏳ Keep for a few days/weeks
```

### Step 5: Clean Up
```
⏳ Delete old backup directories
⏳ System fully migrated ✅
```

---

## File Tree (Complete)

```
Chemparser/
├── MoleculeViewer/
│   ├── server.js                          # Updated ✅
│   ├── generate_svg.py                    # Unchanged
│   ├── nomenclature_to_smiles.py          # Unchanged
│   ├── cache/                             # NEW ✅
│   │   └── moleculeviewer/                # NEW ✅
│   │       └── {hash}.svg                 # Generated SVGs
│   └── svg-cache/                         # OLD (to remove)
│
├── mol2chemfig_server.py                  # Updated ✅
├── cache/                                 # NEW ✅
│   └── mol2chemfig/                       # NEW ✅
│       ├── {hash}.svg                     # Generated SVGs
│       └── {hash}.pdf                     # Generated PDFs
│
├── mol2chemfig_storage/                   # OLD (to remove)
├── backend_data/                          # Docker (unchanged)
├── .gitignore                             # Updated ✅
├── CLAUDE_CONTEXT.md                      # Updated ✅
│
└── Documentation (NEW):
    ├── CACHE_SEPARATION_IMPLEMENTATION.md # Full docs
    ├── CACHE_QUICK_REFERENCE.md           # Quick ref
    ├── CACHE_SEPARATION_SUMMARY.md        # Summary
    ├── CACHE_STRUCTURE_DIAGRAM.md         # This file
    └── test_cache_separation.py           # Test script
```

---

## Quick Reference

| Aspect | MoleculeViewer | mol2chemfig |
|--------|----------------|-------------|
| **Port** | 5000 | 5001 |
| **Cache Location** | `MoleculeViewer/cache/moleculeviewer/` | `cache/mol2chemfig/` |
| **File Format** | `{md5}.svg` | `{sha256}.{svg/pdf}` |
| **Rendering Engine** | RDKit | mol2chemfig (LaTeX) |
| **Cache Info** | `GET /cache-info` | `GET /api/cache/stats` |
| **Clear Cache** | `DELETE /clear-cache` | `POST /api/cache/clear` |
| **Start Command** | `cd MoleculeViewer && node server.js` | `python mol2chemfig_server.py` |

---

**Last Updated**: 2025-11-09
**Status**: Implementation complete, testing pending
