# Cache Separation Implementation

## Overview
This document describes the separated cache folder structure for the Chemparser project, addressing the issue mentioned in `MoleculeViewer/docs/Todolist.md` line 14.

## Problem Statement
Previously, there was potential overlap and confusion in cache organization:
- MoleculeViewer was using `MoleculeViewer/svg-cache/`
- mol2chemfig wrapper was using `mol2chemfig_storage/`
- Docker backend was using `backend_data/`
- Investigation revealed that mol2chemfig files were being mixed into MoleculeViewer's cache

## Solution Implemented

### New Cache Structure
```
Chemparser/
├── MoleculeViewer/
│   └── cache/
│       └── moleculeviewer/      # MoleculeViewer cache (ONLY RDKit-generated SVGs)
│
└── cache/
    └── mol2chemfig/              # mol2chemfig wrapper cache (ONLY mol2chemfig SVGs/PDFs)
```

### Folder Purposes

#### 1. MoleculeViewer Cache: `MoleculeViewer/cache/moleculeviewer/`
- **Server**: MoleculeViewer Node.js server (port 5000)
- **File**: `MoleculeViewer/server.js`
- **Contents**: SVG files generated via RDKit from SMILES and nomenclature
- **File Format**: MD5 hash filenames (e.g., `87970ea81928442541e33f0907da38bc.svg`)
- **Key Features**:
  - Aromatic circles
  - Fancy bonds
  - RDKit rendering engine

#### 2. mol2chemfig Cache: `cache/mol2chemfig/`
- **Server**: mol2chemfig Flask wrapper (port 5001)
- **File**: `mol2chemfig_server.py`
- **Contents**: SVG and PDF files generated via mol2chemfig from SMILES
- **File Format**: SHA256 hash filenames (e.g., `a1b2c3d4e5f6g7h8.svg`, `a1b2c3d4e5f6g7h8.pdf`)
- **Key Features**:
  - ChemFig LaTeX rendering
  - Multiple export formats (SVG, PDF)
  - Custom mol2chemfig options support

#### 3. Docker Backend Data: `backend_data/`
- **Server**: mol2chemfig Docker backend (port 8000)
- **Purpose**: Internal Docker backend temporary/working data
- **Note**: This is NOT a cache for the wrapper server, it's for Docker internal use

## Code Changes

### 1. MoleculeViewer Server (`MoleculeViewer/server.js`)
```javascript
// OLD
const CACHE_DIR = path.join(__dirname, 'svg-cache');

// NEW
const CACHE_DIR = path.join(__dirname, 'cache', 'moleculeviewer');
```

### 2. mol2chemfig Server (`mol2chemfig_server.py`)
```python
# OLD
STORAGE_DIR = Path("mol2chemfig_storage")

# NEW
STORAGE_DIR = Path("cache") / "mol2chemfig"
```

### 3. .gitignore
Added comprehensive cache exclusions:
```gitignore
# Cache - separate folders for each system
cache/
svg-cache/
mol2chemfig_storage/
*.cache
```

### 4. CLAUDE_CONTEXT.md
Updated documentation to reflect new cache locations:
- MoleculeViewer cache: `MoleculeViewer/cache/moleculeviewer/`
- mol2chemfig wrapper cache: `cache/mol2chemfig/`
- Docker backend data: `backend_data/` (internal use only)

## Benefits

### 1. Clear Separation
- Each system has its own dedicated cache folder
- No mixing of files from different rendering engines
- Easy to identify which cache belongs to which system

### 2. No Conflicts
- MoleculeViewer cache is within its own directory structure
- mol2chemfig cache is at the project root level
- Different hashing algorithms prevent filename collisions

### 3. Easier Maintenance
- Clear cache for specific system without affecting others
- Easy to monitor cache size per system
- Debugging is simpler with separated files

### 4. Better Organization
- Logical folder structure
- Self-documenting (folder names indicate purpose)
- Follows best practices for multi-service applications

## Cache Management

### View Cache Info

**MoleculeViewer**:
```bash
curl http://localhost:5000/cache-info
```

**mol2chemfig wrapper**:
```bash
curl http://localhost:5001/api/cache/stats
```

### Clear Cache

**MoleculeViewer**:
```bash
curl -X DELETE http://localhost:5000/clear-cache
```

**mol2chemfig wrapper**:
```bash
curl -X POST http://localhost:5001/api/cache/clear
```

### Manual Cache Cleanup

**MoleculeViewer**:
```bash
# Windows
rmdir /s /q "MoleculeViewer\cache\moleculeviewer"
mkdir "MoleculeViewer\cache\moleculeviewer"

# Linux/Mac
rm -rf MoleculeViewer/cache/moleculeviewer/*
```

**mol2chemfig**:
```bash
# Windows
rmdir /s /q "cache\mol2chemfig"
mkdir "cache\mol2chemfig"

# Linux/Mac
rm -rf cache/mol2chemfig/*
```

## Migration Notes

### Old Cache Cleanup
After confirming the new system works:

1. **Backup old caches** (optional):
   ```bash
   # Windows
   move MoleculeViewer\svg-cache MoleculeViewer\svg-cache.backup
   move mol2chemfig_storage mol2chemfig_storage.backup

   # Linux/Mac
   mv MoleculeViewer/svg-cache MoleculeViewer/svg-cache.backup
   mv mol2chemfig_storage mol2chemfig_storage.backup
   ```

2. **Delete old caches** (after testing):
   ```bash
   # Windows
   rmdir /s /q MoleculeViewer\svg-cache.backup
   rmdir /s /q mol2chemfig_storage.backup

   # Linux/Mac
   rm -rf MoleculeViewer/svg-cache.backup
   rm -rf mol2chemfig_storage.backup
   ```

## Testing

### Test MoleculeViewer Cache
```bash
# Start server
cd MoleculeViewer
node server.js

# Generate SVG (should create cache/moleculeviewer/)
curl "http://localhost:5000/img/smiles?smiles=CCO"

# Verify cache
curl http://localhost:5000/cache-info
```

### Test mol2chemfig Cache
```bash
# Start server
python mol2chemfig_server.py

# Generate SVG (should create cache/mol2chemfig/)
curl -X POST http://localhost:5001/api/generate \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO"}'

# Verify cache
curl http://localhost:5001/api/cache/stats
```

## Verification Checklist

- [x] MoleculeViewer uses `MoleculeViewer/cache/moleculeviewer/`
- [x] mol2chemfig wrapper uses `cache/mol2chemfig/`
- [x] No file mixing between caches
- [x] Both servers create their cache directories automatically
- [x] .gitignore excludes all cache directories
- [x] Documentation updated (CLAUDE_CONTEXT.md)
- [ ] Both servers tested and working
- [ ] Cache operations (read/write) verified
- [ ] Old cache folders cleaned up (after testing)

## Troubleshooting

### Cache Not Created
**Symptom**: Cache directory doesn't exist after server starts

**Solution**:
- Check file permissions
- Ensure server has write access to parent directories
- Verify path separators (use `path.join()` in Node.js, `Path()` in Python)

### Files in Wrong Cache
**Symptom**: mol2chemfig files in MoleculeViewer cache or vice versa

**Solution**:
- Restart both servers to ensure they're using updated code
- Clear both caches
- Regenerate test images
- Verify which server is generating which files

### Permission Errors
**Symptom**: "Permission denied" when writing to cache

**Solution**:
- Check directory ownership
- Ensure cache directories are writable
- On Windows, check if folders are read-only
- Run servers with appropriate permissions

## Related Files
- `MoleculeViewer/server.js` - MoleculeViewer cache configuration
- `mol2chemfig_server.py` - mol2chemfig cache configuration
- `.gitignore` - Cache exclusions
- `CLAUDE_CONTEXT.md` - Project documentation
- `MoleculeViewer/docs/Todolist.md` - Original requirement (line 14)
