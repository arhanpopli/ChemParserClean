# Cache System Quick Reference

## New Cache Structure

```
Chemparser/
├── MoleculeViewer/
│   └── cache/
│       └── moleculeviewer/      ← RDKit SVGs (port 5000)
│
└── cache/
    └── mol2chemfig/              ← mol2chemfig SVGs/PDFs (port 5001)
```

## File Patterns

### MoleculeViewer Cache
- **Location**: `MoleculeViewer/cache/moleculeviewer/`
- **File Format**: `{md5_hash}.svg`
- **Example**: `87970ea81928442541e33f0907da38bc.svg`
- **Hash**: MD5 of `type:value:widthxheight`

### mol2chemfig Cache
- **Location**: `cache/mol2chemfig/`
- **File Format**: `{sha256_hash}.{ext}`
- **Example**: `a1b2c3d4e5f6g7h8.svg`, `a1b2c3d4e5f6g7h8.pdf`
- **Hash**: SHA256 of `smiles:options`

## Quick Commands

### Start Servers
```bash
# MoleculeViewer
cd MoleculeViewer && node server.js

# mol2chemfig wrapper
python mol2chemfig_server.py

# Docker backend
docker-compose up -d
```

### Test Cache Separation
```bash
python test_cache_separation.py
```

### View Cache Stats
```bash
# MoleculeViewer
curl http://localhost:5000/cache-info

# mol2chemfig
curl http://localhost:5001/api/cache/stats
```

### Clear Cache
```bash
# MoleculeViewer
curl -X DELETE http://localhost:5000/clear-cache

# mol2chemfig
curl -X POST http://localhost:5001/api/cache/clear
```

## Migration from Old Caches

### Backup Old Caches (Optional)
```bash
# Windows
move MoleculeViewer\svg-cache MoleculeViewer\svg-cache.backup
move mol2chemfig_storage mol2chemfig_storage.backup

# Linux/Mac
mv MoleculeViewer/svg-cache MoleculeViewer/svg-cache.backup
mv mol2chemfig_storage mol2chemfig_storage.backup
```

### Delete Old Caches (After Testing)
```bash
# Windows
rmdir /s /q MoleculeViewer\svg-cache.backup
rmdir /s /q mol2chemfig_storage.backup

# Linux/Mac
rm -rf MoleculeViewer/svg-cache.backup
rm -rf mol2chemfig_storage.backup
```

## Troubleshooting

### Cache Not Created
1. Check server logs for errors
2. Verify file permissions
3. Ensure servers are using updated code (restart if needed)

### Files in Wrong Location
1. Restart both servers
2. Clear all caches
3. Run test script: `python test_cache_separation.py`

### Permission Errors
1. Check directory ownership
2. Verify write permissions on parent directories
3. On Windows, ensure folders are not read-only

## Related Documentation
- Full details: `CACHE_SEPARATION_IMPLEMENTATION.md`
- Project context: `CLAUDE_CONTEXT.md`
- Original requirement: `MoleculeViewer/docs/Todolist.md` (line 14)
