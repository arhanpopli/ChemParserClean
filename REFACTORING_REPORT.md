# ChemParser Refactoring Report

## Summary of Cleanup Completed

### Files/Folders Removed

#### Server Folders (Completely Removed)
- `mol2chemfig/` - Python mol2chemfig library (no longer used)
- `ChemDoodleWeb-11.0.0/` - Unused ChemDoodle library
- `backend_source_mol2chemfig/` - Old backend source
- `backend_source_chemistry/` - Old backend source
- `chem-extension-unstable/` - Unstable extension copy
- `chemistry/` - Old chemistry folder
- `frontend_source/` - Old frontend source
- `frontend_static/` - Old frontend static
- `backend_data/` - Old backend data
- `docs/` - 100+ old documentation files
- `cache/` - Old cache folders
- `temp/`, `tests/`, `__pycache__/` - Temporary folders

#### Root Files Removed (~50+ files)
- Old batch files: `dev-start-*.bat`, `util-*.bat`, `1-start-all.bat`
- Old Python servers: `mol2chemfig_server.py`, `pubchem_server.py`, `unified_server.py`
- Old configs: `.env`, `docker-compose.yml`, `mcp_config.json`
- Old interfaces: `unified-interface.html`, `quintessen*.css/js`
- Test files: `test_*.py`, `test_*.html`, `test-*.html`
- Documentation: 30+ markdown files

#### Extension Cleanup (chem-extension/)
- 13 old markdown documentation files
- 4 old popup backup files
- 7 test HTML files
- `designs/` and `popup-designs/` folders (26MB of mockups)
- `csp-bypass.js`, `LOAD-ME.txt`, `popup-animation.js`

### Current Project Structure

```
Chemparser/
├── .claude/           # Claude AI configuration
├── .git/              # Git repository
├── .github/           # GitHub workflows
├── chem-extension/    # Chrome extension (20 essential files)
│   ├── manifest.json
│   ├── content.js     # Main content script
│   ├── background.js  # Service worker
│   ├── popup.html/js  # Settings popup
│   ├── integrated-search.js  # Client-side search (no server needed)
│   ├── smiles-drawer.min.js  # Client-side 2D rendering
│   ├── 3Dmol-min.js   # 3D viewer library
│   ├── kekule.min.js  # Chemistry library
│   └── ...
├── Molview/           # MolView server (for search API)
│   └── molview/
│       └── search-server.js  # Port 8001 search API
├── smilesDrawer/      # SmilesDrawer library source
├── MoleculeViewer/    # PARTIALLY REMOVED (pubchem subfolder locked)
├── README.md
├── package.json
└── start-*.bat        # Server startup scripts
```

---

## Refactoring Recommendations

### 1. Remove MoleculeViewer Code from content.js

**Current State:** content.js still has ~50 references to MoleculeViewer API (localhost:5000)

**What to Do:**
- The `loadMoleculeViewerImage()` function and related code can be removed
- Keep only: SmilesDrawer (client-side), CDK Depict, and RCSB/PubChem image fetching
- Remove `MOLECULE_VIEWER_API` constant and fallback code

**Lines to Remove/Refactor:**
- Line 11: `const MOLECULE_VIEWER_API = 'http://localhost:5000';`
- Lines 2959-3334: `loadMoleculeViewerImage()` function (can be simplified)
- Lines 3217-3268: MoleculeViewer API calls
- Lines 5640-5680: `buildMoleculeViewerRequest()` function

**Recommendation:** Keep the structure but change MoleculeViewer to use CDK Depict as fallback:
```javascript
// Replace MOLECULE_VIEWER_API with CDK Depict
const CDK_DEPICT_API = 'https://www.simolecule.com/cdkdepict/depict/bow/svg';
```

### 2. Simplify Rendering Engine Options

**Current State:** popup.js offers multiple engines:
- SmilesDrawer (client-side) ✅ Keep
- CDK Depict (external API) ✅ Keep
- MoleculeViewer (dead server) ❌ Remove
- mol2chemfig (dead server) ❌ Remove
- PubChem images ✅ Keep

**What to Do:**
- Remove MoleculeViewer and mol2chemfig options from popup.html/js
- Default to SmilesDrawer (works offline, no server needed)
- CDK Depict as fallback (external API, reliable)

### 3. Remove Locked MoleculeViewer Folder

**Current State:** `MoleculeViewer/pubchem` subfolder is locked

**What to Do:**
1. Close any file explorers/terminals accessing that folder
2. Run: `rm -rf MoleculeViewer`
3. Or add to `.gitignore` if it can't be deleted

### 4. Consolidate Server Scripts

**Current State:**
- `start-molview.bat` - MolView PHP server
- `start-search-server.bat` - Search API server
- `run-search-server.js` - Node wrapper
- `launcher-server.js` - Old launcher

**Recommendation:**
- Keep only `start-search-server.bat` (essential for search API)
- Remove `start-molview.bat` if PHP server not used
- Remove `launcher-server.js` if not needed

### 5. Update CLAUDE.md

**Current State:** `.claude/CLAUDE.md` references dead servers (ports 5000, 5001, 5002)

**What to Do:**
- Remove MoleculeViewer server docs (port 5000)
- Remove mol2chemfig server docs (port 5001)
- Remove PubChem server docs (port 5002)
- Keep only MolView/search-server.js docs (port 8001)

---

## RCSB 2D Image Status

**Status: Working Correctly**

The RCSB 2D image fetching is already properly implemented:

1. `search-server.js` returns `image_url` from RCSB CDN:
   ```javascript
   result.image_url = `https://cdn.rcsb.org/images/structures/${pdbId.toLowerCase()}_model-1.jpeg`;
   ```

2. `content.js` displays the image with optional white background removal:
   - Line 3036: `previewImg.src = moleculeData.imageUrl;`
   - Lines 3054-3115: Canvas-based background removal

No changes needed for RCSB images.

---

## Priority Order for Future Work

1. **HIGH:** Remove MoleculeViewer references from content.js
2. **HIGH:** Update popup.html/js to remove dead engine options
3. **MEDIUM:** Delete locked MoleculeViewer folder
4. **MEDIUM:** Update .claude/CLAUDE.md
5. **LOW:** Consolidate/remove unused server scripts

---

## Estimated Cleanup Impact

| Metric | Before | After |
|--------|--------|-------|
| Total Files | ~7000+ | ~500 |
| Repository Size | ~150MB | ~15MB |
| Active Servers Needed | 4 | 1 (search-server) |
| Extension Files | 50+ | 20 |
| Documentation Files | 150+ | 2 (README, USAGE) |

The project is now significantly cleaner and more maintainable.
