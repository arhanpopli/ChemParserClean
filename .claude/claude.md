# ChemParser Project Context

## Project Overview
ChemParser is a chemical structure visualization system with Chrome extension integration. It replaces `chem:chemicalname:` patterns in web pages with rendered molecular structure images.

## Architecture

### Servers (4 total)
| Server | Port | Technology | Purpose |
|--------|------|------------|---------|
| MoleculeViewer | 5000 | Node.js | SMILES/nomenclature to SVG conversion |
| Mol2ChemFig Flask | 5001 | Python/Flask | Wrapper for Docker backend |
| Mol2ChemFig Docker | 8000 | Docker | Actual mol2chemfig rendering |
| PubChem | 5002 | Python/Flask | PubChem images and 3D models |

### Key Files
- `MoleculeViewer/server.js` - Node.js server (port 5000)
- `mol2chemfig_server.py` - Flask wrapper (port 5001)
- `pubchem_server.py` - PubChem server (port 5002)
- `docker-compose.yml` - Docker config for mol2chemfig backend
- `chem-extension/content.js` - Chrome extension content script
- `chem-extension/popup.js` - Extension settings popup

### Web Interfaces
- `http://localhost:5000/unified-interface.html` - Main unified interface with all tabs
- `http://localhost:5002/static/viewer-3d.html?name=X&embed=true` - 3D viewer page

## Starting the Servers

```bash
# 1. Start Docker Desktop first (required for mol2chemfig)

# 2. Start Docker backend
docker-compose up -d

# 3. Start all servers (in separate terminals or use launcher)
node MoleculeViewer/server.js    # Port 5000
python mol2chemfig_server.py     # Port 5001
python pubchem_server.py         # Port 5002
```

## Extension Configuration

### API Endpoints (in content.js)
```javascript
const MOLECULE_VIEWER_API = 'http://localhost:5000';
const MOL2CHEMFIG_API = 'http://localhost:5001';  // Flask wrapper, NOT 8000
const PUBCHEM_API = 'http://localhost:5002';
```

### Rendering Engines (user selectable)
1. **MoleculeViewer** (default) - Uses `/img/smiles` or `/img/nomenclature`
2. **mol2chemfig** - Uses `/m2cf/submit` on port 5001
3. **PubChem** - Uses `/pubchem/img/{compound}`

### 3D Viewer Feature
- Already implemented in `content.js` function `show3DViewerInline()`
- Uses MolView.org iframes: `https://embed.molview.org/v1/?mode=balls&smiles=XXX`
- Enable via popup toggle: "Enable 3D Viewer"
- Works with PubChem renderer engine

## Cache Directories
- `MoleculeViewer/cache/moleculeviewer/` - MoleculeViewer cache
- `cache/mol2chemfig/` - Mol2ChemFig cache
- `pubchem-cache/` - PubChem cache

## Common Issues

### "Mol2ChemFig Load Failed" Error
**Cause:** Docker Desktop not running or container not started
**Fix:**
1. Start Docker Desktop
2. Run `docker-compose up -d`
3. Verify: `curl http://localhost:8000/m2cf/reset -X POST`

### "Error: Load Failed" (MoleculeViewer)
**Cause:** Server not running on port 5000
**Fix:** Start `node MoleculeViewer/server.js`

### "PubChem: Not Found"
**Cause:** Server not running on port 5002 OR compound doesn't exist
**Fix:** Start `python pubchem_server.py`

## Health Checks
```bash
# Check all servers
curl http://localhost:5000/health  # MoleculeViewer
curl http://localhost:5001/health  # mol2chemfig Flask
curl http://localhost:8000/m2cf/reset -X POST  # Docker backend
curl http://localhost:5002/health  # PubChem
```

## Test Files
- `test_m2cf_full.html` - mol2chemfig testing
- `test_pubchem.html` - PubChem testing
- `test_mol2chemfig_fixes.html` - Feature testing

## Extension Usage
Replace chemical names in text with: `chem:ethanol:` or `chem:CCO:` (SMILES)

The extension will:
1. Detect `chem:X:` patterns
2. Call selected rendering engine
3. Replace text with molecular structure image
4. Optionally show 3D viewer (if enabled)
