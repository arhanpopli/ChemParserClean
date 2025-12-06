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

---

## MolView Local Instance (Port 8000)

### Directory Structure
```
Molview/molview/
├── embed/
│   ├── v1/          # Lightweight embed (4KB HTML) - limited fallback
│   └── v2/          # Full app embed (50KB HTML) - full fallback chain
├── src/js/
│   ├── Loader.js    # Full loader with 3D fallback
│   ├── Loader.Embed.js  # Embed loader (updated with fallback)
│   ├── Request.js   # API calls to PubChem, RCSB, COD, CIR
│   ├── Model.js     # 3D render engine wrapper
│   └── chem/        # Chemistry libs (MolFile parser, SMILES gen)
├── search.php       # PHP search API
└── search-server.js # Node.js search server (port 8001)
```

### Embed URLs
- V1: `http://localhost:8000/embed/v1/?cid=XXX` - lightweight, limited fallback
- V2: `http://localhost:8000/embed/v2/?cid=XXX` - full features, full fallback

### Search Server (Port 8001)
Unified search across PubChem, RCSB, COD with compact JSON responses:
```
http://localhost:8001/search?q=rhinovirus
http://localhost:8001/search?q=calcite
http://localhost:8001/search?q=aspirin
```

Returns: name, SMILES, SDF availability, embed URLs, source type (compound/biomolecule/mineral)

---

## 3D Fallback Mechanism (CRITICAL)

### Why Some Complex Molecules Fail in Embed V1

PubChem doesn't have 3D coordinates for all molecules (e.g., cardiolipin CID:166177218).

### Fallback Chain (Main App & V2)
```
1. PubChem 3D SDF → if fails:
2. Get SMILES from PubChem → send to CIR → if fails:
3. Parse 2D SDF → extract SMILES → send to CIR → if fails:
4. Load 2D as 3D (flat molecule)
```

### CIR (Chemical Identifier Resolver)
- Service: `https://cactus.nci.nih.gov/chemical/structure/`
- Generates 3D coordinates from SMILES
- Example: `https://cactus.nci.nih.gov/chemical/structure/CCO/file?format=sdf&get3d=True`

### Documentation
See: `MOLVIEW_FALLBACK_RESEARCH.md` for detailed analysis

---

## Future: 3Dmol.js Integration

Goal: Replace MolView iframe dependency with direct 3Dmol.js rendering

### Required Implementation
1. Fetch structure data from PubChem/RCSB/COD APIs
2. Implement CIR fallback for missing 3D coordinates
3. Use 3Dmol.js library for WebGL rendering
4. Handle all molecule types: compounds, proteins, minerals

### Data Sources
| Type | API | Format |
|------|-----|--------|
| Compounds | PubChem | SDF |
| Proteins | RCSB | PDB |
| Minerals | COD | CIF |
| 3D Generation | CIR | SDF from SMILES |

---

## Session Log

### 2024-12-05: MolView Fallback Investigation
- Analyzed why embed/v1 shows FLAT molecules while v2 shows proper 3D bends
- **ROOT CAUSE FOUND:**
  - V1 has CIR code but no Sketcher component to extract SMILES from 2D
  - V2 uses `Sketcher.getSMILES()` to get SMILES, then CIR generates 3D
  - **The `get3d=True` parameter in CIR API generates 3D coordinates via force field optimization**
- V1's try/catch for MolFile silently fails because MolFile class isn't in the bundle
- **Solution for 3Dmol.js:** Get SMILES from PubChem API, send to CIR with `get3d=True`
- Documented findings in `MOLVIEW_FALLBACK_RESEARCH.md`

### Key API for 3D Generation
```
https://cactus.nci.nih.gov/chemical/structure/{SMILES}/file?format=sdf&get3d=True
```
This endpoint takes a SMILES string and returns an SDF with computed 3D coordinates!

### Implementation Complete
- Added `getSMILESFromPubChem(cid)` - fetches SMILES from PubChem
- Added `generate3DFromCIR(smiles)` - generates 3D via CIR with `get3d=True`
- Added `get2DFromPubChem(cid)` - fallback for flat 2D
- Updated `loadMolecule()` with full fallback chain:
  1. PubChem 3D conformer
  2. PubChem computed 3D
  3. **CIR 3D generation from SMILES** (the key for cardiolipin!)
  4. 2D as 3D (flat fallback)
  5. MolView iframe (emergency)

### Test URLs
- Cardiolipin: `3dmol-viewer.html?cid=166177218` (uses CIR fallback for 3D)
- Aspirin: `3dmol-viewer.html?cid=2244` (has PubChem 3D)
- Insulin: `3dmol-viewer.html?pdbid=4ins` (direct RCSB PDB loading)
- Rhinovirus: `3dmol-viewer.html?pdbid=4rhv` (large biomolecule, cartoon style)
- By name: `3dmol-viewer.html?name=hemoglobin` (auto-detects type)

### Update 2: Direct RCSB/Biomolecule Support
- Added `getPDBFromRCSB(pdbId)` - fetches PDB directly from RCSB
- Added `loadBiomolecule(pdbId)` - renders proteins with 3Dmol.js
- Added URL parameter `?pdbid=XXXX` for direct PDB loading
- Added RCSB search fallback for protein names
- Auto-styling: cartoon for large (>500 atoms), ball-stick for small
- **No MolView needed for biomolecules anymore!**

### Improved CIR Validation
- Validates SDF has V2000/V3000 marker (not HTML error page)
- Checks minimum line count
- Properly URL-encodes SMILES with `encodeURIComponent()`
