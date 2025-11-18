# Chemparser Project - Context for AI Agents

## Project Overview
Chemparser is a chemistry visualization project with four main components:
1. **MoleculeViewer** - Node.js server (port 5000) using RDKit to generate SVG images from SMILES
2. **mol2chemfig** - Docker backend (port 8000) converting molecules to ChemFig LaTeX format with SVG output
3. **PubChem** - Python Flask server (port 5002) fetching images and 3D models from PubChem API
4. **chem-extension** - Chrome extension that detects `chem:` syntax and displays molecule images

## Architecture

### Servers
- **MoleculeViewer** (port 5000):
  - Node.js/Express server
  - Path: `MoleculeViewer/server.js`
  - Python backend: `MoleculeViewer/generate_svg.py`
  - Cache: `MoleculeViewer/cache/moleculeviewer/`

- **mol2chemfig Docker Backend** (port 8000):
  - Docker container from `docker-compose.yml`
  - Main file mounted: `m2cf_fixed.py` → `/usr/src/app/src/m2cf.py`
  - API endpoints: `/m2cf/submit`, `/m2cf/search`, `/m2cf/layers`
  - Recently updated: atom separation from 16pt to 28pt, font size from 8pt to 12pt

- **mol2chemfig Wrapper** (port 5001):
  - Optional Flask wrapper: `mol2chemfig_server.py`
  - Cache: `cache/mol2chemfig/`

### Extension
- **Path**: `chem-extension/`
- **Key Files**:
  - `manifest.json` - Extension config
  - `content.js` - Main logic (2500+ lines)
  - `popup.html` + `popup.js` - Settings UI
  - `styles.css` - Styling

- **APIs Used**:
  ```javascript
  const MOLECULE_VIEWER_API = 'http://localhost:5000';
  const MOL2CHEMFIG_API = 'http://localhost:8000';
  ```

## Recent Changes
1. ✅ Test infrastructure created (`test_runner.py`)
2. ✅ mol2chemfig SVG size increased (28pt atom sep, 12pt font)
3. 🔄 Dark mode for mol2chemfig (in progress)

## Key File Locations
- Extension content script: `chem-extension/content.js`
- mol2chemfig backend: `m2cf_fixed.py` (lines 140-151 for LaTeX template)
- MoleculeViewer server: `MoleculeViewer/server.js`
- Docker config: `docker-compose.yml`
- Environment vars: `.env`

## Testing
- Test runner: `test_runner.py`
- Start servers:
  - MoleculeViewer: `cd MoleculeViewer && node server.js`
  - mol2chemfig: `docker-compose up -d`
- Test endpoints:
  - MoleculeViewer: `http://localhost:5000/img/smiles?smiles=CCO`
  - mol2chemfig: `POST http://localhost:8000/m2cf/submit` with `{"textAreaData":"CCO"}`

## Common Patterns

### Dark Mode in Extension
```javascript
const isDarkMode = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
if (isDarkMode) {
  svgContent = svgContent
    .replace(/#000000/gi, '#FFFFFF')
    .replace(/#000\b/gi, '#FFF')
    .replace(/stroke="black"/gi, 'stroke="white"')
    .replace(/fill="black"/gi, 'fill="white"');
}
```

### LaTeX Template (m2cf_fixed.py)
```latex
\documentclass{minimal}
\usepackage{xcolor, mol2chemfig}
\setatomsep{28pt}  % Size control
\renewcommand{\printatom}[1]{\fontsize{12pt}{14pt}\selectfont{\ensuremath{\mathsf{#1}}}}
\begin{document}
%s
\end{document}
```

## Important Notes
- Docker containers must be restarted after changing `m2cf_fixed.py`
- Extension needs reload in `chrome://extensions` after changes
- SVG colors use mix of formats: `#000000`, `stroke='#000000'`, `stroke="black"`
- Cache folders (separated):
  - MoleculeViewer: `MoleculeViewer/cache/moleculeviewer/`
  - mol2chemfig wrapper: `cache/mol2chemfig/`
  - Docker backend: `backend_data/` (internal use only)

## Todos Reference
See `MoleculeViewer/docs/Todolist.md` for detailed requirements.
