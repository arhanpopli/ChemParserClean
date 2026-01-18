# ChemParser Project Context

## Project Overview
ChemParser is a fully client-side Chrome extension that replaces `chem:chemicalname:` patterns in web pages with rendered molecular structure images. **No local servers required!**

## Architecture (v4.0 - Fully Client-Side)

The extension works completely offline using:
- **IntegratedSearch** - Queries PubChem, RCSB, COD APIs directly (no local server)
- **SmilesDrawer** - Client-side 2D molecular rendering
- **3Dmol.js** - Client-side 3D molecular viewing
- **CDK Depict** - Optional external API for 2D rendering

### Key Files (chem-extension/)
| File | Purpose |
|------|---------|
| `manifest.json` | Extension configuration |
| `content.js` | Main content script - pattern detection & rendering |
| `background.js` | Service worker |
| `popup.html/js` | Settings popup |
| `integrated-search.js` | Client-side search (PubChem, RCSB, COD APIs) |
| `smiles-drawer.min.js` | 2D structure rendering |
| `3Dmol-min.js` | 3D viewer library |
| `3dmol-viewer.html/js` | 3D viewer page |
| `2d-viewer.html/js` | 2D viewer page |

## Extension Usage

Replace chemical names in text with: `chem:ethanol:` or `chem:CCO:` (SMILES)

The extension will:
1. Detect `chem:X:` patterns
2. Query IntegratedSearch (PubChem/RCSB/COD)
3. Render structure with SmilesDrawer (client-side)
4. Optionally show 3D viewer with 3Dmol.js

## Rendering Engines

| Engine | Type | Description |
|--------|------|-------------|
| SmilesDrawer | Client-side | Fast, works offline (default) |
| CDK Depict | External API | Advanced options, requires internet |

## Data Sources (via IntegratedSearch)

| Source | Data Type | API |
|--------|-----------|-----|
| PubChem | Compounds (100M+) | pubchem.ncbi.nlm.nih.gov |
| RCSB PDB | Proteins/Biomolecules (250k+) | rcsb.org |
| COD | Minerals/Crystals (500k+) | crystallography.net |

## Development

### Loading the Extension
1. Open Chrome → `chrome://extensions`
2. Enable "Developer mode"
3. Click "Load unpacked"
4. Select the `chem-extension` folder

### Testing
- Open any webpage
- Type `chem:aspirin:` in a text field
- The extension will render the structure

### Settings
- Click the extension icon to open settings
- Choose rendering engine (SmilesDrawer or CDK Depict)
- Configure rendering options
- Toggle 3D viewer

## Project Structure

```
Chemparser/
├── .claude/           # Claude AI configuration
├── .github/           # GitHub workflows
├── chem-extension/    # Chrome extension (fully self-contained)
│   ├── manifest.json
│   ├── content.js
│   ├── background.js
│   ├── popup.html/js
│   ├── integrated-search.js
│   ├── smiles-drawer.min.js
│   ├── 3Dmol-min.js
│   └── ...
└── README.md
```

## No Servers Required!

Previous versions required multiple local servers. Version 4.0+ is fully client-side:

- ~~MoleculeViewer (port 5000)~~ → SmilesDrawer (client-side)
- ~~mol2chemfig (port 5001)~~ → Removed
- ~~PubChem server (port 5002)~~ → IntegratedSearch (direct API)
- ~~MolView server (port 8000/8001)~~ → embed.molview.org (external)
