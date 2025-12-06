# ChemTex

**ChemTex** is a powerful Chrome extension for rendering chemical structures and formulas directly in your browser. It supports mhchem, chemfig, SMILES notation, and molecule names with intelligent search across PubChem, RCSB PDB, and the Crystallography Open Database (COD).

## ✨ Features

- **Client-Side Rendering**: All structures render in your browser using SmilesDrawer - no server required!
- **Intelligent Search**: Automatically detects and renders:
  - **Compounds** from PubChem (100M+ molecules)
  - **Biomolecules/Proteins** from RCSB PDB (250k+ structures)
  - **Minerals** from COD (500k+ crystals)
- **Multiple Input Formats**:
  - `chem:caffeine:` - Molecule names
  - `chem:CCO:` - SMILES notation
  - `\ce{H2O}` - mhchem formulas
  - `\chemfig{...}` - chemfig structures
- **3D Visualization**: Interactive 3D molecular viewers with customizable styles
- **Customizable Rendering**: Control carbons, hydrogens, aromatic rings, atom numbering, rotations, and more
- **White Background Removal**: Clean protein structure images from RCSB

## 🚀 Installation

### From Source
1. Clone this repository:
   ```bash
   git clone https://github.com/arhanpopli/ChemTex.git
   cd ChemTex
   ```

2. Load the extension in Chrome:
   - Open `chrome://extensions/`
   - Enable "Developer mode"
   - Click "Load unpacked"
   - Select the `chem-extension` folder

### From Chrome Web Store
*Coming soon!*

## 📖 Usage

### Basic Syntax
Use the `chem:` prefix followed by a molecule name or SMILES:

```
chem:aspirin:
chem:rhinovirus:
chem:calcite:
chem:CCO:
```

### Advanced Options
Add flags to customize rendering:

```
chem:benzene/d+c+h+o+r45:
```

**Flags:**
- `+c` - Show carbon atoms
- `+h` - Show hydrogen atoms
- `+o` - Aromatic ring circles
- `+m` - Show methyl groups (CH₃)
- `+n` - Atom numbering
- `+p` - Flip horizontal
- `+q` - Flip vertical
- `+r45` - Rotate 45°
- `+s150` - Size 150%
- `+3d` - Show 3D viewer

## 🔧 Configuration

Click the extension icon to access settings:

- **Rendering Engine**: SmilesDrawer (recommended) or PubChem
- **Rendering Options**: Customize appearance (carbons, hydrogens, themes, etc.)
- **3D Viewer Settings**: Choose viewer source, style, size, and colors
- **Protein Options**: White background removal for RCSB images
- **Mineral Options**: Crystallography display settings

## 🏗️ Architecture

### Core Components

- **IntegratedSearch** (`integrated-search.js`): Client-side search module that queries:
  - RCSB PDB Search API for biomolecules
  - COD API for minerals
  - PubChem PUG REST API for compounds
  
- **SmilesDrawer**: Client-side 2D structure renderer
- **3Dmol.js**: Interactive 3D molecular visualization
- **MolView Integration**: Optional local MolView server support

### No Server Required!
ChemTex v6.0+ uses direct API calls to external databases, eliminating the need for a local search server.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 💖 Support

If you find ChemTex useful, consider supporting development:

- **Discord**: [Join our community](https://discord.gg/YOUR_INVITE)
- **Donate**: ETH donations accepted via MetaMask (click the extension icon)

## 📄 License

MIT License - see LICENSE file for details

## 🙏 Acknowledgments

- [SmilesDrawer](https://github.com/reymond-group/smilesDrawer) - 2D structure rendering
- [3Dmol.js](https://3dmol.csb.pitt.edu/) - 3D visualization
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Compound database
- [RCSB PDB](https://www.rcsb.org/) - Protein structures
- [COD](https://www.crystallography.net/) - Mineral structures

## 📝 Changelog

### v6.0 (2025-01-06)
- 🎉 Renamed to ChemTex
- ✨ Added IntegratedSearch module (no server needed!)
- 🔍 Direct API queries to RCSB, COD, and PubChem
- 💎 Enhanced mineral and biomolecule support
- 🎨 Added Discord and MetaMask donate buttons
- 🧹 Cleaned up UI (removed verbose info boxes)
- 🐛 Fixed RCSB API query format
- ⚡ Increased API timeouts for reliability
- 🎯 Improved similarity scoring for better search results

---

**Made with ❤️ by [Arhan Popli](https://github.com/arhanpopli)**
