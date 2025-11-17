# 🔬 PubChem 3D Viewer Integration Guide

## ✨ What's New

Your Chrome extension now supports **interactive 3D molecular viewers** powered by PubChem! View any molecule in stunning 3D with multiple rendering styles and controls.

## 🚀 Quick Start

### 1. Start the PubChem Server

**Windows:**
```bash
cd MoleculeViewer\pubchem
start.bat
```

**macOS/Linux:**
```bash
cd MoleculeViewer/pubchem
npm start
```

Server will start on: **http://localhost:5002**

### 2. Enable 3D Viewer in Extension

1. Click the extension icon in Chrome
2. Go to **Developer Options** tab
3. Scroll to **🎨 3D Viewer Mode**
4. Enable **Enable 3D Viewer** toggle
5. Click "Save Settings"

### 3. View Molecules in 3D

On any webpage, use the chemistry syntax:

```
chem:histamine:
chem:caffeine:
chem:dopamine:
```

You'll now see a **🔮 3D** button on each molecule image. Click it to open the interactive 3D viewer!

## 🎨 3D Viewer Features

### Rendering Styles

- **Ball & Stick** - Traditional molecular model with atoms as spheres and bonds as cylinders
- **Stick** - Bonds shown as sticks only, no atom spheres
- **Space Filling** - Atoms shown at their van der Waals radii
- **Wireframe** - Minimal line representation

### Interactive Controls

- **🔄 Auto Rotate** - Automatically rotate the molecule
- **👁️ Show Hydrogens** - Toggle hydrogen atom visibility
- **🖱️ Mouse Controls:**
  - **Left Click + Drag** - Rotate molecule
  - **Scroll Wheel** - Zoom in/out
  - **Right Click + Drag** - Pan view

### Additional Options

- **🔗 Open in PubChem** - View full PubChem page with more data
- **💾 Download SDF** - Download 3D structure data file
- **✓ Close** - Close the viewer window

## 📐 How It Works

### Architecture

```
Webpage with chem:histamine:
         ↓
Extension detects formula
         ↓
Fetches 2D image from PubChem server
         ↓
Displays image with 🔮 3D button (if enabled)
         ↓
Click 3D button
         ↓
Opens 3D viewer in new window
         ↓
Loads PubChem 3D structure data
         ↓
Renders interactive 3D model
```

### Two Server Options

#### Option 1: Node.js Server (Recommended)
- Location: `MoleculeViewer/pubchem/server.js`
- Port: 5002
- Fast, lightweight, production-ready
- Direct PubChem API integration

#### Option 2: Python Flask Server
- Location: `pubchem_server.py` (root directory)
- Port: 5002
- Alternative if you prefer Python
- Same functionality

**You only need to run ONE of these servers!**

## 🎯 Usage Examples

### Example 1: Simple Molecule
```html
<!-- On your webpage -->
<p>Histamine (chem:histamine:) is a neurotransmitter.</p>
```

**Result:**
- 2D structure image appears inline
- 🔮 3D button appears in top-right corner
- Click to view in 3D

### Example 2: Multiple Molecules
```html
<p>Common neurotransmitters include:</p>
<ul>
  <li>Dopamine: chem:dopamine:</li>
  <li>Serotonin: chem:serotonin:</li>
  <li>Norepinephrine: chem:norepinephrine:</li>
</ul>
```

**Result:**
- Each molecule gets its own 2D image and 3D button
- Click any 3D button to view that molecule

### Example 3: SMILES Notation
```html
<p>Ethanol structure: chem:CCO:</p>
<p>Benzene ring: chem:c1ccccc1:</p>
```

**Result:**
- Works with SMILES notation too!

## 🔧 Configuration

### Extension Settings

**Main Tab:**
- **Renderer Engine** → Select "PubChem"
- **PubChem Image Size** → `large` (recommended), `small`, or custom (e.g., `500x500`)
- **Record Type** → `2d` or `3d` (for 3D projection images)

**Developer Options Tab:**
- **Enable 3D Viewer** → Toggle to show/hide 3D buttons
- **Enable 3D Models** → Additional setting for 3D features

### Server Configuration

Edit `MoleculeViewer/pubchem/server.js`:
```javascript
const PORT = 5002; // Change port if needed
```

## 📊 Supported Molecules

### By Name
Any compound in PubChem database:
- Common names: `histamine`, `caffeine`, `aspirin`, `glucose`
- IUPAC names: `2-aminoethylphenol`, `methyl salicylate`
- Drug names: `ibuprofen`, `acetaminophen`, `morphine`

### By SMILES
Standard SMILES notation:
- `CCO` - Ethanol
- `CC(=O)O` - Acetic acid
- `c1ccccc1` - Benzene
- `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO` - D-Glucose (with stereochemistry)

## 🧪 Testing

### Test File
Open `test_pubchem_3d.html` in your browser:

```bash
# Located in project root
# Simply open the file in Chrome
```

**Features:**
- Search any molecule by name
- View 2D structure
- Open 3D viewer
- Quick access to 6 example molecules

### Manual Testing

1. Start PubChem server
2. Open: http://localhost:5002/pubchem/3d-viewer?name=histamine&embed=true
3. Verify:
   - ✅ 3D model loads
   - ✅ Can rotate with mouse
   - ✅ Style dropdown works
   - ✅ Show hydrogens toggle works
   - ✅ Auto rotate works
   - ✅ Download SDF works

## 🎓 Tips & Tricks

### Tip 1: Better Performance
- Use `large` image size for best quality without performance hit
- PubChem images are cached automatically

### Tip 2: Stereochemistry
Enable **3D Mode** in mol2chemfig settings to preserve stereochemistry:
- D-glucose vs L-glucose will show correct 3D structure
- Chiral centers are preserved

### Tip 3: Complex Molecules
For large biomolecules:
- Use **Stick** or **Wireframe** style for clarity
- Disable **Show Hydrogens** to reduce clutter
- Use **Space Filling** to visualize molecular surface

### Tip 4: Teaching & Presentations
- Open multiple 3D viewers side-by-side to compare molecules
- Use **Auto Rotate** for presentation mode
- Download SDF files for use in other molecular viewers

## 🔍 API Reference

### Quick API Calls

```javascript
// Get 2D image
fetch('http://localhost:5002/pubchem/img/histamine')

// Get compound info
fetch('http://localhost:5002/pubchem/info?name=histamine')

// Get 3D model (SDF)
fetch('http://localhost:5002/pubchem/3d-model?name=histamine')

// Open 3D viewer programmatically
window.open('http://localhost:5002/pubchem/3d-viewer?name=histamine&embed=true', 
            '_blank', 'width=1000,height=700');
```

See `MoleculeViewer/pubchem/README.md` for complete API documentation.

## 🐛 Troubleshooting

### 3D Button Not Appearing

**Solution:**
1. Check extension popup → Developer Options
2. Verify **Enable 3D Viewer** is ON
3. Reload the webpage (Ctrl+R)
4. Check console for errors (F12)

### 3D Viewer Opens But Blank

**Solution:**
1. Check internet connection (loads data from PubChem)
2. Try a different molecule
3. Check browser console for CORS errors
4. Verify server is running on port 5002

### Molecule Not Found

**Solution:**
1. Check spelling of molecule name
2. Try SMILES notation instead
3. Search on pubchem.ncbi.nlm.nih.gov to verify compound exists
4. Try common name vs IUPAC name

### Server Won't Start

**Solution:**
```bash
# Check if Node.js is installed
node --version

# Install dependencies
cd MoleculeViewer/pubchem
npm install

# Check if port 5002 is in use
netstat -ano | findstr :5002

# Kill process on port 5002 if needed
# Then restart server
```

## 📚 Additional Resources

### PubChem Database
- Main site: https://pubchem.ncbi.nlm.nih.gov
- Search compounds: https://pubchem.ncbi.nlm.nih.gov/search
- API docs: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

### Project Documentation
- `MoleculeViewer/pubchem/README.md` - Server API documentation
- `PUBCHEM_INTEGRATION.md` - Technical implementation details
- `chem-extension/README.md` - Extension usage guide

## 🌟 Feature Comparison

| Feature | 2D Image | 3D Viewer |
|---------|----------|-----------|
| Shows structure | ✅ | ✅ |
| Interactive rotation | ❌ | ✅ |
| Multiple render styles | ❌ | ✅ (4 styles) |
| Show/hide hydrogens | ❌ | ✅ |
| Auto rotation | ❌ | ✅ |
| Stereochemistry | ⚠️ Limited | ✅ Full 3D |
| Download data | ❌ | ✅ (SDF) |
| Load time | Fast | Medium |
| Works offline | ✅ (cached) | ❌ (needs PubChem) |

## 🎉 What's Next?

Future enhancements planned:
- [ ] Offline 3D rendering (local library)
- [ ] Custom color schemes
- [ ] Measure distances and angles
- [ ] Export as PNG/SVG
- [ ] Animation of conformers
- [ ] Protein structure support

## 💡 Examples Gallery

### Common Molecules to Try

**Simple:**
- `chem:water:` or `chem:H2O:`
- `chem:ethanol:` or `chem:CCO:`
- `chem:acetone:` or `chem:CC(=O)C:`

**Aromatics:**
- `chem:benzene:`
- `chem:toluene:`
- `chem:naphthalene:`

**Biomolecules:**
- `chem:glucose:`
- `chem:alanine:`
- `chem:cholesterol:`

**Drugs:**
- `chem:aspirin:`
- `chem:ibuprofen:`
- `chem:penicillin:`

**Neurotransmitters:**
- `chem:dopamine:`
- `chem:serotonin:`
- `chem:acetylcholine:`

## 📝 Summary

You now have full 3D molecular viewing capabilities integrated into your Chrome extension! 🎉

**Key Points:**
1. ✅ PubChem server provides both 2D images and 3D viewers
2. ✅ Enable 3D Viewer in Developer Options to see 🔮 3D buttons
3. ✅ Click 3D button to open interactive 3D molecular viewer
4. ✅ Works with both molecule names and SMILES notation
5. ✅ Multiple rendering styles and interactive controls

**Start exploring molecules in 3D today!** 🧪🔬

---

Made with ❤️ for chemistry education and research
