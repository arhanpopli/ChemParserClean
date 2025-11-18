# Mol2chemfig - Complete Molecule Visualization System

## 🚀 Quick Start (30 seconds)

### Windows Users
1. Go to: `C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer`
2. **Double-click `START_ALL.bat`**
3. Your browser opens automatically
4. That's it! The system is running

## 📚 What You Have

### MoleculeViewer Backend
- Flask server that converts SMILES/chemical names to SVG molecules
- Runs on: **http://localhost:5000**
- Web interface for direct testing

### Chrome Extension
- Works in ChatGPT and other sites
- Renders molecules inline from `chem:text:` tags

## 💬 How to Use

### In ChatGPT with Extension
Type in ChatGPT:
```
Show me chem:benzene:

Or: chem:1-chloro-benzene:

Or SMILES: chem:CC(=O)C:
```

The extension will replace the `chem:text:` with a rendered molecule image.

### Direct Web Interface
1. Go to: http://localhost:5000
2. Enter a SMILES string or chemical name
3. Click "Convert to SVG"
4. Your molecule appears instantly

## 🛠️ File Structure

```
MoleculeViewer/
├── START_ALL.bat           ← Use this to start everything
├── START.bat               ← Just starts Flask server
├── app/
│   └── api.py              ← Flask backend
├── templates/
│   └── index.html          ← Web interface
└── QUICKSTART.md           ← Quick reference

chem-extension/
├── manifest.json           ← Extension config
├── content.js              ← Main extension code
└── USAGE.md                ← Extension guide
```

## 🔧 Advanced Usage

### Start Just the Server
```bash
cd MoleculeViewer
python -m flask --app app.api run --host=0.0.0.0 --port=5000
```

Then go to: http://localhost:5000

### Use on Different Port
Edit `START_ALL.bat` and change `5000` to your desired port number.

## 📋 Requirements

- Python 3.8+
- RDKit (for molecule structure)
- Flask (web server)
- Chrome browser with the extension installed

## ✨ Features

### Molecule Rendering
- SMILES to SVG conversion
- Chemical name lookup
- Support for complex nomenclature (1-chloro-benzene, etc.)
- Customizable visualization options

### Web Interface
- Dark mode support
- Visualization options (show carbons, methyls, etc.)
- Aromatic ring detection
- 24-hour caching

### Chrome Extension
- Inline rendering in ChatGPT
- Automatic SMILES/name detection
- Double-colon syntax for clarity: `chem:name:`

## 🐛 Troubleshooting

| Problem | Solution |
|---------|----------|
| "Port 5000 already in use" | Close other Flask processes: `taskkill /F /IM python.exe` |
| "Failed to fetch" | Restart with START_ALL.bat |
| Extension not showing molecules | Reload ChatGPT page (F5) |
| Complex names not working | Use double colons: `chem:1-chloro-benzene:` |

## 📝 Examples

### Working Syntax (ALWAYS use colons on both sides)
```
chem:benzene:           ✅
chem:CCO:              ✅
chem:aspirin:          ✅
chem:1-methyl-heptane: ✅
```

### Incorrect Syntax (won't work)
```
chem:benzene           ❌ Missing second colon
chem:1-chloro-benzene  ❌ Missing second colon
```

## 🎯 Next Steps

1. **First time setup:**
   - Double-click `START_ALL.bat`
   - Your browser opens automatically

2. **Test in MoleculeViewer web interface:**
   - Try: `CCO`, `c1ccccc1`, `benzene`, `aspirin`

3. **Test in ChatGPT with extension:**
   - Type: `Look up chem:benzene:`
   - Extension renders it as an image

## 💡 Tips

- Use double colons for **all** molecule tags: `chem:anything:`
- Names can have numbers, hyphens, and spaces inside the colons
- SMILES strings should use valid SMILES notation
- Molecule images are cached for 24 hours
- Each molecule gets a unique cache file

---

**That's all! Questions? See QUICKSTART.md in the MoleculeViewer folder.**
