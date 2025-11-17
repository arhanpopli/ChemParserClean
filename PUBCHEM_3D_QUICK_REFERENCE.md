# 🚀 PubChem 3D Viewer - Quick Reference

## ⚡ Quick Start (30 seconds)

```bash
# 1. Start server
cd MoleculeViewer\pubchem
start.bat

# 2. Enable in extension
#    Popup → Developer Options → Enable 3D Viewer ✓

# 3. Use on webpage
#    Type: chem:histamine:
#    Click: 🔮 3D button
```

## 📍 Key Files

| File | Purpose | Location |
|------|---------|----------|
| `server.js` | PubChem Node.js server | `MoleculeViewer/pubchem/` |
| `start.bat` | Start server (Windows) | `MoleculeViewer/pubchem/` |
| `content.js` | Extension integration | `chem-extension/` |
| `test_pubchem_3d.html` | Test interface | Project root |

## 🎨 3D Viewer Controls

| Control | Action |
|---------|--------|
| 🎨 Style | Ball & Stick, Stick, Space Filling, Wireframe |
| 👁️ Hydrogens | Show/Hide hydrogen atoms |
| 🔄 Rotate | Auto-rotate molecule |
| 🖱️ Mouse | Drag=rotate, Scroll=zoom |
| 🔗 PubChem | Open full page |
| 💾 Download | Get SDF file |

## 🧪 Example Molecules

```
chem:histamine:
chem:caffeine:
chem:dopamine:
chem:aspirin:
chem:glucose:
chem:ethanol:
```

## 🔗 API Endpoints

```
# 2D Image
http://localhost:5002/pubchem/img/histamine

# 3D Viewer
http://localhost:5002/pubchem/3d-viewer?name=histamine&embed=true

# Compound Info
http://localhost:5002/pubchem/info?name=histamine

# 3D Model (SDF)
http://localhost:5002/pubchem/3d-model?name=histamine
```

## ⚙️ Extension Settings

**Main Tab:**
- Renderer Engine: **PubChem**
- Image Size: **large**
- Record Type: **2d**

**Developer Options:**
- Enable 3D Viewer: **ON** ✓

## 🐛 Quick Fixes

| Problem | Solution |
|---------|----------|
| No 3D button | Enable in Developer Options |
| Button doesn't work | Check server running (port 5002) |
| Molecule not found | Check spelling or try SMILES |
| Server won't start | Run `npm install` in pubchem folder |

## 📚 Documentation

- **User Guide**: `PUBCHEM_3D_GUIDE.md`
- **API Docs**: `MoleculeViewer/pubchem/README.md`
- **Implementation**: `PUBCHEM_3D_IMPLEMENTATION_COMPLETE.md`

## ✅ Checklist

- [ ] Server started on port 5002
- [ ] Extension installed in Chrome
- [ ] 3D Viewer enabled in settings
- [ ] Test page opened
- [ ] 3D button appears
- [ ] 3D viewer opens
- [ ] Controls work

## 🎯 Status

**Implementation:** ✅ COMPLETE  
**Testing:** ⏳ YOUR TURN  
**Server:** ✅ Running on http://localhost:5002

---

**Everything is ready! Just start the server and test it!** 🚀
