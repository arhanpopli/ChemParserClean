# 🚀 MoleculeViewer Quick Reference

## ⚡ Quick Start (30 seconds)

```powershell
# 1. Install dependencies (one time only)
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer
npm install
pip install rdkit requests

# 2. Start server (keep running)
npm start

# 3. Test in browser
http://localhost:5000/img/smiles?smiles=CCO

# 4. Reload Chrome extension & test in ChatGPT
# Type: chem:acetone
```

---

## 🔗 Image URLs (Copy & Paste)

### SMILES Examples
```
Ethanol:
http://localhost:5000/img/smiles?smiles=CCO

Benzene:
http://localhost:5000/img/smiles?smiles=c1ccccc1

Acetone:
http://localhost:5000/img/smiles?smiles=CC(=O)C

Methane:
http://localhost:5000/img/smiles?smiles=C
```

### Nomenclature Examples
```
Acetone:
http://localhost:5000/img/nomenclature?nomenclature=acetone

Benzene:
http://localhost:5000/img/nomenclature?nomenclature=benzene

Ethanol:
http://localhost:5000/img/nomenclature?nomenclature=ethanol

Formaldehyde:
http://localhost:5000/img/nomenclature?nomenclature=formaldehyde
```

---

## 📊 How Extension Routes Requests

| Input | Detected As | Endpoint | URL |
|-------|------------|----------|-----|
| `chem:acetone` | Nomenclature | `/img/nomenclature` | `?nomenclature=acetone` |
| `chem:CCO` | SMILES | `/img/smiles` | `?smiles=CCO` |
| `chem:benzene` | Nomenclature | `/img/nomenclature` | `?nomenclature=benzene` |
| `chem:c1ccccc1` | SMILES | `/img/smiles` | `?smiles=c1ccccc1` |

---

## 🧪 Test Commands

```powershell
# Health check
curl http://localhost:5000/health

# Get cache info
curl http://localhost:5000/cache-info

# Clear cache
curl -X DELETE http://localhost:5000/clear-cache

# Test SMILES (Ethanol)
curl "http://localhost:5000/img/smiles?smiles=CCO" > ethanol.svg

# Test Nomenclature (Acetone)
curl "http://localhost:5000/img/nomenclature?nomenclature=acetone" > acetone.svg
```

---

## 📁 Important Directories

```
MoleculeViewer/
├── server.js              ← Main server code
├── svg-cache/             ← Generated images stored here
├── package.json           ← Node dependencies
└── node_modules/          ← Installed packages
```

---

## 🔧 Configuration

### Change Port (if 5000 is taken):
Edit `server.js` line 8:
```javascript
const PORT = 5000;  // Change to 3000, 8000, etc
```

### Change Cache Size/Cleanup:
Edit `server.js` to modify cache management

### Change Image Size Defaults:
Query parameters control this:
```
&width=300&height=200
```

---

## 🐛 Frequent Issues & Fixes

| Issue | Fix |
|-------|-----|
| `Cannot find module 'express'` | `npm install` |
| `RDKit not installed` | `pip install rdkit` |
| `requests not installed` | `pip install requests` |
| Port 5000 in use | Change PORT in server.js |
| Images not loading | Restart server, reload extension |
| Nomenclature fails | Check spelling, wait for PubChem |

---

## 📋 Checklist

- [ ] `npm install` completed
- [ ] `pip install rdkit requests` completed
- [ ] `npm start` running (terminal shows green checkmark)
- [ ] Can open `http://localhost:5000/health` in browser
- [ ] Extension reloaded in chrome://extensions
- [ ] Tested `chem:acetone` in ChatGPT
- [ ] Images display inline ✅

---

## 🎯 Architecture Diagram

```
ChatGPT
   ↓
Extension detects: chem:acetone
   ↓
Regex matches: plain text = nomenclature
   ↓
Creates URL: http://localhost:5000/img/nomenclature?nomenclature=acetone
   ↓
Sets: <img src="...">
   ↓
Browser requests: GET /img/nomenclature?nomenclature=acetone
   ↓
Node.js Server
   ├→ Check cache: ❌ Not found
   ├→ Call Python: nomenclature → SMILES
   ├→ Call Python: SMILES → SVG
   ├→ Save to cache
   └→ Return SVG (Content-Type: image/svg+xml)
   ↓
Browser receives: SVG image
   ↓
Browser renders: Inline molecule structure
   ↓
User sees: 🧪 Molecule in ChatGPT! ✅
```

---

## 💡 Pro Tips

1. **Faster second load:** Images are cached, so repeating `chem:acetone` is instant
2. **Shareable URLs:** You can send URLs like `http://localhost:5000/img/smiles?smiles=CCO` to others (same network)
3. **Clear cache:** If server acts weird, `curl -X DELETE http://localhost:5000/clear-cache`
4. **Custom sizes:** Add `&width=400&height=300` to change image size
5. **Check logs:** Watch terminal while testing to see requests coming through

---

## 🔗 All Endpoints

```
GET  /                          → Info page
GET  /health                    → Server status
GET  /img/smiles                → Render SMILES
GET  /img/nomenclature          → Render nomenclature name
GET  /cache-info                → Cache statistics
DELETE /clear-cache             → Clear all cache
```

---

## 📞 Support

**Server won't start?**
- Check Node.js is installed: `node -v`
- Check npm is installed: `npm -v`
- Delete `node_modules/` and run `npm install` again

**Python errors?**
- Check Python is installed: `python --version`
- Check RDKit: `python -c "import rdkit; print(rdkit.__version__)"`
- Check requests: `python -c "import requests; print(requests.__version__)"`

**Still stuck?**
- Check terminal output for error messages
- Look in `/cache-info` endpoint for diagnostics
- Restart both server (`npm start`) and extension

---

**Server Status:** ✅ Running on http://localhost:5000

Ready to rock! 🎸
