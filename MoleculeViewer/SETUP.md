# MoleculeViewer Node.js Server

## ✅ DONE - Server Running

Your Node.js server is **running on http://localhost:5000** and ready to serve molecule images.

---

## 🚀 Quick Start

### Already Installed:
- ✅ Node.js server (`server.js`)
- ✅ Python helpers for chemistry
- ✅ npm dependencies
- ✅ Python packages (rdkit, requests)
- ✅ Extension updated

### Next: Reload Extension

1. Go to **chrome://extensions**
2. Click **reload** icon on your extension
3. Open **ChatGPT** and type:
   ```
   chem:acetone
   chem:benzene
   chem:CCO
   ```

You should see **inline molecule structures** 🧪

---

## 🔗 API Endpoints

| Format | URL |
|--------|-----|
| **SMILES** | `http://localhost:5000/img/smiles?smiles=CCO` |
| **Name** | `http://localhost:5000/img/nomenclature?nomenclature=acetone` |
| **Health** | `http://localhost:5000/health` |

---

## 📊 How It Works

```
chem:acetone (typed in ChatGPT)
    ↓
Extension creates: http://localhost:5000/img/nomenclature?nomenclature=acetone
    ↓
Server checks cache → renders if needed
    ↓
Returns SVG image
    ↓
Browser displays inline 🧪
```

---

## ⚡ Files Created

| File | Purpose |
|------|---------|
| `server.js` | Main Node.js server |
| `generate_svg.py` | SMILES → SVG |
| `nomenclature_to_smiles.py` | Name → SMILES |
| `svg-cache/` | Image cache |

---

## 🎯 That's It!

Reload extension → Test in ChatGPT → See molecules! ✅
