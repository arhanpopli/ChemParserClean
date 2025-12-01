# ChemParser Quick Start Guide

## 🚀 Start Everything

```bash
1-start-all.bat
```

This starts **5 servers** automatically:
1. **Port 8000** - MolView PHP (viewer embeds)
2. **Port 8001** - Search API (**REQUIRED** - autocorrect for ALL engines)
3. **Port 5000** - MoleculeViewer (SVG rendering)
4. **Port 5001** - Mol2ChemFig (LaTeX chemistry)
5. **Port 5002** - PubChem Server (PubChem images)

## ✅ Verify Everything Works

### Test Search API (Most Important!)
Open in browser:
```
http://localhost:8001/search?q=aspirin
```

Should return JSON with compound data.

### Test Each Server
```
http://localhost:8000/                    ← MolView main page
http://localhost:5000/unified-interface.html  ← MoleculeViewer interface
http://localhost:5001/health              ← Mol2ChemFig health check
http://localhost:5002/health              ← PubChem health check
```

## 🧪 Test the Extension

### 1. Load Extension in Chrome
1. Open Chrome
2. Go to `chrome://extensions`
3. Enable "Developer mode"
4. Click "Load unpacked"
5. Select: `C:\Users\Kapil\Personal\PROJECTS\Chemparser\chem-extension`

### 2. Test Autocorrect (The Key Feature!)

Open any webpage with a text input, type:

**Typo Test (Should Autocorrect):**
```
chem:booxite:
```
✅ Should show: "✓ Autocorrected: booxite → Brookite"
✅ Should display mineral structure

**Protein Test:**
```
chem:rhinovirus:
```
✅ Should display protein structure (no autocorrect notice, name was correct)

**Regular Compound:**
```
chem:aspirin:
```
✅ Should display aspirin structure

**SMILES Test:**
```
chem:CCO:
```
✅ Should display ethanol structure

## 🎨 Select Rendering Engine

1. Click extension icon (puzzle piece in Chrome)
2. Select one of:
   - **MoleculeViewer** (default, best for most uses)
   - **mol2chemfig** (LaTeX-style)
   - **PubChem** (direct images)
   - **Client-Side** (offline, no server needed)

**NOTE:** ALL engines now use the Search API (port 8001) automatically!
Autocorrect works with every engine.

## ⚙️ Settings

### Enable 3D Stereochemistry
1. Open extension popup
2. Expand "MoleculeViewer Options"
3. Enable "Use 3D SMILES for stereochemistry"
4. Now queries return isomeric SMILES with stereo info

### Example Settings
```
✅ Extension Enabled
✅ Render mhchem LaTeX
✅ Render chemfig LaTeX
✅ Performance Mode
⬜ Developer Mode

Rendering Engine: MoleculeViewer ← Select this for best experience
```

## 🔍 How It Works

### Every Query Goes Through Search API
```
User: chem:booxite:
  ↓
Search API (port 8001):
  - Autocorrect typos
  - Identify type (compound/protein/mineral)
  - Return SMILES + metadata
  ↓
Selected Engine (MoleculeViewer/PubChem/etc.):
  - Render using corrected SMILES
  ↓
Show result with autocorrect notice (if typo was fixed)
```

## 🐛 Troubleshooting

### Problem: "localhost refused to connect"
**Solution**: Search API not running
```bash
cd Molview\molview
node search-server.js
```

### Problem: Extension not working
**Check:**
1. All servers running? (check console windows)
2. Extension loaded in Chrome?
3. Try reloading the webpage

### Problem: Wrong compound displayed
**Reason**: Search API found a different match
**Solution**:
- Use exact SMILES: `chem:CCO:`
- Check search API response: `http://localhost:8001/search?q=yourquery`

## 🎯 Test Cases

### Minerals (COD Database)
```
chem:quartz:          ← Should work
chem:brookite:        ← Should work
chem:booxite:         ← Should autocorrect to Brookite
```

### Proteins (PDB Database)
```
chem:rhinovirus:      ← Should find protein
chem:hemoglobin:      ← Should find protein
chem:insulin:         ← Should find protein
```

### Compounds (PubChem)
```
chem:aspirin:         ← Should work
chem:ethanol:         ← Should work
chem:caffeine:        ← Should work
```

### SMILES (Direct)
```
chem:CCO:             ← Ethanol
chem:CC(=O)OC1=CC=CC=C1C(=O)O:  ← Aspirin
```

## 🎓 Advanced Usage

### Force Specific Renderer
Add class to the text:
```html
<span class="chemfig-pubchem">chem:aspirin:</span>  ← Force PubChem
<span class="chemfig-mol2chemfig">chem:benzene:</span>  ← Force mol2chemfig
```

### Disable Autocorrect (Not Recommended)
Currently no way to disable - it's a core feature!
If you need exact matching, use SMILES instead of names.

## 📊 Console Logs

Open browser console (F12) to see:
```
🔍 UNIVERSAL SEARCH API - Preprocessing query
🔎 Query: "booxite"
✅ Search API Result: { ... }
🎯 Autocorrect: "booxite" → "Brookite"
🧬 SMILES: ...
📦 Type: mineral
```

## ✨ Cool Features

### 1. Typo Tolerance
```
chem:aspriin:  → Works! (autocorrects to Aspirin)
chem:etanol:   → Works! (autocorrects to Ethanol)
chem:benzen:   → Works! (autocorrects to Benzene)
```

### 2. Multi-Database Search
- PubChem: 100M+ compounds
- PDB: 200K+ protein structures
- COD: 500K+ mineral structures

### 3. Visual Feedback
Purple banner shows when autocorrect happens:
```
┌────────────────────────────────────┐
│ ✓ Autocorrected: booxite → Brookite │
└────────────────────────────────────┘
[Structure displayed below]
```

## 🛑 Stopping Servers

### Stop All Servers
```bash
util-stop-all.bat
```

### Manual Stop
Close all the console windows that opened from `1-start-all.bat`

## 📚 Documentation

- `UNIFIED_SEARCH_IMPLEMENTATION.md` - Complete technical docs
- `UNIFIED_SEARCH_GUIDE.md` - User guide
- `MOLVIEW_API_GUIDE.md` - MolView API reference

---

**Need Help?** Check the console logs (F12) for detailed error messages.
**Found a Bug?** The autocorrect similarity threshold can be adjusted in `search-server.js`
