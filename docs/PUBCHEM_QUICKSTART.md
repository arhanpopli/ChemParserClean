# PubChem Integration - Quick Start Guide

## 🚀 Get Started in 3 Steps

### Step 1: Start the PubChem Server

**Windows:**
```bash
cd C:\Users\Kapil\Personal\PROJECTS\Chemparser
start_pubchem.bat
```

**Linux/Mac:**
```bash
cd ~/Personal/PROJECTS/Chemparser
python pubchem_server.py
```

You should see:
```
======================================================================
✅ PubChem Server running on http://localhost:5002
======================================================================

📍 API Endpoints:
   Image (direct):  http://localhost:5002/pubchem/img/histamine
   Image (info):    http://localhost:5002/pubchem/image?name=histamine
   3D Model:        http://localhost:5002/pubchem/3d-model?name=histamine
   3D Viewer:       http://localhost:5002/pubchem/3d-viewer?name=histamine
   Compound Info:   http://localhost:5002/pubchem/info?name=histamine

💾 Cache Directory: C:\Users\Kapil\Personal\PROJECTS\Chemparser\pubchem-cache
======================================================================
```

### Step 2: Configure the Chrome Extension

1. Click the Chemistry Extension icon in Chrome
2. Under **"🧪 Rendering Engine"**, select **"🌐 PubChem"**
3. Configure PubChem options (optional):
   - Image Size: `large` (recommended)
   - Enable 3D Models: ✅ ON
   - Default View Type: `2D Structure`
4. Click away from the popup to save

### Step 3: Test It!

Open the test file:
```
C:\Users\Kapil\Personal\PROJECTS\Chemparser\test_pubchem.html
```

Or create a simple test page with:
```html
<p>Histamine: chem:histamine:</p>
<p>Caffeine: chem:caffeine:</p>
<p>Ethanol: chem:CCO:</p>
<p>Benzene: chem:c1ccccc1:</p>
```

Reload the page and you should see:
- Molecule images from PubChem
- **🔮 3D** buttons in the top-right corner
- Click 3D buttons to view interactive 3D models

## ✅ Verification Checklist

- [ ] PubChem server running on port 5002
- [ ] Extension shows "🌐 PubChem Server (localhost:5002)" in settings
- [ ] Test page shows molecule images
- [ ] 3D buttons appear on images
- [ ] Clicking 3D button opens viewer window
- [ ] Browser console shows no errors (F12)

## 🔍 Quick Tests

### Test 1: Direct API Access
Open in browser:
```
http://localhost:5002/pubchem/img/histamine
```
Should show a PNG image of histamine.

### Test 2: Compound Info
Open in browser:
```
http://localhost:5002/pubchem/info?name=caffeine
```
Should show JSON data with CID, URLs, and conformer info.

### Test 3: Extension Integration
1. Open any webpage
2. Add: `chem:aspirin:`
3. Reload page
4. Should see aspirin structure from PubChem

## 🐛 Troubleshooting

### Server won't start
```bash
# Check if dependencies are installed
pip install -r requirements_pubchem.txt

# Check if port 5002 is in use
netstat -ano | findstr :5002
```

### Images not showing
1. Check browser console (F12) for errors
2. Verify server is running: http://localhost:5002/health
3. Check extension is enabled
4. Verify "PubChem" mode is selected in extension settings

### 3D button not appearing
1. Open extension settings
2. Verify "Enable 3D Models" is ON
3. Reload the page

## 📊 Features Overview

| Feature | Status |
|---------|--------|
| 2D Structure Images | ✅ Working |
| 3D Model Viewer | ✅ Working |
| Chemical Nomenclature | ✅ Working |
| SMILES Notation | ✅ Working |
| Image Caching | ✅ Working |
| 3D Button | ✅ Working |
| Size Controls | ✅ Working |
| Dark Mode Support | ✅ Working |

## 🎯 Example Usage

### Chemical Names
```
chem:histamine:
chem:dopamine:
chem:serotonin:
chem:caffeine:
chem:aspirin:
```

### SMILES
```
chem:CCO:          (Ethanol)
chem:c1ccccc1:     (Benzene)
chem:CC(=O)O:      (Acetic Acid)
chem:Cc1ccccc1:    (Toluene)
```

### Complex Molecules
```
chem:adrenaline:
chem:penicillin:
chem:morphine:
chem:ascorbic acid:
```

## 🔧 Advanced Configuration

### Custom Image Sizes
In extension settings, select:
- `small` (100x100)
- `large` (300x300) - recommended
- `500x500` (extra large)

### 3D Model Options
- Enable/disable 3D buttons globally
- Choose between 2D or 3D default view
- Access different conformers via API

### Cache Management
View cache stats:
```
http://localhost:5002/cache-info
```

Clear cache:
```bash
curl -X DELETE http://localhost:5002/clear-cache
```

## 📚 More Information

For detailed documentation, see:
- `PUBCHEM_INTEGRATION.md` - Full documentation
- `pubchem_server.py` - Server source code
- `test_pubchem.html` - Comprehensive test suite

## 🎉 Success!

If you can see molecule images and 3D buttons, you're all set! Enjoy using PubChem integration for chemistry visualization.

For questions or issues, refer to the main Chemparser documentation.
