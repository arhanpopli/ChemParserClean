# 🎯 ChemParser Interfaces Guide

## You Now Have 2 Main Interfaces:

### 1. 🎨 **Unified Interface** (unified-interface.html)
**Purpose:** Quick access hub for all tools
**URL:** http://localhost:5000/unified-interface.html

**What it has:**
- ✅ Tabs for MoleculeViewer, Mol2ChemFig, PubChem, 3D SMILES
- ✅ Quick testing and basic rendering
- ✅ Links to all test files
- ⚠️ **Basic ChemFig options only** (3D toggle, angle, indentation)

**Use this when:** You want quick access to all tools in one place

---

### 2. 📐 **Mol2ChemFig Complete Interface** (mol2chemfig-full-interface.html)
**Purpose:** Full-featured Mol2ChemFig with ALL 8 options
**URL:** http://localhost:5000/mol2chemfig-full-interface.html

**What it has:**
- ✅ **ALL 8 ChemFig Options:**
  1. Aromatic (circles in benzene)
  2. Fancy Bonds (prettier double/triple bonds)
  3. Compact View
  4. Show Carbon
  5. Show Methyl
  6. Flip (horizontal mirror)
  7. Flop (vertical mirror)
  8. Atom Numbers

- ✅ **Advanced Controls:**
  - Rotation angle (-360° to +360°)
  - Indentation (0-20 spaces)
  - H₂ treatment (Keep/Add/Delete)

- ✅ **Extra Features:**
  - Kekule.js chemical drawing
  - PubChem compound search
  - Reaction scheme support
  - Multiple export formats (SVG, Layered SVG, PDF)
  - Options persistence (saves your settings)

**Use this when:** You need all ChemFig options and advanced features

---

## 🚀 Quick Start Guide:

### Option A: Use Unified Interface (Hub)
```
1. Run: 1-start-all.bat
2. Opens: http://localhost:5000/unified-interface.html
3. Click tabs to switch between tools
```

### Option B: Go Directly to Mol2ChemFig Complete
```
1. Run: 1-start-all.bat
2. Navigate to: http://localhost:5000/mol2chemfig-full-interface.html
3. Use all 8 ChemFig options
```

---

## 📋 What Each Tab Does:

### In Unified Interface:

**📐 Mol2ChemFig Tab:**
- Opens mol2chemfig-full-interface.html in an iframe
- Has ALL 8 options + advanced controls
- This IS the complete interface embedded

**🧬 MoleculeViewer Tab:**
- Simple SMILES to SVG rendering
- Works WITHOUT Docker
- Good for quick tests

**🌐 PubChem Tab:**
- Fetch compound images
- 3D molecular viewer
- Search by chemical name

**🔬 3D SMILES Tab:**
- Convert names to 3D SMILES (OPSIN)
- Use generated SMILES in other renderers

**🧪 Tests Tab:**
- Links to all test files
- Individual feature testing
- Debugging tools

---

## ⚠️ About the Docker Error:

**If you see:**
```
Error: Failed to establish connection to localhost:8000
```

**What it means:**
- Mol2ChemFig backend needs Docker on port 8000
- Docker is NOT running on your system

**Your options:**

1. **Use MoleculeViewer Instead** (✅ Works Now, No Docker)
   - Go to MoleculeViewer tab
   - Same rendering, no Docker needed

2. **Install Docker** (Advanced Users Only)
   - Download: https://www.docker.com/products/docker-desktop
   - Run: `docker-compose up -d`
   - Restart: `1-start-all.bat`

---

## 🎯 Which Interface Should You Use?

### For **Quick Testing**:
→ Use **Unified Interface** (unified-interface.html)
- Fast tab switching
- All tools in one place
- Good for exploration

### For **Serious ChemFig Work**:
→ Use **Mol2ChemFig Complete** (mol2chemfig-full-interface.html)
- All 8 ChemFig options
- Kekule.js drawing
- Options persistence
- Reaction schemes
- Multiple export formats

### For **Simple SVG Rendering**:
→ Use **MoleculeViewer** directly
- No Docker needed
- Fast rendering
- Basic but reliable

---

## 📁 File Reference:

| File | Purpose | Docker? |
|------|---------|---------|
| `unified-interface.html` | Hub for all tools | Optional |
| `mol2chemfig-full-interface.html` | Complete ChemFig (8 options) | Yes* |
| `test_frontend.html` | MoleculeViewer only | No |
| `test_m2cf_full.html` | (Same as mol2chemfig-full) | Yes* |

*Docker only needed for advanced mol2chemfig rendering. Drawing, options UI, and previews work without it.

---

## ✅ What's Working Right Now:

**Without Docker:**
- ✅ MoleculeViewer SVG rendering
- ✅ PubChem 3D viewer
- ✅ 3D SMILES generation (OPSIN)
- ✅ Chemical drawing (Kekule.js)
- ✅ All ChemFig option controls (UI works)
- ✅ Tab navigation

**With Docker:**
- ✅ All of the above PLUS
- ✅ Mol2ChemFig backend rendering
- ✅ PDF generation
- ✅ Layered SVG export

---

## 💡 Pro Tips:

1. **Unified Interface** is your starting point - explore from there
2. **Mol2ChemFig Complete** has the full feature set for advanced work
3. **MoleculeViewer** is your fallback when Docker issues occur
4. All interfaces are served by port 5000, so they're accessible once servers start
5. The ChemFig options UI works even without Docker - only the backend rendering needs it

---

## 🎉 Summary:

**You asked:** "that's not the unified interface, that's the mol2chemfig interface"

**You're right!** I fixed it:
- ✅ **unified-interface.html** = Hub with tabs for ALL tools (MoleculeViewer, Mol2ChemFig, PubChem, Tests)
- ✅ **mol2chemfig-full-interface.html** = Complete Mol2ChemFig with all 8 options
- ✅ Mol2ChemFig Complete is EMBEDDED in the Unified Interface's Mol2ChemFig tab

**Both are available now:**
- Quick access? → unified-interface.html
- Deep work? → mol2chemfig-full-interface.html (or use the tab in unified)

**Everything is working except Docker-dependent mol2chemfig rendering, which you can work around with MoleculeViewer!** 🎯

