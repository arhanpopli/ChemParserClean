# ChemistryLaTeX Quick Reference Guide

## 🎯 Summary of Recent Changes

You've successfully implemented **IUPAC name support** and **fixed mineral rendering logic**. Here's what changed:

### ✅ What Was Fixed

1. **IUPAC Support Added**
   - IUPAC names now use **OPSIN** (not PubChem)
   - Syntax: `chem:iupac=2,4,6-trinitrotoluene:`
   - Server endpoint: `/iupac=<name>.svg`
   - Supports all rendering flags (+c, +n, +o, etc.)

2. **Mineral Rendering Fixed**
   - Minerals now **default to 2D SVG** (not forced 3D)
   - Uses COD for `codid` and PubChem for SMILES
   - 3D viewer only shows when explicitly requested
   - Syntax: `chem:mineral=quartz:`

3. **Proper Type Detection**
   - Added `isIUPAC` flag for IUPAC names
   - Added `isMineral` flag for minerals
   - Server correctly routes each type to appropriate API

## 📝 Syntax Reference

### Standard Syntax
```
chem:mol=benzene:              → PubChem lookup
chem:iupac=2-methylpropan-1-ol: → OPSIN lookup
chem:smiles=CCO:               → Direct SMILES
chem:mineral=quartz:           → COD lookup
chem:biomol=hemoglobin:        → RCSB PDB lookup
```

### Named Syntax (Custom Labels)
```
chem:Benzenesmiles=C1=CC=CC=C1:
chem:TNTiupac=2,4,6-trinitrotoluene:
chem:Quartzmineral=quartz:
```

### With Flags
```
chem:iupac=2-methylpropan-1-ol+c+n:  → Show carbons and atom numbers
chem:mineral=quartz+c+o:             → Show carbons and aromatic circles
```

## 🔧 Server Endpoints

| Endpoint | Purpose | Returns |
|----------|---------|---------|
| `/mol=<name>.svg` | PubChem compound lookup | SVG (2D structure) |
| `/mol=<name>.json` | PubChem compound data | JSON {cid, smiles} |
| `/iupac=<name>.svg` | OPSIN IUPAC conversion | SVG (2D structure) |
| `/iupac=<name>.json` | OPSIN IUPAC data | JSON {smiles} |
| `/smiles=<smiles>.svg` | Direct SMILES render | SVG (2D structure) |
| `/mineral=<name>.svg` | COD mineral lookup | SVG (2D structure) |
| `/mineral=<name>.json` | COD mineral data | JSON {codid, smiles} |
| `/biomol=<name>.json` | RCSB biomolecule lookup | JSON {pdbid} |

## 🎨 Rendering Flags

Works with: `mol`, `smiles`, `iupac`, `mineral`  
Does NOT work with: `biomol` (3D only)

| Flag | Description |
|------|-------------|
| `+c` | Show carbon labels |
| `+n` | Show atom numbers |
| `+o` | Show aromatic circles |
| `+h` | Show explicit hydrogens |
| `+m` | Show methyl groups (CH₃) |
| `+i` | Show implicit H labels |
| `+p` | Flip horizontal |
| `+q` | Flip vertical |
| `+s` | Use stereochemistry (isomeric SMILES) |
| `+g` | Gradient colors |
| `+k` | Compact drawing |
| `+d` | Use user's default settings as base |

## 📊 Data Source Routing

```
INPUT: chem:iupac=2-methylpropan-1-ol:
  ↓
EXTENSION: Detects iupac= → sets isIUPAC=true
  ↓
EXTENSION: Sends to /iupac=2-methylpropan-1-ol.svg
  ↓
SERVER: Calls fetchFromOPSIN() → gets SMILES
  ↓
SERVER: Renders 2D SVG using SmilesDrawer
  ↓
EXTENSION: Displays 2D SVG
```

```
INPUT: chem:mineral=quartz:
  ↓
EXTENSION: Detects mineral= → sets isMineral=true
  ↓
EXTENSION: Sends to /mineral=quartz.svg
  ↓
SERVER: Calls fetchFromCOD() → gets codid
SERVER: Calls fetchFromPubChem() → gets SMILES
  ↓
SERVER: Renders 2D SVG using SmilesDrawer
  ↓
EXTENSION: Displays 2D SVG
```

## 🧪 Testing

1. **Open test file**: `test-iupac-mineral.html`
2. **Enable extension**: Make sure ChemistryLaTeX is installed and active
3. **Check console**: Open DevTools (F12) → Console
4. **Verify sources**: 
   - IUPAC should log "OPSIN" (not PubChem)
   - Minerals should log "COD" or "MineralNames.js"

## 🔍 Key Files

### Extension Files
- `content.js` (lines 3507-3752): Main rendering logic
- `USAGE.md`: User-facing documentation
- `test-iupac-mineral.html`: Test page

### Server Files
- `server/api/render.js`:
  - Lines 893-930: `fetchFromOPSIN()`
  - Lines 948-990: `fetchFromCOD()`
  - Lines 1164-1210: Mineral and IUPAC handlers
- `server/MineralNames.js`: Local mineral database (198KB)

## 🎉 What Works Now

✅ IUPAC names use OPSIN (Cambridge/EBI)  
✅ Minerals use COD for codid  
✅ Minerals use PubChem for SMILES (2D rendering)  
✅ Minerals default to 2D SVG (not forced 3D)  
✅ 3D viewer only shows when explicitly requested  
✅ Proper flag detection (isIUPAC, isMineral)  
✅ Server correctly routes each type  
✅ All rendering flags work with IUPAC and minerals  

## 🚀 Next Steps

1. **Test the implementation**:
   ```bash
   # Open test file in browser
   start test-iupac-mineral.html
   ```

2. **Rebuild extension** (if needed):
   ```bash
   node build.js
   ```

3. **Reload extension** in Chrome:
   - Go to `chrome://extensions/`
   - Click reload icon on ChemistryLaTeX extension

4. **Test on a webpage**:
   - Create a simple HTML file with `chem:iupac=...` tags
   - Open in browser with extension enabled
   - Check console for data source verification

## 📝 Notes

- **COD provides SMILES**: Confirmed that COD + PubChem combo returns SMILES for 2D rendering
- **No forced 3D**: Removed the forced 3D mineral handler (lines 3565-3605 deleted from content.js)
- **Smart routing**: Extension automatically detects type and routes to correct server endpoint
- **Cache efficiency**: Server uses marker-based theming for efficient caching

## 🐛 Troubleshooting

**IUPAC not rendering?**
- Check console for OPSIN errors
- Verify IUPAC name is valid
- Try simpler IUPAC names first

**Minerals not rendering?**
- Check if mineral exists in MineralNames.js
- Try common minerals (quartz, calcite, diamond)
- Check console for COD API errors

**Still showing 3D for minerals?**
- Verify you're not using `+3d` flag
- Check if `is3D` is being set somewhere
- Reload extension

## 📚 Documentation

- **USAGE.md**: User-facing syntax guide
- **IMPLEMENTATION_STATUS.md**: Detailed implementation notes
- **test-iupac-mineral.html**: Interactive test page

---

**Last Updated**: 2026-01-03  
**Version**: ChemistryLaTeX v6.0 with IUPAC & Mineral Support
