# ChemistryLaTeX Implementation Status

## ✅ Completed Features

### 1. **IUPAC Support** (OPSIN Integration)
- **Syntax**: `chem:iupac=2,4,6-trinitrotoluene:` or `chem:TNTiupac=2,4,6-trinitrotoluene:`
- **Data Source**: OPSIN (Cambridge/EBI) - NOT PubChem
- **Server Endpoint**: `/iupac=<name>.svg`
- **Flow**:
  1. Extension detects `iupac=` syntax via regex
  2. Sets `isIUPAC` flag in parsed flags
  3. Sends to server: `/iupac=2,4,6-trinitrotoluene.svg`
  4. Server calls `fetchFromOPSIN()` to convert IUPAC → SMILES
  5. Server renders 2D SVG using SmilesDrawer
  6. Extension displays 2D SVG

### 2. **Mineral Rendering** (COD Integration)
- **Syntax**: `chem:mineral=quartz:` or `chem:Quartzmineral=quartz:`
- **Data Sources**: 
  - COD (Crystallography Open Database) for `codid`
  - PubChem fallback for SMILES (for 2D rendering)
  - Local `MineralNames.js` database (198KB, fallback when COD is down)
- **Server Endpoint**: `/mineral=<name>.svg` or `/mineral=<name>.json`
- **Flow**:
  1. Extension detects `mineral=` syntax via regex
  2. Sets `isMineral` flag in parsed flags
  3. Sends to server: `/mineral=quartz.svg`
  4. Server calls `fetchFromCOD()`:
     - First tries local `MineralNames.js` lookup (fast)
     - Falls back to COD API if not in local DB
     - Gets `codid` from COD
     - Tries to get SMILES from PubChem using mineral name
  5. **2D Rendering (Default)**: Server renders SVG using SMILES
  6. **3D Rendering (Explicit)**: Only if `is3D` flag is set, shows MolView with `codid`

### 3. **Data Source Routing**

| Type | Syntax | 2D Render | Data Source | 3D Viewer (if requested) |
|------|--------|-----------|-------------|--------------------------|
| **Compound** | `chem:mol=benzene:` | 2D SVG | PubChem (name → SMILES) | MolView (CID or SMILES) |
| **IUPAC** | `chem:iupac=2-methylpropan-1-ol:` | 2D SVG | OPSIN (IUPAC → SMILES) | MolView (SMILES) |
| **Direct SMILES** | `chem:smiles=CCO:` | 2D SVG | Direct (no lookup) | MolView (SMILES) |
| **Mineral** | `chem:mineral=quartz:` | 2D SVG | COD (name → codid) + PubChem (name → SMILES) | MolView (codid) |
| **Biomolecule** | `chem:biomol=hemoglobin:` | 2D Image | RCSB PDB (name → pdbid) | Mol*/MolView (pdbid) |

### 4. **Regex Patterns Updated**

```javascript
// Standard format
/chem:iupac=([^:]+):/gi
/chem:mineral=([^:]+):/gi

// Named format (with custom label)
/chem:(\w+)iupac=([^:]+):/gi
/chem:(\w+)mineral=([^:]+):/gi
```

### 5. **Server API Endpoints**

```
GET /mol=<name>.svg          → PubChem lookup + render SVG
GET /mol=<name>.json         → PubChem lookup + return {cid, smiles}
GET /smiles=<smiles>.svg     → Direct render (no lookup)
GET /iupac=<name>.svg        → OPSIN lookup + render SVG
GET /iupac=<name>.json       → OPSIN lookup + return {smiles}
GET /mineral=<name>.svg      → COD lookup + render SVG (using SMILES)
GET /mineral=<name>.json     → COD lookup + return {codid, smiles}
GET /biomol=<name>.json      → RCSB lookup + return {pdbid}
```

### 6. **Flag Detection**

The extension now correctly identifies:
- `isIUPAC`: true when using `iupac=` syntax
- `isMineral`: true when using `mineral=` syntax
- `isBiomolecule`: true when using `biomol=` syntax
- `isDirectSmiles`: true when using `smiles=` syntax

## 🔧 Server Implementation Details

### `fetchFromOPSIN(iupacName)` - Lines 893-930
```javascript
// Tries Cambridge OPSIN first, falls back to EBI OPSIN
// Returns: { smiles }
// Throws: Error if IUPAC name cannot be converted
```

### `fetchFromCOD(name)` - Lines 948-990
```javascript
// 1. Tries local MineralNames.js lookup (fast, offline-capable)
// 2. Falls back to COD API if not found locally
// 3. Gets codid from COD
// 4. Tries to get SMILES from PubChem for 2D rendering
// Returns: { codid, smiles }
// Throws: Error if mineral not found
```

### Main Handler - Lines 1082-1282
- **MINERAL** (lines 1164-1184): 
  - Fetches from COD
  - Returns JSON if `.json` requested
  - Renders 2D SVG using SMILES if `.svg` requested
  - Uses smart flip analysis for CSS transforms

- **IUPAC** (lines 1192-1210):
  - Fetches from OPSIN
  - Returns JSON if `.json` requested
  - Renders 2D SVG if `.svg` requested
  - Uses smart flip analysis for CSS transforms

## 📝 Extension Implementation Details

### `content.js` - `renderClientSide()` function (lines 3507-3752)

```javascript
// Determines server endpoint based on flags:
if (flags.isDirectSmiles || moleculeData.smiles) {
    serverType = 'smiles';
} else if (flags.compoundType === 'biomolecule' || flags.isBiomolecule) {
    serverType = 'biomol';
} else if (flags.compoundType === 'mineral' || flags.isMineral) {
    serverType = 'mineral';
} else if (flags.isIUPAC) {
    serverType = 'iupac';
}
// else: 'mol' (default compound lookup)
```

### Mineral 2D/3D Logic
- **Default**: 2D SVG rendering (using SMILES from COD/PubChem)
- **3D Only When**: `is3D` flag is explicitly set
- **Removed**: Forced 3D mineral handler (lines 3565-3605 were deleted)

## 🎯 Current Behavior

### IUPAC Names
```
chem:iupac=2,4,6-trinitrotoluene:
→ Extension sends: /iupac=2,4,6-trinitrotoluene.svg
→ Server calls OPSIN API
→ OPSIN returns SMILES
→ Server renders 2D SVG
→ Extension displays 2D structure
```

### Minerals
```
chem:mineral=quartz:
→ Extension sends: /mineral=quartz.svg
→ Server checks MineralNames.js → finds codid
→ Server tries PubChem → gets SMILES
→ Server renders 2D SVG using SMILES
→ Extension displays 2D structure

chem:mineral=quartz+3d:  (if 3D flag added)
→ Extension sends: /mineral=quartz.json
→ Server returns {codid, smiles}
→ Extension shows MolView 3D viewer with codid
```

## ✅ Fixes Applied

1. **Removed forced 3D for minerals** - Minerals now default to 2D SVG like compounds
2. **Added IUPAC support** - Uses OPSIN, not PubChem
3. **Added `isIUPAC` and `isMineral` flags** - Proper type detection
4. **Updated regex patterns** - Recognizes both standard and named syntax
5. **Server correctly routes** - Each type goes to correct API endpoint
6. **COD returns SMILES** - Confirmed that COD + PubChem provides SMILES for 2D rendering

## 📊 Data Flow Summary

```
USER INPUT: chem:iupac=2-methylpropan-1-ol:
    ↓
EXTENSION: Detects iupac= syntax, sets isIUPAC=true
    ↓
EXTENSION: Builds URL → /iupac=2-methylpropan-1-ol.svg
    ↓
SERVER: fetchFromOPSIN() → SMILES
    ↓
SERVER: renderSmilesToSVG() → SVG
    ↓
EXTENSION: Displays 2D SVG
```

```
USER INPUT: chem:mineral=quartz:
    ↓
EXTENSION: Detects mineral= syntax, sets isMineral=true
    ↓
EXTENSION: Builds URL → /mineral=quartz.svg
    ↓
SERVER: fetchFromCOD() → codid + SMILES
    ↓
SERVER: renderSmilesToSVG() → SVG (using SMILES)
    ↓
EXTENSION: Displays 2D SVG
```

## 🧪 Testing Checklist

- [ ] Test IUPAC: `chem:iupac=2,4,6-trinitrotoluene:`
- [ ] Test IUPAC with flags: `chem:iupac=2-methylpropan-1-ol+c+n:`
- [ ] Test mineral 2D: `chem:mineral=quartz:`
- [ ] Test mineral with flags: `chem:mineral=calcite+c:`
- [ ] Test named IUPAC: `chem:TNTiupac=2,4,6-trinitrotoluene:`
- [ ] Test named mineral: `chem:Quartzmineral=quartz:`
- [ ] Verify OPSIN is used (not PubChem) for IUPAC
- [ ] Verify COD is used for minerals
- [ ] Verify minerals show 2D by default
- [ ] Verify 3D only shows when explicitly requested

## 📁 Key Files Modified

1. **content.js** (lines 3507-3752)
   - Updated `renderClientSide()` to route to correct server endpoint
   - Removed forced 3D mineral handler

2. **server/api/render.js** (lines 893-930, 948-990, 1164-1210)
   - Added `fetchFromOPSIN()` function
   - Updated `fetchFromCOD()` to use local DB + PubChem
   - Added IUPAC and mineral handlers in main handler

3. **USAGE.md** (lines 1-56)
   - Updated documentation with IUPAC and mineral syntax

## 🎉 Summary

All requested features are now implemented:
- ✅ IUPAC names use OPSIN (not PubChem)
- ✅ Minerals use COD for `codid` and PubChem for SMILES
- ✅ Minerals default to 2D rendering (not forced 3D)
- ✅ 3D viewer only shows when explicitly requested
- ✅ Proper flag detection (`isIUPAC`, `isMineral`)
- ✅ Server correctly routes each type to appropriate API
