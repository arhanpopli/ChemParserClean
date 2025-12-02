# Mineral 2D Rendering Implementation

## Overview
Updated the Chrome extension to render 2D images for minerals using SMILES when available from the search API, with formula-based fallback for minerals without SMILES.

## Implementation Details

### 1. Minerals with SMILES (e.g., Quartz - COD 5000035)
**Flow:**
1. Search API returns SMILES: `[Si]1([O])([O])O[Si]([O])([O])O[Si]([O])([O])O[Si]([O])([O])O[Si]([O])([O])O1`
2. Extension detects mineral has SMILES
3. Updates `moleculeData.smiles` with the mineral's SMILES
4. **Falls through to compound rendering logic** (STEP 3)
5. Renders using selected engine (MoleculeViewer, mol2chemfig, PubChem, or client-side)

**Result:** Minerals with SMILES get proper 2D chemical structure diagrams, just like regular compounds!

### 2. Minerals without SMILES (e.g., Niggliite - COD 9008913)
**Flow:**
1. Search API returns formula: `Pt Sn` (no SMILES available)
2. Extension detects mineral has formula but no SMILES
3. Creates a simple SVG displaying:
   - Mineral name (e.g., "Niggliite")
   - Chemical formula (e.g., "Pt Sn")
   - Instruction text: "Click 3D button for structure"
4. Adds a "3D" button overlay
5. Clicking 3D button loads the MolView 3D viewer

**Result:** Minerals without SMILES get a clean formula display with easy access to 3D structure.

## Code Changes

### File: `content.js`
**Location:** Lines 2847-3021 (mineral/biomolecule handling in `loadMoleculeViewerImage`)

**Key Logic:**
```javascript
// For minerals with SMILES, render them like compounds
if (searchData.compoundType === 'mineral' && searchData.smiles) {
  // Update moleculeData with SMILES
  moleculeData.smiles = searchData.smiles;
  moleculeData.nomenclature = searchData.correctedName;
  moleculeData.formula = searchData.searchResult.formula;
  
  // Fall through to compound rendering logic
  // (Renders using MoleculeViewer, mol2chemfig, etc.)
}
// For minerals without SMILES but with formula
else if (searchData.compoundType === 'mineral' && searchData.searchResult.formula) {
  // Create simple SVG with formula
  const svgContent = `
    <svg>
      <text>${mineralName}</text>
      <text>${formula}</text>
      <text>Click 3D button for structure</text>
    </svg>
  `;
  // Display with 3D button overlay
}
```

## Rendering Servers Support

All rendering servers now support minerals with SMILES:

### 1. **MoleculeViewer** (localhost:5000)
- Receives SMILES from extension
- Renders 2D structure using SmilesDrawer
- Works for minerals just like compounds

### 2. **mol2chemfig** (localhost:1000)
- Receives SMILES from extension
- Converts to chemfig/SVG
- Works for minerals just like compounds

### 3. **PubChem** (localhost:5002)
- For compounds only (minerals don't have PubChem CIDs)
- Not used for minerals

### 4. **Client-Side Renderer**
- Uses SmilesDrawer directly in browser
- Receives SMILES from extension
- Works for minerals just like compounds

## Example Outputs

### Quartz (with SMILES)
```
Input: "quartz"
Search API: Returns SMILES from COD
Extension: Renders 2D structure using selected engine
Display: Chemical structure diagram with bonds and atoms
```

### Niggliite (without SMILES)
```
Input: "niggliite"
Search API: Returns formula "Pt Sn" (no SMILES)
Extension: Creates formula SVG
Display: 
  ┌─────────────────┐
  │   Niggliite     │
  │     Pt Sn       │
  │ Click 3D for... │
  └─────────────────┘
  [3D Button]
```

## Testing

### Test Minerals with SMILES:
- Quartz: `chem:quartz:/d`
- Silicon dioxide: `chem:silicon dioxide:/d`

### Test Minerals without SMILES:
- Niggliite: `chem:niggliite:/d`
- Gold: `chem:gold:/d`

### Expected Behavior:
1. **With SMILES**: Should render a proper 2D chemical structure
2. **Without SMILES**: Should show formula display with 3D button
3. **All cases**: 3D button should work and load MolView 3D viewer

## Future Enhancements

1. **Better Formula Rendering**: Use subscripts for formulas (H₂O instead of H2O)
2. **Crystal Structure Hints**: Add small icons indicating crystal system
3. **Color Coding**: Different colors for different mineral classes
4. **SMILES Generation**: Attempt to generate SMILES from formulas for simple minerals
5. **Fallback to PubChem**: For minerals that might have PubChem entries
