# 3D Viewer Fallback Implementation Summary

## Problem
When trying to view large molecules like cardiolipin (CID 166177218) in 3D, PubChem returns a 404 error because no SDF 3D structure data exists for that compound.

**Error:**
```
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/166177218/sdf?record_type=3d 404 (Not Found)
Code: PUGREST.NotFound
Message: No records found for the given CID(s)
```

## Solution Implemented

### 1. Updated `pubchem-3d-wrapper.js`

Added automatic fallback logic that:
1. **Checks SDF availability** before loading PubChem viewer
2. **Falls back to MolView** with SMILES parameter when SDF data is unavailable
3. **Uses isomeric/stereochemical SMILES** for accurate 3D representation

**Key Changes:**
- Added `checkSDFAvailability(cid)` function to verify if 3D data exists
- Modified `loadPubChemViewer()` to accept `fallbackSmiles` parameter
- Automatically switches to MolView URL format: `https://embed.molview.org/v1/?smiles=<SMILES>`
- Updated loading message when using fallback

### 2. Next Step Required

**Update content.js** to pass SMILES parameter when creating wrapper URLs:

When calling `pubchem-3d-wrapper.html`, need to pass `smiles` parameter:
```javascript
const wrapperUrl = chrome.runtime.getURL('pubchem-3d-wrapper.html') + 
                   `?cid=${cid}&name=${encodeURIComponent(name)}&smiles=${encodeURIComponent(smiles)}`;
```

This ensures the wrapper has isomeric SMILES available for fallback.

## Fallback Flow

```
User clicks 3D button for "cardiolipin"
    ↓
pubchem-3d-wrapper loads with CID + SMILES
    ↓
checkSDFAvailability(166177218)
    ↓
HEAD request to PubChem SDF API → 404
    ↓
hasSDFData = false
    ↓
Fallback to MolView with SMILES:
https://embed.molview.org/v1/?smiles=CCCCC/C=C\C/C=C\CCCCCCCC(=O)OC[C@@H]...
    ↓
✅ 3D viewer loads successfully!
```

## Testing

1. **Test with compound that has 3D data** (e.g., aspirin CID 2244):
   - Should load PubChem 3D viewer normally
   
2. **Test with compound lacking 3D data** (e.g., cardiolipin CID 166177218):
   - Should automatically switch to MolView
   - Loading message should change to "Loading MolView 3D Viewer..."
   - MolView should render using SMILES

## Files Modified

- ✅ `pubchem-3d-wrapper.js` - Added fallback logic
- ⏳ `content.js` - Need to pass `smiles` parameter in wrapper URL

## Benefits

- **Seamless fallback**: Works transparently without user intervention
- **Better coverage**: Can view 3D for molecules without PubChem 3D data
- **Uses accurate SMILES**: Preserves stereochemistry with isomeric SMILES
- **Clear feedback**: Console logs show which viewer is being used
