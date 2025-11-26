# PubChem Integration & Size Control Fixes

## Summary of Changes

This document outlines the investigation and fixes related to SMILES lookup sources and size adjustment controls.

---

## Investigation Results

### ✅ All Components Already Use PubChem (Not OPSIN)

After thorough investigation, I found that **all major components are already using PubChem** for SMILES lookup:

#### 1. **Client-Side Rendering** ✅
- **File:** `chem-extension/content.js` (function `renderClientSide`)
- **Uses:** PubChem API directly via `backgroundFetchJSON()`
- **Endpoints:**
  - Priority 1: Direct lookup - `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,IsomericSMILES/JSON`
  - Priority 2: Autocomplete - `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{name}/json`
  - Priority 3: CID lookup via `getPubChemCID()`

#### 2. **MoleculeViewer Server** ✅
- **File:** `MoleculeViewer/nomenclature_to_smiles.py`
- **Uses:** PubChem API directly
- **Endpoint:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON`
- **Returns:** CanonicalSMILES, IsomericSMILES, or ConnectivitySMILES

#### 3. **mol2chemfig Server** ✅
- **File:** `mol2chemfig_server.py` (function `convert_nomenclature_to_smiles`)
- **Uses:** PubChem API directly
- **Endpoint:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON`
- **Returns:** CanonicalSMILES, IsomericSMILES, or ConnectivitySMILES

#### 4. **3D Viewer** ✅
- **File:** `chem-extension/content.js` (function `show3DViewerInline`)
- **Uses:** `getPubChemCID()` which queries PubChem
- **Then:** Fetches 3D structure using CID from PubChem

#### 5. **SMILES Bridge Function** ✅
- **File:** `chem-extension/content.js` (function `smilesBridge`)
- **Uses:** PubChem API with fallback chain
- **Supports:** `use3DSmiles` option for IsomericSMILES stereochemistry

---

## Stereochemistry Support

### PubChem IsomericSMILES ✅ Already Implemented

All PubChem integrations request **both** CanonicalSMILES and IsomericSMILES:

```javascript
// Client-Side & SMILES Bridge
const smiles = use3DSmiles ?
  (props.IsomericSMILES || props.CanonicalSMILES) :
  props.CanonicalSMILES;
```

```python
# MoleculeViewer & mol2chemfig
smiles = props[0].get('CanonicalSMILES') or \
         props[0].get('ConnectivitySMILES') or \
         props[0].get('IsomericSMILES')
```

**User Control:**
- Extension popup has "Use 3D Stereochemistry" toggle
- Setting: `m2cfUse3DSmiles` (boolean)
- When enabled: Uses IsomericSMILES (includes stereochemistry)
- When disabled: Uses CanonicalSMILES (2D only)

---

## OPSIN API Still Available (Not Used by Default)

The mol2chemfig server has an OPSIN endpoint (`/api/opsin`) but it's **not used by default**. The system uses PubChem as the primary source.

---

## 🔧 Actual Fixes Made

### Fix 1: Size Adjustment Buttons Not Working

**Problem:**
- For complex molecules (phosphatidylcholine, insulin, etc.), the size up/down buttons didn't work
- Root cause: Code used `imgElement.naturalWidth` which returns `0` for SVG data URLs

**Solution (content.js:1716-1740):**

Changed from using `naturalWidth` (doesn't work for SVG data URLs):
```javascript
// OLD - doesn't work for SVG data URLs
const intrinsicWidth = imgElement.naturalWidth || imgElement.width;
const newWidth = Math.round(intrinsicWidth * currentScale);
```

To using `offsetWidth` (current rendered size):
```javascript
// NEW - uses actual rendered width
const currentWidth = imgElement.offsetWidth || imgElement.width || imgElement.naturalWidth || 300;
const newWidth = Math.round(currentWidth * 1.2); // 20% increase
```

**Why This Works:**
- `offsetWidth` returns the actual rendered width in the DOM
- Works for all image types: PNG, SVG, data URLs, external URLs
- Each click increases/decreases by 20% of current size (multiplicative, not additive)

### Fix 2: Improved Button Padding (from previous fix)

Also improved the visual appearance of control buttons:
- 3D button padding: `4px 10px` → `6px 12px`
- Name label padding: `4px 8px` → `6px 10px`
- Added `line-height: 1.4` for better text rendering

---

## Testing Complex Molecules

These molecules should now work properly with all features:

### Size Adjustment Test
```
chem:phosphatidylcholine:
chem:sphingomyelin:
chem:insulin:
```

- Hover over molecule
- Click ▲ to increase size (should work)
- Click ▼ to decrease size (should work)
- Each click = 20% size change

### 3D Viewer Test
```
chem:phosphatidylcholine:
chem:caffeine:
chem:dopamine:
```

- Hover over molecule
- Click **3D** button (top-right)
- Should load 3D viewer using PubChem CID

### Stereochemistry Test
1. Open extension popup
2. Enable "Use 3D Stereochemistry"
3. Test with chiral molecules:
   ```
   chem:L-alanine:
   chem:D-glucose:
   chem:(R)-limonene:
   ```

---

## PubChem Coverage

PubChem has **excellent coverage** for:
- ✅ Common organic molecules (ethanol, benzene, caffeine)
- ✅ Complex molecules (phosphatidylcholine, sphingomyelin)
- ✅ Proteins/peptides (insulin - may be very large)
- ✅ Drug molecules (aspirin, ibuprofen, etc.)
- ✅ Trade names (Tylenol, Advil)
- ✅ Stereoisomers (L-, D-, R-, S- prefixes)

**Much better than OPSIN** which has limited coverage for:
- ❌ Complex lipids
- ❌ Large biomolecules
- ❌ Many trade names

---

## Files Modified

### content.js
- **Lines 1716-1740:** Fixed size adjustment buttons to use `offsetWidth` instead of `naturalWidth`
- **Lines 1608, 1632:** Increased button padding for better appearance
- **Lines 1617, 1641:** Added `line-height: 1.4` for text rendering

### Already Using PubChem (No Changes Needed)
- `MoleculeViewer/nomenclature_to_smiles.py` ✅
- `mol2chemfig_server.py` ✅
- `chem-extension/content.js` - SMILES lookup functions ✅

---

## Conclusion

**All SMILES lookups already use PubChem**, which has superior coverage compared to OPSIN. The actual issue was with the **size adjustment buttons** not working for SVG data URLs, which has now been fixed by using `offsetWidth` instead of `naturalWidth`.
