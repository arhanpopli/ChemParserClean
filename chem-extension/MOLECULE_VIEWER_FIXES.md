# MoleculeViewer Fixes - Summary

## Issues Fixed

### 1. Removed Autocorrection Notice ✓
**Problem**: The extension was displaying "✓ Autocorrected: alpha-glucose → Alpha-glucose" messages when molecules were autocorrected.

**Solution**: Commented out all `showAutocorrectNotice()` calls in `content.js`:
- Line 2724: Disabled in `loadMoleculeViewerImage()` for server-side rendering
- Line 2904: Disabled in `loadPubChemImage()` for PubChem images

The autocorrection still happens internally (for accuracy), but the user no longer sees the notification.

---

### 2. Non-Linear (Parabolic) Molecule Sizing ✓
**Problem**: Molecule size was increasing linearly with molecule complexity, causing very large molecules to become excessively large.

**Solution**: Implemented parabolic scaling in `size-controls.js`:

#### New Sizing Formula
```javascript
size = BASE_SIZE + sqrt(complexity) * SCALE_FACTOR * 50
```

#### How It Works
- **Small molecules**: Size increases quickly (steep part of parabola)
- **Large molecules**: Size increases slowly (flattening part of parabola)
- **Maximum cap**: 400px to prevent excessive sizes

#### Configuration
```javascript
BASE_SIZE = 150px          // Starting size for small molecules
SIZE_SCALE_FACTOR = 0.5    // Controls growth rate (lower = slower)
MAX_DEFAULT_SIZE = 400px   // Maximum default size
```

#### Complexity Estimation
- Uses SMILES string length as primary metric
- Falls back to nomenclature length (with 0.5x factor)
- Square root function creates the parabolic curve

#### Example Scaling
| Molecule Type | Complexity | Old Size | New Size |
|--------------|------------|----------|----------|
| Methane (C)  | ~10        | 300px    | ~165px   |
| Glucose      | ~30        | 300px    | ~215px   |
| Cholesterol  | ~100       | 300px    | ~300px   |
| Large Protein| ~500       | 300px    | ~400px (capped) |

---

## Files Modified

1. **`chem-extension/content.js`**
   - Commented out autocorrection notices (lines 2724, 2904)

2. **`chem-extension/size-controls.js`**
   - Added `calculateDefaultSize()` function with parabolic scaling
   - Updated `loadImageSize()` to use new calculation
   - Replaced fixed DEFAULT_WIDTH/HEIGHT with dynamic sizing

---

## Testing Recommendations

1. **Test small molecules** (e.g., "methane", "water") - should be reasonably sized
2. **Test medium molecules** (e.g., "glucose", "aspirin") - should scale moderately
3. **Test large molecules** (e.g., "cholesterol", "adrenaline") - should not be excessively large
4. **Test autocorrection** (e.g., "alpha-glucose", "rauxite") - should NOT show notice
5. **Verify manual sizing** - up/down arrows should still work correctly

---

## Technical Details

### Parabolic Curve Explanation
The square root function creates a curve where:
- **Early growth**: `sqrt(10) = 3.16` → rapid increase
- **Later growth**: `sqrt(100) = 10` → slower increase (only 3x despite 10x complexity)
- **Late growth**: `sqrt(500) = 22.4` → very slow increase (only 2.2x despite 5x complexity)

This matches the user's requirement: "increases then its rate of increase in size with respect to the actual size of the molecule decreases"

### Why Square Root?
- Creates the "first half of parabola" effect
- Natural diminishing returns
- Prevents extreme sizes while maintaining good scaling for small/medium molecules
- Mathematically simple and predictable
