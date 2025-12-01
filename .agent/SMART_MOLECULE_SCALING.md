# Smart Molecule Scaling Implementation

## Problem
When switching to MoleculeViewer, large molecules (like insulin) were being scaled linearly, causing them to overflow the page boundaries. The previous implementation used a fixed 1.5x scale for all molecules regardless of their size.

## Solution
Implemented an intelligent, non-linear scaling system that automatically adjusts the display scale based on the intrinsic dimensions of each molecule, **with hard maximum size limits** to prevent page overflow.

### Key Changes

#### 1. Smart Scaling Algorithm (`wrapImageWithSizeControls`)
- **Location**: `content.js` lines ~1024-1098
- **Algorithm**: Uses a parabolic (square root) scaling curve
- **Formula**: `scale = baseScale × (referenceSize / moleculeSize)^0.5`
  - `baseScale`: 1.5 (reference scale for medium molecules)
  - `referenceSize`: 300px (diagonal reference size)
  - `moleculeSize`: Calculated as diagonal (√(width² + height²))
  - `exponent`: 0.5 (square root for parabolic curve)

#### 2. Hard Maximum Size Limits
**CRITICAL**: All molecules are constrained to maximum display dimensions:
- **Max Width**: 600px
- **Max Height**: 500px

These limits are enforced via CSS `max-width` and `max-height` properties and apply **regardless of calculated scale**. This ensures that even if the smart scaling algorithm calculates a large scale, the molecule will never overflow the page.

**Locations**:
- `applyScaleToImage()` - Lines 1182-1193 (applies to all scaled images)
- `loadMoleculeViewerImage()` - Lines 2926-2927 (MoleculeViewer SVG images)
- `renderClientSide()` - Lines 2612-2613 (SmilesDrawer rendering)

#### 3. Scale Ranges by Molecule Size
- **Very small molecules** (< 150px diagonal): 1.7x - 2.0x scale
  - Examples: Water, methane, benzene
  - **Capped at**: 600x500px max
- **Small molecules** (150-250px): 1.5x - 1.7x scale
  - Examples: Glucose, ethanol
  - **Capped at**: 600x500px max
- **Medium molecules** (250-500px): 1.0x - 1.5x scale
  - Examples: Caffeine, aspirin
  - **Capped at**: 600x500px max
- **Large molecules** (500-800px): 0.6x - 1.0x scale
  - Examples: Cholesterol, complex steroids
  - **Capped at**: 600x500px max
- **Very large molecules** (> 800px): 0.4x - 0.6x scale
  - Examples: Insulin, proteins, large polymers
  - **Capped at**: 600x500px max

#### 4. Dimension Extraction
The algorithm extracts molecule dimensions from:
1. `naturalWidth` and `naturalHeight` properties
2. SVG `width` and `height` attributes (from data URL)
3. SVG `viewBox` attribute (fallback)

#### 5. Updated Functions

**`loadImageSize()`** - Lines 812-846
- Changed to return empty object `{}` instead of `{ scale: 1.5 }`
- Allows smart scaling algorithm to calculate appropriate defaults
- Still respects saved user preferences

**`adjustImageSize()`** - Lines 939-941
- Updated fallback from 1.5x to 1.0x for more neutral default
- Should rarely be used since scale is set by `wrapImageWithSizeControls`

**`applyScaleToImage()`** - Lines 1180-1195
- Added hard maximum size limits (600x500px)
- Ensures molecules never overflow page

## Benefits

1. **No Overflow**: Hard limits guarantee molecules never exceed 600x500px
2. **Better Visibility**: Small molecules scale up for easier viewing
3. **Smooth Transition**: Parabolic curve provides gradual scaling changes
4. **User Control**: Users can still manually adjust size with +/- buttons (within limits)
5. **Persistent Preferences**: Manual adjustments are saved per molecule
6. **Consistent Sizing**: All rendering methods (MoleculeViewer, SmilesDrawer, PubChem) use same limits

## Testing Recommendations

Test with molecules of varying sizes:
- **Small**: `chem:water:`, `chem:methane:`, `chem:benzene:`
- **Medium**: `chem:glucose:`, `chem:caffeine:`, `chem:aspirin:`
- **Large**: `chem:cholesterol:`, `chem:testosterone:`
- **Very Large**: `chem:insulin:`, `chem:hemoglobin:` ← Should now fit within 600x500px!

## Technical Notes

- The scaling is calculated once when the molecule is first rendered
- Saved user preferences override the smart scaling (but still respect hard limits)
- The algorithm uses logarithmic dampening (square root) to prevent aggressive scaling
- Minimum scale: 0.4x (prevents molecules from becoming too small)
- Maximum scale: 2.0x (prevents excessive enlargement)
- **Hard maximum display size: 600px × 500px** (enforced via CSS, cannot be exceeded)
