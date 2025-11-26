# PubChem Migration & Fixes Plan

## Current Status Analysis

### ✅ Already Using PubChem:
1. **Extension (content.js)** - Uses `smilesBridge()` function which fetches from PubChem API
   - Priority 1: Direct PubChem API
   - Priority 2: PubChem Autocomplete (for class names like "sphingomyelin")
   
2. **MoleculeViewer** - Uses PubChem via `nomenclature_to_smiles.py`
   - Endpoint: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON`

### ❌ Issues Found:

1. **Missing Functions in content.js**:
   - `getPubChemCID()` - Referenced but not defined (lines 2712, 2729, 1978)
   - `addHoverControls()` - Referenced but not defined (lines 733, 2785)

2. **3D Viewer**:
   - Currently tries to get CID using missing `getPubChemCID()` function
   - Should use PubChem to get SMILES first, then get CID from SMILES
   - For complex molecules like phosphatidylcholine, needs autocomplete fallback

3. **Size Controls**:
   - May not work for large molecules due to missing intrinsic width/height
   - Need to investigate `applyScaleToImage()` function

4. **mol2chemfig Backend**:
   - Currently only uses PubChem for numeric CID lookups
   - Should use PubChem for name→SMILES conversion like the extension does

## Implementation Plan

### Step 1: Add Missing Helper Functions

Add `getPubChemCID()` function that:
- Uses the existing `smilesBridge()` to get SMILES first
- Then queries PubChem to get CID from SMILES or name
- Handles autocomplete for complex molecules

Add `addHoverControls()` function that:
- Creates molecule name label
- Creates 3D viewer button
- Handles hover states

### Step 2: Update 3D Viewer Logic

Modify `show3DViewerInline()` to:
- Use `smilesBridge()` to get SMILES for complex molecules
- Fall back to direct name lookup for simple molecules
- Handle both IsomericSMILES (3D) and CanonicalSMILES (2D)

### Step 3: Fix Size Controls

Investigate and fix:
- Why `naturalWidth` might be undefined for large SVGs
- Add fallback dimensions for complex molecules
- Ensure scale calculations work with all molecule sizes

### Step 4: Update mol2chemfig Backend (Optional)

Add PubChem name→SMILES conversion similar to MoleculeViewer:
- Create Python helper function
- Use same PubChem API endpoints
- Add autocomplete fallback

## Code Changes Needed

### 1. Add `getPubChemCID()` function to content.js

```javascript
/**
 * Get PubChem CID for a molecule name or SMILES
 * Uses smilesBridge for complex molecules, direct lookup for simple ones
 */
async function getPubChemCID(nameOrSmiles) {
  console.log('%c🔍 Getting PubChem CID for:', 'color: #0088FF; font-weight: bold;', nameOrSmiles);
  
  try {
    // Try direct CID lookup first (works for names and SMILES)
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(nameOrSmiles)}/cids/JSON`;
    const data = await backgroundFetchJSON(url);
    
    if (data && data.IdentifierList && data.IdentifierList.CID && data.IdentifierList.CID.length > 0) {
      const cid = data.IdentifierList.CID[0];
      console.log('%c✅ Found CID:', 'color: #00FF00; font-weight: bold;', cid);
      return cid;
    }
  } catch (e) {
    console.warn('⚠️ Direct CID lookup failed:', e.message);
  }
  
  // Fallback: Use smilesBridge to get SMILES, then get CID from SMILES
  console.log('%c🌉 Using SMILES Bridge fallback...', 'color: #FF8800;');
  const bridgeResult = await smilesBridge(nameOrSmiles, { use3DSmiles: false });
  
  if (bridgeResult && bridgeResult.smiles) {
    console.log('%c✅ Got SMILES from bridge:', 'color: #00FF00;', bridgeResult.smiles);
    
    // Now get CID from SMILES
    try {
      const smilesUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(bridgeResult.smiles)}/cids/JSON`;
      const smilesData = await backgroundFetchJSON(smilesUrl);
      
      if (smilesData && smilesData.IdentifierList && smilesData.IdentifierList.CID && smilesData.IdentifierList.CID.length > 0) {
        const cid = smilesData.IdentifierList.CID[0];
        console.log('%c✅ Found CID from SMILES:', 'color: #00FF00; font-weight: bold;', cid);
        return cid;
      }
    } catch (e) {
      console.warn('⚠️ CID lookup from SMILES failed:', e.message);
    }
  }
  
  console.error('%c❌ Could not find CID for:', 'color: #FF0000; font-weight: bold;', nameOrSmiles);
  return null;
}
```

### 2. Add `addHoverControls()` function to content.js

```javascript
/**
 * Add hover controls (name label + 3D button) to molecule image container
 */
function addHoverControls(container, moleculeName, moleculeData) {
  // Create hover controls container
  const hoverControls = document.createElement('div');
  hoverControls.className = 'chem-hover-controls';
  hoverControls.style.cssText = `
    position: absolute;
    top: 4px;
    right: 4px;
    display: flex;
    flex-direction: column;
    gap: 4px;
    opacity: 0;
    transition: opacity 0.2s;
    pointer-events: auto;
    z-index: 1000;
  `;
  
  // Create molecule name label
  const nameLabel = document.createElement('div');
  nameLabel.className = 'chem-name-label';
  nameLabel.textContent = moleculeName;
  nameLabel.style.cssText = `
    background: rgba(0, 0, 0, 0.8);
    color: white;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 11px;
    font-family: monospace;
    max-width: 200px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  `;
  
  // Create 3D viewer button
  const viewer3DBtn = document.createElement('button');
  viewer3DBtn.className = 'chem-3d-btn';
  viewer3DBtn.innerHTML = '🔮 3D';
  viewer3DBtn.title = 'View in 3D';
  viewer3DBtn.style.cssText = `
    background: rgba(0, 0, 0, 0.8);
    color: white;
    border: none;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 11px;
    cursor: pointer;
    transition: background 0.2s;
  `;
  
  viewer3DBtn.addEventListener('mouseenter', () => {
    viewer3DBtn.style.background = 'rgba(0, 0, 0, 0.95)';
  });
  
  viewer3DBtn.addEventListener('mouseleave', () => {
    viewer3DBtn.style.background = 'rgba(0, 0, 0, 0.8)';
  });
  
  viewer3DBtn.addEventListener('click', async (e) => {
    e.stopPropagation();
    console.log('%c🔮 3D Viewer button clicked', 'color: #764ba2; font-weight: bold;');
    await show3DViewerInline(moleculeData, container);
  });
  
  hoverControls.appendChild(nameLabel);
  hoverControls.appendChild(viewer3DBtn);
  container.appendChild(hoverControls);
  
  // Show/hide on hover
  container.addEventListener('mouseenter', () => {
    hoverControls.style.opacity = '1';
  });
  
  container.addEventListener('mouseleave', () => {
    hoverControls.style.opacity = '0';
  });
}
```

### 3. Fix Size Controls for Large Molecules

Update `applyScaleToImage()` to handle missing naturalWidth:

```javascript
function applyScaleToImage(svgImg, scale) {
  // Try to get intrinsic width, with multiple fallbacks
  let intrinsicWidth = svgImg.naturalWidth || svgImg.width;
  
  // If still no width (can happen with complex SVGs), parse from SVG content
  if (!intrinsicWidth || intrinsicWidth === 0) {
    try {
      // Try to extract width from SVG src if it's a data URL
      if (svgImg.src && svgImg.src.startsWith('data:image/svg+xml')) {
        const svgContent = atob(svgImg.src.split(',')[1]);
        const widthMatch = svgContent.match(/width="(\d+)"/);
        if (widthMatch) {
          intrinsicWidth = parseInt(widthMatch[1]);
          console.log('%c📏 Extracted width from SVG:', 'color: #9c88ff;', intrinsicWidth);
        }
      }
    } catch (e) {
      console.warn('Could not extract SVG width:', e);
    }
  }
  
  // Final fallback: use default width
  if (!intrinsicWidth || intrinsicWidth === 0) {
    intrinsicWidth = 300;
    console.warn('%c⚠️ Using default width:', 'color: orange;', intrinsicWidth);
  }
  
  const newWidth = Math.round(intrinsicWidth * scale);

  // Only set width, let height: auto maintain aspect ratio
  svgImg.style.width = `${newWidth}px`;
  svgImg.style.height = 'auto';
  svgImg.style.maxWidth = 'none';
  svgImg.style.maxHeight = 'none';

  console.log(`Applied scale ${scale}x: intrinsic width ${intrinsicWidth}px → ${newWidth}px (height: auto)`);
}
```

## Testing Plan

1. Test with simple molecules (e.g., "aspirin", "caffeine")
2. Test with complex molecules (e.g., "phosphatidylcholine", "sphingomyelin", "insulin")
3. Test 3D viewer with both simple and complex molecules
4. Test size controls with various molecule sizes
5. Test dark mode inversion
6. Test stereochemistry (3D SMILES) option

## Priority

1. **HIGH**: Add missing `getPubChemCID()` and `addHoverControls()` functions
2. **HIGH**: Fix size controls for large molecules
3. **MEDIUM**: Test 3D viewer with complex molecules
4. **LOW**: Update mol2chemfig backend (already works via server-side conversion)
