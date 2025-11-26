# Client-Side Rendering Button Controls Fix

## Problem
When using **Client-Side** rendering mode, the molecule images didn't show the interactive control buttons (up/down arrows for size adjustment and 3D viewer button) that were visible when using MoleculeViewer or mol2chemfig renderers.

## Root Cause
The `wrapImageWithRotationContainer()` function had a conditional check:
```javascript
if (!rotation || rotation === 0 || rotation === 360) {
  addHoverControls(img, moleculeName, moleculeData);
  return img;  // ❌ Returns bare <img> element
}
```

When there was no rotation (the common case), it would:
1. Try to add hover controls directly to the `<img>` element
2. Return the bare `<img>` element

**The issue:** HTML `<img>` elements **cannot have child elements**. The `addHoverControls()` function appends buttons and overlays as children, which silently fails on `<img>` elements.

## Solution
Modified `wrapImageWithRotationContainer()` to **always create a wrapper div**, regardless of rotation:

```javascript
function wrapImageWithRotationContainer(img, rotation, moleculeName, moleculeData) {
  // Always create wrapper container (img elements cannot have children for hover controls)
  const wrapper = document.createElement('div');
  wrapper.className = 'chemfig-rotation-wrapper chemfig-diagram';
  wrapper.style.cssText = `
    display: inline-block !important;
    position: relative !important;
    vertical-align: middle !important;
    margin: 0 12px 0 0 !important;
    padding: 0 !important;
    visibility: visible !important;
    opacity: 1 !important;
    line-height: 0 !important;
    width: fit-content !important;
    height: fit-content !important;
  `;

  // Add the image to wrapper
  wrapper.appendChild(img);

  // Add hover controls to wrapper (not to img, since img can't have children)
  addHoverControls(wrapper, moleculeName, moleculeData);

  return wrapper;
}
```

## What This Fixes
Now when using **Client-Side** rendering, molecules will have the same interactive controls as other renderers:

### 🎯 Bottom-Left: Size Controls
- **▲ button** - Increase molecule size
- **▼ button** - Decrease molecule size
- Hidden by default, visible on hover

### 🎯 Top-Right: 3D Viewer Button
- **3D button** - Open interactive 3D molecular viewer
- Hidden by default, visible on hover
- Blue background with white text

### 🎯 Bottom-Right: Molecule Name
- Always visible label showing the molecule name
- Semi-transparent black background

## Testing
Use `test-client-side-buttons.html` to verify the fix:

1. Load the Chrome extension
2. Set rendering engine to "Client-Side"
3. Open `test-client-side-buttons.html`
4. Hover over molecules - all controls should appear
5. Test the size adjustment buttons (▲▼)
6. Test the 3D viewer button

## Files Changed
- **content.js** (lines 1561-1585) - Modified `wrapImageWithRotationContainer()` function

## Compatibility
This change maintains full compatibility with:
- MoleculeViewer renderer ✅
- mol2chemfig renderer ✅
- PubChem renderer ✅
- All molecular flags (+r, +fh, +fv, +s, etc.) ✅
- 3D viewer integration ✅

## Related Functions
- `wrapImageWithRotationContainer()` - Creates wrapper with controls
- `addHoverControls()` - Adds the actual button elements
- `renderClientSide()` - Main client-side rendering function that calls the wrapper
