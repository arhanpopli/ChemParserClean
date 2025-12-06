# Image Size Controls Feature

## Overview

The Chemistry Extension now includes interactive size controls for all rendered molecule images. Users can adjust the size of each molecule image using intuitive up/down arrow buttons that appear on hover.

## Features

### 1. Visual Size Controls
- **Up/Down Arrows**: Appear in the bottom-left corner of each molecule image when you hover over it
- **Smooth Transitions**: Controls fade in/out smoothly for a polished user experience
- **Responsive**: Each click increases or decreases the image size by 20 pixels
- **Size Constraints**: Images are constrained between 100px and 800px to maintain usability

### 2. Persistent Size Storage

Two developer options allow you to control how size preferences are saved:

#### Save Size Per Page
- When enabled, the extension remembers different sizes for the same molecule on different pages
- Perfect for customizing your view on specific pages
- Storage key includes both the page URL and the molecule identifier

#### Save Size By SMILES
- When enabled, all instances of the same molecule (same SMILES) use the same size across all pages
- Ideal for consistent molecule sizing across your entire browsing session
- Overrides the "Save Size Per Page" setting when both are enabled

### 3. Smart Storage System
- Uses `chrome.storage.local` for persistent storage
- Automatic size loading when molecules are rendered
- Instant saving when you adjust sizes
- Efficient storage using molecule identifiers (SMILES or nomenclature)

## Usage

### Basic Usage
1. **Enable Size Saving**: Open the extension popup and enable either:
   - "Save Size Per Page" (for page-specific sizes)
   - "Save Size By SMILES" (for global molecule sizes)

2. **Adjust Size**: Hover over any molecule image to reveal the size controls in the bottom-left corner

3. **Click Controls**:
   - Click ▲ (up arrow) to increase size
   - Click ▼ (down arrow) to decrease size

4. **Verify Persistence**: Reload the page to confirm your size preferences are saved

### Testing
Use the included test file to verify the feature:
```
file:///C:/Users/Kapil/Personal/PROJECTS/Chemparser/chem-extension/test-size-controls.html
```

The test file includes:
- Multiple instances of the same molecule to test SMILES-based sizing
- Different molecule types (SMILES, nomenclature)
- Complex structures to ensure broad compatibility
- Testing checklist for comprehensive validation

## Implementation Details

### Files Modified

#### 1. `content.js`
- **Added**: Complete size control system (lines 144-324)
  - `getImageKey()` - Generate unique keys for molecules
  - `getPageImageKey()` - Generate page-specific keys
  - `loadImageSize()` - Load saved size from storage
  - `saveImageSize()` - Save size to storage
  - `createSizeControls()` - Create UI buttons
  - `adjustImageSize()` - Handle size adjustments
  - `wrapImageWithSizeControls()` - Wrap images with controls

- **Modified**: All image replacement calls to use `wrapImageWithSizeControls()`
  - Line 1090-1095: MoleculeViewer images
  - Line 1304-1309: mol2chemfig chemfig fallback
  - Line 1406-1411: mol2chemfig data URI images
  - Line 1458-1463: mol2chemfig raw SVG images
  - Line 1533-1560: mol2chemfig Blob URL handling

#### 2. `popup.html`
- **Added**: New "Image Size Controls" section (lines 501-526)
  - "Save Size Per Page" toggle
  - "Save Size By SMILES" toggle
  - Explanatory info box

#### 3. `popup.js`
- **Added**: DOM element references (lines 45-46)
- **Added**: Default settings (lines 82-83)
- **Added**: Settings loading (lines 122-123)
- **Added**: Event handlers (lines 178-193)

### Storage Schema

#### Per-Page Storage
```javascript
{
  "page:https://example.com:smiles:CCO": {
    "width": 400,
    "height": 267
  }
}
```

#### SMILES-Based Storage
```javascript
{
  "smiles:CCO": {
    "width": 400,
    "height": 267
  }
}
```

### Configuration Constants

Defined in `content.js`:
```javascript
const SIZE_STEP = 20;        // Pixels per click
const MIN_SIZE = 100;        // Minimum image size
const MAX_SIZE = 800;        // Maximum image size
const DEFAULT_WIDTH = 300;   // Default width
const DEFAULT_HEIGHT = 200;  // Default height
```

## User Experience

### Visual Feedback
- Controls are hidden by default to keep the UI clean
- Smooth opacity transitions when hovering
- Dark semi-transparent buttons that work on light and dark backgrounds
- Clear up/down arrow indicators (▲ ▼)

### Size Adjustment Behavior
- Aspect ratio is maintained when resizing
- Size changes are applied immediately
- Saves happen automatically in the background
- No page reload required for size changes (but reload shows persisted sizes)

## Browser Compatibility

- Chrome (Manifest V3)
- Edge (Manifest V3)
- Other Chromium-based browsers supporting chrome.storage API

## Storage Limits

Chrome's storage API has the following limits:
- **sync.storage**: 100KB (used for settings)
- **local.storage**: 5MB (used for size data)

With our implementation:
- Each size entry is ~40-80 bytes
- Can store approximately 60,000+ molecule sizes before hitting limits
- Storage is efficient and automatically managed

## Future Enhancements

Potential improvements:
1. Add a "Reset to Default" button for each image
2. Add a "Clear All Sizes" option in the popup
3. Add keyboard shortcuts for size adjustment
4. Add a size slider in addition to buttons
5. Add size presets (small, medium, large)
6. Export/import size configurations

## Troubleshooting

### Controls Not Appearing
- Ensure the extension is loaded and enabled
- Check that at least one save option is enabled in the popup
- Verify that images are being rendered by the extension

### Sizes Not Persisting
- Make sure either "Save Size Per Page" or "Save Size By SMILES" is enabled
- Check browser console for any storage errors
- Verify chrome.storage permissions in manifest.json

### Images Too Large/Small
- Use the size controls to adjust to your preference
- Default size is 300x200 pixels
- Minimum size is 100x100 pixels
- Maximum size is 800x800 pixels

## Developer Notes

### Adding to Existing Projects
To add this feature to other Chrome extensions:

1. Copy the size control functions from `content.js` (lines 144-324)
2. Add the UI elements to your popup HTML
3. Add the event handlers to your popup.js
4. Modify your image replacement code to use `wrapImageWithSizeControls()`
5. Ensure `chrome.storage` permission is in manifest.json

### Customization
Easily customize the feature by modifying the constants:
- Change `SIZE_STEP` for larger/smaller increments
- Adjust `MIN_SIZE` and `MAX_SIZE` for different size ranges
- Modify button styles in `createSizeControls()` function
- Change storage behavior in `loadImageSize()` and `saveImageSize()`

## Credits

Implemented for the Chemparser project as specified in `MoleculeViewer/docs/Todolist.md` lines 4-6.

## License

Same license as the parent Chemparser project.
