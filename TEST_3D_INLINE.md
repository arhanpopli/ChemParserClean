# Testing 3D Inline Viewer

## Steps to Test:

1. **Reload Extension**:
   - Go to `chrome://extensions/`
   - Click the reload button on your Chemistry extension

2. **Enable 3D Viewer**:
   - Click the extension icon
   - Go to "Developer Options" tab
   - Toggle ON "Enable 3D Viewer"

3. **Test on a Page**:
   - Create an HTML file or use any webpage
   - Type: `chem:histamine:`
   - The 3D viewer should appear **inline** (not popup)

4. **Test Toggle**:
   - Click the "📷 2D" button → switches to 2D image
   - Click the "🔮 3D" button → switches back to 3D viewer

## Expected Result:

```
Before: chem:histamine: (text)
After:  [3D VIEWER IN A BOX WITH PURPLE BORDER]
        ↑ Contains rotating 3D molecule
        ↑ Has toggle button top-right
        ↑ Has compound name bottom-left
```

## Features:

✅ **Inline embedding** - appears where you typed the formula
✅ **2D/3D toggle** - switch views without leaving the page
✅ **Styled container** - purple border, dark theme
✅ **Compound label** - shows molecule name
✅ **Size controls** - 600x400px default

## Troubleshooting:

If nothing appears:
1. Check browser console (F12) for errors
2. Verify server is running: http://localhost:5002/health
3. Try a different molecule: `chem:aspirin:` or `chem:caffeine:`
4. Check Developer Options has "Enable 3D Viewer" ON

## Test Molecules:

- `chem:histamine:` - small molecule
- `chem:aspirin:`
- `chem:caffeine:`
- `chem:dopamine:`
- `chem:glucose:`
