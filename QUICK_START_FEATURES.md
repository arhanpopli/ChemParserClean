# Quick Start: New Features

## Three Features Fixed/Added

### 1. Size Control Arrows ✅ FIXED

**What:** Resize molecules with up/down arrow buttons

**How to use:**
1. Hover over any molecule image
2. Arrows appear in bottom-left corner
3. Click ▲ to make bigger (20px step)
4. Click ▼ to make smaller (20px step)
5. Size persists on page reload

**Default:** ENABLED (saves globally by SMILES)

**Settings:**
- `Save Size by SMILES` - Same size across all pages (enabled by default)
- `Save Size per Image` - Different size on each page

---

### 2. OPSIN 3D Stereochemistry ✅ WORKING

**What:** Get 3D SMILES with stereochemistry for chiral molecules

**How to use:**
1. Open extension popup
2. Enable "3D Stereochemistry (OPSIN)"
3. Reload page
4. Open console (F12) to see stereochemistry markers

**Default:** DISABLED (enable in popup)

**Example:**
```
Without OPSIN: chem:glucose: → C(C1C(C(C(C(O1)O)O)O)O)O
With OPSIN:    chem:glucose: → O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO
                                    ↑     ↑      ↑     ↑
                                    3D stereochemistry markers
```

**Console logs to look for:**
- 🔮 Purple: "Trying OPSIN 3D for name→3D SMILES conversion..."
- ✅ Green: "OPSIN 3D conversion SUCCESS: [SMILES]"

---

### 3. MolView.org 3D Viewer ✅ IMPLEMENTED

**What:** Interactive 3D molecular viewer using MolView.org

**How to use:**
1. Open extension popup
2. Enable "Enable 3D Viewer"
3. Reload page
4. Click 3D button on molecules
5. MolView iframe appears

**Default:** DISABLED (enable in popup)

**Features:**
- Interactive 3D rotation
- Multiple display modes (balls, sticks, vdw, wireframe)
- Toggle between 2D image and 3D model
- Uses SMILES when available (more accurate)

**URL format:**
```
https://embed.molview.org/v1/?mode=balls&smiles=CCO
```

---

## Testing All Features (2 Minutes)

### Quick Test

1. **Load test page:**
   ```
   file:///C:/Users/Kapil/Personal/PROJECTS/Chemparser/test_features_fixes.html
   ```

2. **Open extension popup and enable:**
   - ✅ "Save Size by SMILES" (already enabled)
   - ✅ "3D Stereochemistry (OPSIN)"
   - ✅ "Enable 3D Viewer"

3. **Reload page** (F5)

4. **Test size controls:**
   - Hover over ethanol → arrows appear
   - Click up arrow 3x → molecule grows
   - Reload page → size persists

5. **Test OPSIN 3D:**
   - Open console (F12)
   - Look for purple OPSIN logs
   - Check glucose SMILES has @ symbols

6. **Test MolView:**
   - Click 3D button on caffeine
   - Verify MolView iframe appears
   - Check console: URL should contain "embed.molview.org"

---

## Troubleshooting

### Size Controls Not Appearing

**Problem:** Arrows don't show on hover

**Solutions:**
1. Check popup settings: "Save Size by SMILES" should be ON
2. Hard reload page (Ctrl+Shift+R)
3. Check console for errors
4. Verify extension is enabled

### OPSIN 3D Not Working

**Problem:** Console shows errors or no OPSIN logs

**Solutions:**
1. Check backend is running: `curl http://localhost:8000/m2cf/opsin-3d?name=glucose`
2. Verify toggle is enabled in popup
3. Check console for red error messages
4. Try simpler molecule like "glucose"

**Backend test:**
```bash
curl "http://localhost:8000/m2cf/opsin-3d?name=glucose"
# Should return: {"smiles": "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", ...}
```

### MolView Not Loading

**Problem:** 3D button doesn't work or shows blank

**Solutions:**
1. Check internet connection (MolView is external service)
2. Open console and look for "📍 MolView URL" log
3. Verify URL contains "embed.molview.org"
4. Try opening MolView URL directly in browser
5. Check for CORS errors in console

---

## Files Modified

- `chem-extension/popup.js` - Changed `saveSizeBySMILES` default to `true`
- `chem-extension/content.js` - Updated 10 lines for size controls and MolView

---

## Console Logs Reference

### Size Controls
```
✅ Saved size for smiles:CCO: {width: 340, height: 220}
✅ Adjusted size: 340x220
```

### OPSIN 3D
```
🔮 Trying OPSIN 3D for name→3D SMILES conversion...
✅ OPSIN 3D conversion SUCCESS: O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO [3D]
```

### MolView
```
🔮 SHOWING MOLVIEW 3D VIEWER INLINE
📍 MolView URL: https://embed.molview.org/v1/?mode=balls&smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
✅ 3D viewer embedded inline
```

---

## Performance Notes

- **Size controls:** ~2KB added to content.js, minimal performance impact
- **OPSIN 3D:** Only active when enabled, adds ~100ms for name lookup
- **MolView:** External service, requires internet, ~1-2s load time

---

## Next Steps

1. Test in browser with real websites
2. Verify no regressions in existing features
3. Consider adding:
   - Keyboard shortcuts for size controls
   - Visual indicator when OPSIN 3D is active
   - Offline fallback for MolView
   - Size reset button

---

**For detailed technical documentation, see `FEATURES_FIXES_SUMMARY.md`**
