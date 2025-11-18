# ✅ Final Checklist: Make Options Work

## The Issue
The extension finds SMILES but **doesn't tell the server what options are enabled**, so it generates images with default settings.

## The Fix Applied
✅ **Backend** (`MoleculeViewer/app/api.py`) - Now reads options from URL parameters  
✅ **Extension** (`chem-extension/content.js`) - Now sends options as URL parameters  
✅ **Logging** - Added detailed console logs to debug the flow

---

## Step-by-Step Testing (Do This Now!)

### 1️⃣ Reload Chrome Extension
```
chrome://extensions
→ Find "Chemistry Formula Renderer"
→ Click the circular RELOAD button
```
**WHY**: Extension needs to load the new `content.js` code

### 2️⃣ Verify Settings Are Saved
```
1. Click extension icon (popup opens)
2. Select "MoleculeViewer" renderer
3. Enable "Aromatic Circles" → Should see "Aromatic circles enabled. Reload page to apply."
4. Close popup
5. Open popup again → "Aromatic Circles" should still be checked
```
**IF NOT**: Settings aren't saving. Check popup.js event listeners.

### 3️⃣ Test on Debug Page
```
1. Open: file:///C:/Users/Kapil/Personal/PROJECTS/Mol2chemfig/extension_debug_test.html
2. Press F12 (DevTools)
3. Go to Console tab
4. Press F5 (reload page)
5. Look for colored console logs
```

### 4️⃣ What You Should See in Console

**Good Output (Working):**
```
🧪 Using MOLECULEVIEWER renderer engine
📤 Using SMILES endpoint
⚙️ Rendering Options: {aromaticCircles: true, showCarbons: false, ...}
API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true&...
🌐 Fetching from backend...
✅ Got response from backend
📍 Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_xyz.svg
✅ Image loaded successfully
```

**Bad Output (Not Working):**
```
🧪 Using MOLECULEVIEWER renderer engine
📤 Using SMILES endpoint
⚙️ Rendering Options: {aromaticCircles: false, ...}  ← All false!
API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=false&...
📍 Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_xyz.svg  ← No "aromatic_circles" in URL!
```

---

## Debugging Based on Console Output

### Problem: ⚙️ Rendering Options shows all `false`
**Cause**: Settings not loaded from chrome.storage  
**Fix**:
1. Open DevTools → Application tab → Storage → chrome.storage.sync
2. Check if `aromaticCircles` key exists and is `true`
3. If not, extension popup isn't saving correctly
4. Check popup.js has event listeners (lines ~160-220)

### Problem: API URL missing `aromatic_circles=true`
**Cause**: Extension using old code  
**Fix**:
1. Go to `chrome://extensions`
2. Click RELOAD on extension
3. Hard refresh test page: `Ctrl+Shift+R`

### Problem: Cache URL doesn't have `aromatic_circles` in filename
**Cause**: Backend using old code  
**Fix**:
```powershell
# Kill Flask
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force

# Clear cache
cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"
Remove-Item app/__pycache__ -Recurse -Force -ErrorAction SilentlyContinue

# Restart
python run_server.py
```

### Problem: No console logs at all
**Cause**: Extension not running  
**Fix**:
1. Check extension is **enabled** in `chrome://extensions`
2. Check content script is injected (look for "Extension Debug Test Page Loaded" log)
3. Try reloading extension and page

---

## Expected Visual Differences

### Benzene (`c1ccccc1`)
- **Without Aromatic Circles**: ![Benzene Kekule](Alternating single/double bonds)
- **With Aromatic Circles**: ![Benzene Circle](Circle inside hexagon)

### Propane (`CCC`)
- **Without Show Methyls**: ![Propane line](Just a line structure)
- **With Show Methyls**: ![Propane labeled](CH₃-CH₂-CH₃ labels visible)

---

## Files Modified Summary

```
✅ MoleculeViewer/app/api.py
   Lines 657-680: img_smiles() - parse options from query params
   Lines 722-745: img_nomenclature() - parse options from query params

✅ chem-extension/content.js
   Lines 806-820: Build URLSearchParams with options (nomenclature)
   Lines 827-842: Build URLSearchParams with options (SMILES)
   Lines 845-860: Log options being used (NEW)
   Lines 912: Log cache URL (already existed)

✅ chem-extension/popup.js
   Lines 24-33: mol2chemfig DOM elements
   Lines 85-93: Load mol2chemfig settings
   Lines 215-295: Event listeners for mol2chemfig options
   
✅ chem-extension/popup.html
   Lines 350-450: mol2chemfig options UI section
```

---

## Quick Commands

### Restart Everything Fresh
```powershell
# Kill all Python
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Sleep -Seconds 2

# Clear Python cache
cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"
Remove-Item app/__pycache__ -Recurse -Force -ErrorAction SilentlyContinue

# Start server
python run_server.py
```

### Check if Files Have New Code
```powershell
# Backend should return matches
Select-String -Path "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\app\api.py" -Pattern "aromatic_circles.*request.args"

# Extension should return matches
Select-String -Path "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\chem-extension\content.js" -Pattern "URLSearchParams"
```

Both should show line numbers if code is present.

---

## Success Criteria

✅ Console shows: `⚙️ Rendering Options: {aromaticCircles: true, ...}` (matches popup settings)  
✅ Console shows: `API URL: ...&aromatic_circles=true&...` (options in URL)  
✅ Console shows: `Cache URL: ...aromatic_circles_...` (option in filename)  
✅ Visual difference: Benzene with circles vs without  
✅ Cache files: Different filenames for different options in `svg-cache/` folder

---

## If STILL Not Working After All This

1. **Export extension logs**:
   - Open console on debug page
   - Right-click → Save as... → console.log
   - Send me the file

2. **Check chrome.storage**:
   - F12 → Application → Storage → chrome.storage.sync
   - Take screenshot
   - Should see: `aromaticCircles: true`, `rendererEngine: "moleculeviewer"`, etc.

3. **Check server is receiving requests**:
   - Look at Flask terminal output
   - Should see: `SMILES endpoint: c1ccccc1`
   - Should see: `Cache URL: ...`

4. **Verify extension is injecting**:
   - Console should show extension logs
   - If no logs: extension content script not running
   - Check manifest.json `matches` patterns

---

## Current Status

✅ Code changes completed  
✅ Debug logging added  
✅ Test page created  
⏳ **YOUR TURN**: Reload extension → Test → Check console logs

**Next Action**: Open `extension_debug_test.html`, press F12, reload page, screenshot console logs
