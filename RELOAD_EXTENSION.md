# ⚠️ CRITICAL: You MUST Reload the Chrome Extension

## The Problem
You updated the code in:
- ✅ `MoleculeViewer/app/api.py` - API now accepts options
- ✅ `chem-extension/content.js` - Extension now sends options

BUT the **Chrome extension is still running old code** because it hasn't been reloaded.

## Solution: Reload Extension in Chrome

### Method 1: Using Chrome Extension Manager
1. Open Chrome
2. Type in address bar: `chrome://extensions`
3. Find "Chemistry Formula Renderer"
4. Click the **RELOAD** button (circular arrow icon)

### Method 2: Using DevTools Shortcut
1. Open any page with chemical formulas
2. Press: `Ctrl+Shift+R`  (Hard refresh)
3. Open DevTools: `F12`
4. Check Console for new API URLs

### What to Look For After Reload
Open DevTools Console (F12) and look for messages like:

```
API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true&show_carbons=true&...
Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_show_carbons_...
```

If you see these, it's working!

## Quick Test Workflow

1. **Kill old Flask server**
   ```powershell
   Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force
   sleep 3
   ```

2. **Start fresh Flask server**
   ```powershell
   cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"
   python run_server.py
   ```

3. **Reload Chrome Extension**
   - Go to `chrome://extensions`
   - Click RELOAD on "Chemistry Formula Renderer"

4. **Test Extension**
   - Open extension popup
   - Select MoleculeViewer
   - Enable "Aromatic Circles"
   - Reload the page with formulas
   - Check DevTools Console for new API URLs

## Debug Log Example

### With Options Working:
```
📤 Using SMILES endpoint
API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1&show_carbons=true&aromatic_circles=true&...
✅ Got response from backend
📊 Backend data: {success: true, cache_url: "http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_show_carbons_abc123.svg"}
📍 Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_show_carbons_abc123.svg
✅ Image loaded successfully
```

### WITHOUT Options (Old Code):
```
📤 Using SMILES endpoint
API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1
✅ Got response from backend
📊 Backend data: {success: true, cache_url: "http://localhost:5000/cache/smiles_c1ccccc1_abc123.svg"}
📍 Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_abc123.svg
✅ Image loaded successfully
```

Notice the difference in the URL and cache_url!

## After Reload, You Should See:

### 1. Different Cache URLs for Different Options
```
Default:            http://localhost:5000/cache/smiles_c1ccccc1_xyz.svg
Aromatic Circles:   http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_xyz.svg
Show Carbons:       http://localhost:5000/cache/smiles_c1ccccc1_show_carbons_xyz.svg
Both Options:       http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_show_carbons_xyz.svg
```

### 2. Visual Differences in Structures
- **Aromatic Circles**: Benzene shows a circle inside the ring
- **Show Carbons**: Carbon atoms labeled with "C"
- **Show Methyls**: Methyl groups labeled with "CH3"

### 3. Cache Files in Disk
```powershell
ls "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\svg-cache\" | Sort-Object LastWriteTime -Descending | Select-Object -First 10
```

Should show files with different names encoding the options

## Verify Changes Were Applied

Run these commands to confirm files have new code:

```powershell
# Check backend has options
Select-String -Path "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\app\api.py" -Pattern "aromatic_circles.*request.args.get" | Select-Object -First 1

# Check extension sends options
Select-String -Path "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\chem-extension\content.js" -Pattern "URLSearchParams" | Select-Object -First 1
```

Both should return a result with line numbers.

## Next Steps

1. ✅ Files have been updated (verified above)
2. ⏳ **Reload Chrome Extension** (this is YOUR action)
3. ⏳ Restart Flask server with fresh Python
4. ⏳ Test in browser and check DevTools Console
5. ⏳ Verify cache URLs include option names
6. ⏳ Visual differences in structures

**The code changes are done. NOW you need to reload the extension!**
