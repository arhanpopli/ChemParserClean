# How to Test the Options System

## Step 1: Kill Old Processes
```powershell
# Kill all Python processes
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Sleep -Seconds 3
```

## Step 2: Clear Python Cache
```powershell
cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"
rm -r app/__pycache__ -Force -ErrorAction SilentlyContinue
rm -r __pycache__ -Force -ErrorAction SilentlyContinue
```

## Step 3: Start Fresh Server
```powershell
cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"
python run_server.py
```
Should see:
```
==================================================
  MoleculeViewer Flask Server
==================================================
Server available at:
  http://localhost:5000
```

## Step 4: Verify API Works
In a **SEPARATE** PowerShell terminal (don't use the one running the server):

```powershell
# Test with default options
Invoke-WebRequest -Uri "http://localhost:5000/img/smiles?smiles=c1ccccc1" -UseBasicParsing | ConvertFrom-Json | Select-Object cache_url

# Test with aromatic_circles option
Invoke-WebRequest -Uri "http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true" -UseBasicParsing | ConvertFrom-Json | Select-Object cache_url
```

You should see **different cache URLs**:
- Without options: `smiles_c1ccccc1_....svg`
- With aromatic_circles: `smiles_c1ccccc1_aromatic_circles_....svg`

## Step 5: Test in Browser

### Test Page 1: Quick Visual Test
```
http://localhost:5000/quick_test.html
```
You should see 3 different benzene structures (visually different if options work)

### Test Page 2: Comprehensive Test
```
http://localhost:5000/test_options_cache.html
```
Click "Run Tests" and watch the cache URLs - each should be different

## Step 6: Reload Chrome Extension

**IMPORTANT**: The extension also needs to be reloaded!

1. Open Chrome
2. Go to: `chrome://extensions`
3. Find "Chemistry Formula Renderer" extension
4. Click the **RELOAD** button (circular arrow icon)
5. Or: Press `Ctrl+Shift+R` to reload the extension

## Step 7: Test Extension

1. Open a page with chemical formulas (or create a test page with SMILES)
2. Open extension popup
3. Select **MoleculeViewer** renderer
4. Toggle **"Aromatic Circles"** option ON
5. **Reload the page** containing the formula
6. Open DevTools Console (F12)
7. Look for logs showing:
   ```
   API URL: http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true&...
   Cache URL: http://localhost:5000/cache/smiles_c1ccccc1_aromatic_circles_...
   ```

## What Should Happen

### If Working:
- Cache URLs contain option names: `aromatic_circles`, `show_carbons`, etc.
- Each option combination creates a **different cached file**
- Visual differences in the rendered structures

### If NOT Working:
- All cache URLs look the same (no options in filename)
- All structures render identically
- Check: Is extension reloaded? Is server running with new code?

## Debug Checklist

- [ ] Flask server shows: "Running on 0.0.0.0:5000"
- [ ] API returns different cache_url for each options combination
- [ ] Chrome extension is **RELOADED** (step 6)
- [ ] Browser page is **RELOADED** after changing options
- [ ] DevTools Console shows correct API URLs with options
- [ ] Cache files created in `MoleculeViewer/svg-cache/` with option names

## Quick Command to Test Everything

Run all these in order:
```powershell
# 1. Kill old processes
Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force; sleep 2

# 2. Start fresh server (in background)
Start-Process python -ArgumentList "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\run_server.py" -WorkingDirectory "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"

# 3. Wait for server
sleep 3

# 4. Test API
$url1 = "http://localhost:5000/img/smiles?smiles=c1ccccc1"
$url2 = "http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true"

$r1 = (Invoke-WebRequest -Uri $url1 -UseBasicParsing | ConvertFrom-Json).cache_url
$r2 = (Invoke-WebRequest -Uri $url2 -UseBasicParsing | ConvertFrom-Json).cache_url

Write-Host "Default:            $r1"
Write-Host "Aromatic Circles:   $r2"
Write-Host "Are they different? $(if ($r1 -eq $r2) { 'NO - PROBLEM!' } else { 'YES - WORKING!' })"
```

## Cache File Location
Check if options are being saved:
```powershell
ls "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\svg-cache\" | Sort-Object LastWriteTime -Descending | Select-Object -First 10
```

Filenames should include option names like:
- `smiles_c1ccccc1_abc123.svg` (no options)
- `smiles_c1ccccc1_aromatic_circles_def456.svg` (with aromatic_circles)

## Common Issues

### Issue: "Cannot connect to server"
- Is Flask server running? Check terminal
- Is port 5000 in use? `netstat -an | findstr :5000`

### Issue: Cache URLs all the same
- Flask server loaded old code (kill and restart)
- Check `app/__pycache__` is deleted before restarting

### Issue: Extension not showing options
- Did you reload the extension? (chrome://extensions)
- Is extension code updated? Check `chem-extension/content.js` has new options code

### Issue: Options work on web but not extension
- Extension might be using old cached API responses
- Hard refresh: `Ctrl+Shift+R` on page with formulas
- Clear browser cache: `Ctrl+Shift+Delete`

## Files to Check

Make sure these were modified:
```
✓ MoleculeViewer/app/api.py - Lines 657-680 (img_smiles options)
✓ MoleculeViewer/app/api.py - Lines 722-745 (img_nomenclature options)
✓ chem-extension/content.js - Lines 800-845 (build options query string)
```

## Verify Changes
```powershell
# Check api.py has options parsing
Select-String -Path "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer\app\api.py" -Pattern "aromatic_circles.*request.args.get"

# Check content.js has URLSearchParams
Select-String -Path "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\chem-extension\content.js" -Pattern "URLSearchParams"
```

Both should return results if changes were applied.
