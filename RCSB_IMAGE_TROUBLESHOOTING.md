# RCSB Image Display Fix

## Problem
The extension console shows:
```
[ChemRenderer] Background fetch FAILED: Failed to fetch
❌ Error: Failed to fetch
```

This means the extension cannot connect to the search server on localhost:8001.

## Root Causes

1. **Search server not running** - The server must be started before loading pages with chemistry formulas
2. **CORS/CSP blocking** - Some websites (like ChatGPT) may block localhost connections
3. **Background script not receiving messages** - Chrome extension messaging might fail

## Solution

### 1. Start the Search Server FIRST

**Always start the server before opening pages with chemistry formulas:**

```powershell
cd C:\Users\Kapil\Personal\STUFF\Chemparser
node Molview\molview\search-server.js
```

Or use the batch file if you have one:
```powershell
.\1-start-all.bat
```

### 2. Verify Server is Running

Test the endpoint:
```powershell
Invoke-RestMethod -Uri "http://localhost:8001/search?q=hemoglobin" -Method Get
```

Should return JSON with:
- `pdbid`: "4hhb"
- `image_url`: "https://cdn.rcsb.org/images/structures/4hhb_model-1.jpeg"
- `embed_url`: "https://molview.org/?pdbid=4hhb"

### 3. Test Locally First

Open the test file in Chrome:
```
file:///C:/Users/Kapil/Personal/STUFF/Chemparser/test-rcsb-images.html
```

This eliminates CSP/CORS issues from external sites.

### 4. Check Extension Permissions

The manifest.json must have:
```json
"host_permissions": ["<all_urls>"],
"content_security_policy": {
  "extension_pages": "... connect-src 'self' http://localhost:* ..."
}
```

✅ This is already configured correctly in your extension.

## How RCSB Images Work

1. **User types:** `chem:hemoglobin/+3d:`
2. **Extension detects** the formula
3. **Background script fetches** from: `http://localhost:8001/search?q=hemoglobin&format=compact`
4. **Server returns:**
   ```json
   {
     "name": "Hemoglobin",
     "pdbid": "4hhb",
     "image_url": "https://cdn.rcsb.org/images/structures/4hhb_model-1.jpeg",
     "embed_url": "https://molview.org/?pdbid=4hhb"
   }
   ```
5. **Extension creates:**
   - An `<img>` with `src="https://cdn.rcsb.org/images/structures/4hhb_model-1.jpeg"`
   - A "View 3D" button that opens `https://molview.org/?pdbid=4hhb`

## Debugging Steps

### Step 1: Check if server is running
```powershell
Get-Process -Name "node" | Select-Object Id, ProcessName
```

### Step 2: Check if port 8001 is listening
```powershell
netstat -ano | findstr ":8001"
```

### Step 3: Test the endpoint
```powershell
curl http://localhost:8001/search?q=hemoglobin
```

### Step 4: Check Chrome Extension Console

1. Open Chrome DevTools (F12)
2. Go to Console tab
3. Look for these messages:
   - ✅ `[Background] FETCH_API success` - Server is accessible
   - ❌ `[Background] FETCH_API error: Failed to fetch` - Server not accessible
   - ❌ `Background fetch FAILED` - Extension cannot connect

### Step 5: Check Background Script Console

1. Go to `chrome://extensions/`
2. Enable "Developer mode"
3. Click "Service worker" under your extension
4. Check for errors

## Expected Behavior

When working correctly, you should see:

1. **Console logs:**
   ```
   [ChemRenderer] Background fetch SUCCESS
   [Background] FETCH_API success
   ```

2. **Visual result:**
   - A molecule viewer image appears
   - It shows the RCSB protein structure image
   - A "View 3D" button is visible
   - Clicking the button opens MolView with the 3D structure

## Testing Examples

Try these formulas in the test HTML file:

- `chem:hemoglobin/+3d:` → Shows red blood cell protein
- `chem:insulin/+3d:` → Shows insulin structure
- `chem:lysozyme/+3d:` → Shows enzyme structure
- `chem:myoglobin/+3d:` → Shows muscle protein
- `chem:dna/+3d:` → Shows DNA double helix

## Common Issues

### Issue: "Failed to fetch"
**Solution:** Start the search server first

### Issue: "Connection refused"
**Solution:** Server crashed or stopped - restart it

### Issue: No images appear
**Solution:** Check if the `/+3d:` flag is present - this triggers the 3D/image mode

### Issue: Works locally but not on ChatGPT
**Solution:** ChatGPT's CSP might block localhost - this is a browser security feature and expected behavior

## Server URLs Changed

The search server now returns official MolView URLs:

| Type | Old URL | New URL |
|------|---------|---------|
| Proteins | `http://localhost:8000/embed/v2/?pdbid=4hhb` | `https://molview.org/?pdbid=4hhb` |
| Minerals | `http://localhost:8000/?codid=1010928` | `https://embed.molview.org/v1/?codid=1010928` |
| Chemicals | `http://localhost:8000/embed/v2/?cid=16617` | `https://embed.molview.org/v1/?cid=16617` |

This means you no longer need to run the MolView PHP server on port 8000!

## Quick Start Checklist

- [ ] Search server running on port 8001
- [ ] Extension installed and enabled  
- [ ] Test HTML file opens without errors
- [ ] Console shows "Background fetch SUCCESS"
- [ ] RCSB images appear in the test page
- [ ] "View 3D" buttons work

If all checkmarks are complete, the extension is working correctly!
