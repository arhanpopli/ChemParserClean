# Console Debug Logs Guide

## 🔍 How to Debug

When you type `chem:booxite:` and press F12 in Chrome, you should see this exact sequence of logs:

## ✅ Expected Log Sequence

### 1. Extension Initialization
```
[ChemRenderer] [INFO] 🔍 Universal Search API: Port 8001 (autocorrect & intelligent filtering enabled for ALL engines!)
[ChemRenderer] [SUCCESS] ✅ Settings loaded
[ChemRenderer] [INFO] Renderer Engine: 🧪 MoleculeViewer (localhost:5000)
```

### 2. Pattern Detection
```
🔍 Universal Search API enabled - autocorrect active for ALL engines
🧪 Using MOLECULEVIEWER renderer engine
🧪 LOADMOLECULEVIEWERIMAGE CALLED!
```

### 3. Search API Query (THE IMPORTANT PART!)
```
═══════════════════════════════════════════════════════════
🔍 UNIVERSAL SEARCH API - Preprocessing query
═══════════════════════════════════════════════════════════
🔎 Query String: "booxite"
📡 Calling Search API:
   URL: http://localhost:8001/search?q=booxite&format=compact
   Method: GET
```

### 4. Search API Response
```
📥 Search API Response Received:
{
  "query": "booxite",
  "corrected_query": "Brookite",
  "name": "Brookite",
  "canonical_smiles": null,
  "isomeric_smiles": null,
  "sdf": {
    "available": true,
    "size_bytes": 230,
    ...
  },
  "source_url": "http://www.crystallography.net/cod/9004137.html",
  "primary_type": "mineral",
  "codid": "9004137",
  "embed_url": "http://localhost:8000/embed/v2/?codid=9004137"   ← THIS IS KEY!
}
```

### 5. Processed Data
```
📊 Processed Data:
   Original Query: "booxite"
   Corrected Name: "Brookite"
   Was Corrected: true
   SMILES: N/A
   Type: mineral
   CID: N/A
   PDBID: N/A
   CODID: 9004137
   🔗 EMBED URL: http://localhost:8000/embed/v2/?codid=9004137   ← VERIFY THIS!
```

### 6. Autocorrect Banner
```
🎯 AUTOCORRECT: "booxite" → "Brookite"
```

### 7. Data Returned
```
📋 Data returned from Search API:
   - Corrected Name: Brookite
   - SMILES: null
   - Compound Type: mineral
   - Was Autocorrected: true
🎨 Showing autocorrect banner in UI...
```

### 8. Rendering Step
```
✅ STEP 1 COMPLETE - Search API data processed
🎨 STEP 2: Rendering with MoleculeViewer...
📤 Using SMILES endpoint for MoleculeViewer
   SMILES: [whatever SMILES was found]
🌐 MoleculeViewer API URL:
   http://localhost:5000/img/smiles?smiles=XXX&width=300&height=200&json=true&t=...
```

## 🐛 Debugging Issues

### Issue 1: "localhost refused to connect"

**Look for this in console:**
```
❌ Error: Failed to fetch from http://localhost:8001/search?q=booxite
```

**OR:**
```
❌ Search API returned error: [error message]
```

**Diagnosis:**
- Search API (port 8001) is NOT running
- Solution: Run `start-servers.bat` or `start-molview.bat`

### Issue 2: Wrong embed URL

**Look for:**
```
🔗 EMBED URL: http://localhost:8000/?codid=9004137   ← WRONG! Missing /embed/v2/
```

**Should be:**
```
🔗 EMBED URL: http://localhost:8000/embed/v2/?codid=9004137   ← CORRECT!
```

**Diagnosis:**
- search-server.js has wrong `localMolViewUrl`
- Solution: Check line 250 in search-server.js:
  ```javascript
  const localMolViewUrl = "http://localhost:8000/embed/v2/";  // Must end with /embed/v2/
  ```

### Issue 3: No autocorrect

**Look for:**
```
📊 Processed Data:
   Was Corrected: false   ← No typo detected
```

**Diagnosis:**
- Name was spelled correctly OR
- Search API couldn't find a match
- This is NORMAL if query is correct

### Issue 4: Search API not being called

**Missing logs:**
- No "🔍 UNIVERSAL SEARCH API" logs
- No "📡 Calling Search API" logs

**Diagnosis:**
- Extension not loaded OR
- Rendering engine set to "client-side" (which skips search API)
- Solution: Select MoleculeViewer/PubChem/mol2chemfig engine

## 📊 Quick Diagnosis Checklist

Copy/paste this to test each step:

```
# 1. Is extension loaded?
Look for: [ChemRenderer] [INFO] in console

# 2. Is Search API being called?
Look for: 📡 Calling Search API:
          URL: http://localhost:8001/search?q=...

# 3. Is Search API responding?
Look for: 📥 Search API Response Received:
          {JSON data}

# 4. Is embed URL correct?
Look for: 🔗 EMBED URL: http://localhost:8000/embed/v2/?...
                                                ↑↑↑↑↑↑↑↑
                                         Must have /embed/v2/

# 5. Is autocorrect working?
Look for: 🎯 AUTOCORRECT: "booxite" → "Brookite"

# 6. Is MoleculeViewer being called?
Look for: 🌐 MoleculeViewer API URL:
          http://localhost:5000/img/smiles?...
```

## 🎯 Test Commands

### Test 1: Search API Directly
Open in browser:
```
http://localhost:8001/search?q=booxite
```

Should see JSON with:
```json
{
  "corrected_query": "Brookite",
  "embed_url": "http://localhost:8000/embed/v2/?codid=9004137"
}
```

### Test 2: Embed URL
Copy the `embed_url` from Test 1, paste in browser.

Should see: 3D mineral viewer

### Test 3: Extension
Type in any webpage:
```
chem:booxite:
```

**Expected console sequence:**
1. ✅ Pattern detected
2. ✅ Search API called
3. ✅ Response received with embed URL
4. ✅ Autocorrect banner shown
5. ✅ MoleculeViewer called with SMILES

## 🔥 Common Error Messages

### "Failed to fetch"
```
❌ Error: Failed to fetch from http://localhost:8001/...
```
**Fix:** Start Search API (port 8001)

### "CORS error"
```
❌ Access to fetch at 'http://localhost:8001' ... has been blocked by CORS
```
**Fix:** Search API should have CORS headers. Check search-server.js line 351:
```javascript
res.setHeader('Access-Control-Allow-Origin', '*');
```

### "No results found"
```
📥 Search API Response Received:
{
  "error": "No results found",
  "query": "xyzabc"
}
```
**Fix:** This is normal for unknown compounds. Try a real compound name.

### "Autocomplete only gives a name"
This is not an error - it's a comment in the code. Ignore it.

## 📸 Screenshot What to Look For

When debugging, take a screenshot showing:

1. **The query you typed:** `chem:booxite:`
2. **Console with these visible:**
   - 🔍 UNIVERSAL SEARCH API header
   - 📡 Calling Search API URL
   - 📥 Search API Response with JSON
   - 🔗 EMBED URL line (most important!)
3. **Any error messages in red**

Share this screenshot for debugging help.

## ✅ Success Indicators

You know it's working when you see:

1. ✅ Purple autocorrect banner in webpage
2. ✅ Console shows: `🎯 AUTOCORRECT: "booxite" → "Brookite"`
3. ✅ Embed URL contains `/embed/v2/`
4. ✅ Molecule structure displays

## 🚀 Full Example Log (Working)

```
🔍 Universal Search API enabled - autocorrect active for ALL engines
🧪 Using MOLECULEVIEWER renderer engine
═══════════════════════════════════════════════════════════
🔍 UNIVERSAL SEARCH API - Preprocessing query
═══════════════════════════════════════════════════════════
🔎 Query String: "booxite"
📡 Calling Search API:
   URL: http://localhost:8001/search?q=booxite&format=compact
   Method: GET
📥 Search API Response Received:
{
  "query": "booxite",
  "corrected_query": "Brookite",
  "embed_url": "http://localhost:8000/embed/v2/?codid=9004137",
  ...
}
📊 Processed Data:
   🔗 EMBED URL: http://localhost:8000/embed/v2/?codid=9004137
🎯 AUTOCORRECT: "booxite" → "Brookite"
✅ STEP 1 COMPLETE - Search API data processed
🎨 STEP 2: Rendering with MoleculeViewer...
📤 Using SMILES endpoint for MoleculeViewer
🌐 MoleculeViewer API URL: http://localhost:5000/img/smiles?...
✅ MoleculeViewer image loaded successfully
```

---

**Last Updated:** 2025-01-28
**Purpose:** Debugging "localhost refused to connect" issues
