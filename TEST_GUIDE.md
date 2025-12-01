# ChemParser Testing Guide

## 🚀 Quick Start

### Start All Servers
```bash
start-servers.bat
```

OR start just MolView servers:
```bash
start-molview.bat
```

## ✅ Step-by-Step Testing

### 1. Verify MolView PHP Server (Port 8000)

Open these URLs in your browser:

**Main viewer:**
```
http://localhost:8000/
```
✅ Should see MolView homepage

**Embed URLs (these are what the search API returns):**
```
http://localhost:8000/embed/v2/?cid=2244
```
✅ Should show Aspirin in 3D viewer

```
http://localhost:8000/embed/v2/?pdbid=1RHV
```
✅ Should show Rhinovirus protein structure

```
http://localhost:8000/embed/v2/?codid=9004137
```
✅ Should show Brookite mineral structure

### 2. Verify Search API (Port 8001)

**Test compound search:**
```
http://localhost:8001/search?q=aspirin
```

✅ Expected response:
```json
{
    "query": "aspirin",
    "corrected_query": null,
    "name": "Aspirin",
    "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "isomeric_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "sdf": {
        "available": true,
        ...
    },
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/2244",
    "primary_type": "compound",
    "cid": 2244,
    "embed_url": "http://localhost:8000/embed/v2/?cid=2244"
}
```

**Test autocorrect (typo):**
```
http://localhost:8001/search?q=booxite
```

✅ Expected response:
```json
{
    "query": "booxite",
    "corrected_query": "Brookite",  ← AUTOCORRECTED!
    "name": "Brookite",
    ...
    "primary_type": "mineral",
    "codid": "9004137",
    "embed_url": "http://localhost:8000/embed/v2/?codid=9004137"
}
```

**Test protein:**
```
http://localhost:8001/search?q=rhinovirus
```

✅ Expected response:
```json
{
    "query": "rhinovirus",
    "name": "Rhinovirus",
    ...
    "primary_type": "biomolecule",
    "pdbid": "1RHV",
    "embed_url": "http://localhost:8000/embed/v2/?pdbid=1RHV"
}
```

### 3. Test Extension

#### Load Extension
1. Open Chrome
2. Go to `chrome://extensions`
3. Enable "Developer mode"
4. Click "Load unpacked"
5. Select: `C:\Users\Kapil\Personal\PROJECTS\Chemparser\chem-extension`

#### Test Cases

**Test 1: Compound (No Typo)**
```
Type in any webpage: chem:aspirin:
```
✅ Expected:
- No autocorrect banner
- Molecule structure displayed
- Console shows: "Query: aspirin" → "No autocorrect needed"

**Test 2: Typo Correction**
```
Type: chem:booxite:
```
✅ Expected:
- Purple banner: "✓ Autocorrected: booxite → Brookite"
- Mineral structure displayed
- Console shows: "🎯 Autocorrect: booxite → Brookite"

**Test 3: Protein**
```
Type: chem:rhinovirus:
```
✅ Expected:
- No autocorrect banner (name is correct)
- Protein structure displayed
- Console shows: "Type: biomolecule"

**Test 4: More Typos**
```
chem:aspriin:   → Should autocorrect to Aspirin
chem:etanol:    → Should autocorrect to Ethanol
chem:benzen:    → Should autocorrect to Benzene
```

**Test 5: SMILES (Direct)**
```
chem:CCO:       → Should work (Ethanol)
chem:CC(=O)OC1=CC=CC=C1C(=O)O:  → Should work (Aspirin)
```

### 4. Test All Rendering Engines

Open extension popup, try each engine:

#### MoleculeViewer (Default)
1. Select "MoleculeViewer" in popup
2. Type: `chem:booxite:`
3. ✅ Should show autocorrect + SVG diagram

#### PubChem
1. Select "PubChem" in popup
2. Reload page
3. Type: `chem:booxite:`
4. ✅ Should show autocorrect + PubChem image

#### mol2chemfig
1. Select "mol2chemfig" in popup
2. Reload page
3. Type: `chem:aspirin:`
4. ✅ Should show LaTeX-style diagram

#### Client-Side
1. Select "Client-Side" in popup
2. Reload page
3. Type: `chem:ethanol:`
4. ✅ Should show offline-rendered SVG

### 5. Console Log Verification

Press F12 to open console, should see:

```
🔍 Universal Search API enabled - autocorrect active for ALL engines
🔍 UNIVERSAL SEARCH API - Preprocessing query
🔎 Query: "booxite"
✅ Search API Result: {...}
🎯 Autocorrect: "booxite" → "Brookite"
🧬 SMILES: ...
📦 Type: mineral
🧪 Using MOLECULEVIEWER renderer engine
```

## 🐛 Troubleshooting

### Problem: "localhost refused to connect"

**Check Server Status:**

1. **Is MolView PHP running?**
   ```
   Open: http://localhost:8000/
   Should see: MolView homepage
   ```

2. **Is Search API running?**
   ```
   Open: http://localhost:8001/search?q=test
   Should see: JSON response
   ```

3. **Check console windows**
   - Should have 5 console windows open after running `start-servers.bat`
   - Look for errors in red text

### Problem: Wrong compound displayed

**Reason:** Search API similarity matching found different result

**Debug:**
```
1. Open: http://localhost:8001/search?q=yourquery
2. Check the "corrected_query" field
3. Check the "primary_type" field (compound/biomolecule/mineral)
```

**Solution:**
- Use exact SMILES: `chem:CCO:`
- Or use more specific name: `chem:ethyl alcohol:` instead of `chem:ethanol:`

### Problem: No autocorrect banner shown

**Expected:** This is normal if the name was spelled correctly!

The banner only appears when:
- `corrected_query` is not null
- Query was actually autocorrected

### Problem: Extension not finding molecules

**Check:**
1. All servers running? (especially ports 8000 and 8001)
2. Extension loaded in Chrome?
3. Try reloading the webpage
4. Check browser console (F12) for errors

## 📊 Performance Benchmarks

**Expected timings:**

| Operation | Time |
|-----------|------|
| Search API query | 50-150ms |
| Autocorrect match | 10-50ms |
| MoleculeViewer render | 200-500ms |
| PubChem image fetch | 300-800ms |
| Total (first load) | 500-1500ms |
| Cached result | 50-200ms |

## 🎓 Advanced Tests

### Test Stereochemistry

1. Enable "Use 3D SMILES" in popup
2. Type: `chem:glucose:`
3. ✅ Should use isomeric SMILES with stereochemistry

### Test Different Compounds

**Minerals:**
- `chem:quartz:`
- `chem:calcite:`
- `chem:halite:`

**Proteins:**
- `chem:insulin:`
- `chem:hemoglobin:`
- `chem:lysozyme:`

**Drugs:**
- `chem:ibuprofen:`
- `chem:paracetamol:`
- `chem:morphine:`

### Test Edge Cases

**Very short query:**
```
chem:H2O:  → Should work (water)
```

**Very long SMILES:**
```
chem:CC(C)CC1=CC=C(C=C1)C(C)C(=O)O:  → Should work (ibuprofen)
```

**Unknown compound:**
```
chem:xyzabc123:  → Should return error or no results
```

## ✅ Success Criteria

All tests pass if:

1. ✅ Both MolView servers start without errors
2. ✅ Search API returns correct JSON for all test queries
3. ✅ Autocorrect works for typos (shows purple banner)
4. ✅ All 4 rendering engines work with autocorrect
5. ✅ Console logs show search API preprocessing
6. ✅ Embed URLs use `/embed/v2/?cid=` format
7. ✅ Minerals, proteins, and compounds all work

## 🎯 Example Test Session

```bash
# 1. Start servers
start-servers.bat

# 2. Test search API
# Open: http://localhost:8001/search?q=booxite
# Verify: corrected_query = "Brookite"

# 3. Test embed URL
# Copy embed_url from step 2
# Open in browser
# Verify: Mineral structure loads

# 4. Test extension
# Load extension in Chrome
# Open any webpage
# Type: chem:booxite:
# Verify: Purple banner + structure

# 5. Test all engines
# Try MoleculeViewer, PubChem, mol2chemfig, Client-Side
# Verify: All show autocorrect + render correctly
```

## 📝 Notes

- **Autocorrect threshold:** 40% similarity (configurable in search-server.js)
- **Cache:** Results cached by browser and search API
- **Port 8001:** REQUIRED for ALL queries (universal preprocessor)
- **Port 8000:** REQUIRED for embed URLs
- **Embed format:** Always `/embed/v2/?cid=X` or `?pdbid=X` or `?codid=X`

---

**Last Updated:** 2025-01-28
**Version:** 3.0 (Universal Search with Fixed Embeds)
