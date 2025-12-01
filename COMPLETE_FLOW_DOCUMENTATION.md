# ChemParser Extension - Complete Flow Documentation

## 🎯 Overview

The ChemParser extension now uses **MolView Search API (port 8001)** as the PRIMARY data source for all molecule queries, with intelligent type detection and autocorrection.

---

## 📊 Complete Data Flow

### Step 1: User Input
User types in webpage: `chem:histamine:`

### Step 2: Extension Pattern Detection
- Content script detects pattern: `/chem:([^:]+):/g`
- Extracts nomenclature: `"histamine"`
- Creates placeholder `<img>` element

### Step 3: Query MolView Search API
```javascript
const searchUrl = `http://localhost:8001/search?q=histamine`;
const searchResult = await backgroundFetchJSON(searchUrl);
```

### Step 4: MolView Search Response
```json
{
  "query": "histamine",
  "corrected_query": "Histamine",
  "name": "Histamine",
  "canonical_smiles": "C1=C(NC=N1)CCN",
  "isomeric_smiles": "C1=C(NC=N1)CCN",
  "sdf": {
    "available": true,
    "size_bytes": 3031,
    "format": "SDF",
    "download_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/774/SDF"
  },
  "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/774",
  "primary_type": "compound",
  "embed_url": "http://localhost:8000/embed/v2/?cid=774",
  "cid": 774
}
```

### Step 5A: Compound Rendering (Client-Side Mode)
**If** `settings.rendererEngine === 'client-side'`:

1. Extension uses `canonical_smiles` from search result
2. Passes SMILES to `SmilesDrawer` library
3. Generates SVG directly in browser
4. Replaces `<img>` with rendered SVG

**Console Output**:
```
🌐 [Client] Priority 0: Using MolView Search API (PRIMARY)
✅ [Client] MolView Search API SUCCESS: C1=C(NC=N1)CCN
📊 [Client] MolView Data: {corrected: "Histamine", type: "compound", has_sdf: true}
🎨 Rendering with SmilesDrawer:
📊 SMILES: C1=C(NC=N1)CCN
✅ SVG rendered as img with standard wrapper
```

### Step 5B: Compound Rendering (Server-Side Mode)
**If** `settings.rendererEngine === 'moleculeviewer'` (or other server):

1. Extension uses `smilesBridge()` function
2. `smilesBridge()` queries MolView Search API (Priority 0)
3. Returns SMILES to renderer
4. Sends SMILES to MoleculeViewer server (`http://localhost:5000/img/smiles?smiles=...`)
5. Replaces `<img>` with server-rendered SVG

**Console Output**:
```
🌉 SMILES BRIDGE: Converting name→SMILES histamine
🌐 [Bridge] Priority 0: Trying MolView Search API...
✅ [Bridge] MolView Search SUCCESS: C1=C(NC=N1)CCN
📊 MolView Data: {corrected: "Histamine", type: "compound", has_sdf: true}
✏️ Autocorrected: histamine → Histamine
```

---

## 🧬 Protein Flow (e.g., `chem:rhinovirus:`)

### Step 1-3: Same as above

### Step 4: MolView Search Response (Protein)
```json
{
  "query": "rhinovirus",
  "corrected_query": "Rhinovirus",
  "name": "Rhinovirus",
  "canonical_smiles": null,
  "isomeric_smiles": null,
  "primary_type": "biomolecule",
  "pdbid": "4rhv",
  "source_url": "https://www.rcsb.org/structure/4rhv",
  "embed_url": "http://localhost:8000/embed/v2/?pdbid=4rhv",
  "pdb_url": "https://files.rcsb.org/view/4rhv.pdb"
}
```

### Step 5: Protein Rendering (ALL Modes)
**When** `primary_type === 'biomolecule'`:

1. Extension detects protein via `bridgeResult.is_protein === true`
2. **Skips SMILES rendering** (no SMILES available)
3. Creates `<iframe>` element
4. Sets `iframe.src = searchResult.embed_url` → `http://localhost:8000/embed/v2/?pdbid=4rhv`
5. Adds metadata: `Type: biomolecule | PDB ID: 4rhv`
6. Replaces `<img>` with `<iframe>` + metadata

**Console Output**:
```
🌉 SMILES BRIDGE: Converting name→SMILES rhinovirus
🌐 [Bridge] Priority 0: Trying MolView Search API...
⚠️ [Bridge] Detected biomolecule/protein - skipping SMILES conversion
🧬 Detected protein/biomolecule - using MolView embed
```

**User sees**: Interactive 3D protein viewer in webpage

---

## 💎 Mineral Flow (e.g., `chem:quartz:`)

### Step 1-3: Same as above

### Step 4: MolView Search Response (Mineral)
```json
{
  "query": "quartz",
  "corrected_query": "Quartz",
  "name": "Quartz",
  "canonical_smiles": null,
  "isomeric_smiles": null,
  "primary_type": "mineral",
  "codid": "5000035",
  "source_url": "http://www.crystallography.net/cod/5000035.html",
  "embed_url": "http://localhost:8000/embed/v2/?codid=5000035",
  "cif_url": "http://localhost:8000/php/cif.php?codid=5000035"
}
```

### Step 5: Mineral Rendering (ALL Modes)
**When** `primary_type === 'mineral'`:

1. Extension detects mineral via `bridgeResult.is_mineral === true`
2. **Skips SMILES rendering** (no SMILES for minerals)
3. Creates `<iframe>` element
4. Sets `iframe.src = searchResult.embed_url` → `http://localhost:8000/embed/v2/?codid=5000035`
5. Adds metadata: `Type: mineral | COD ID: 5000035`
6. Replaces `<img>` with `<iframe>` + metadata

**Console Output**:
```
🌉 SMILES BRIDGE: Converting name→SMILES quartz
🌐 [Bridge] Priority 0: Trying MolView Search API...
⚠️ [Bridge] Detected mineral - skipping SMILES conversion
💎 Detected mineral - using MolView embed
```

**User sees**: Interactive crystal structure viewer in webpage

---

## 🔄 Fallback Chain

### Priority 0: MolView Search API (PRIMARY)
- **URL**: `http://localhost:8001/search?q=NOMENCLATURE`
- **Returns**: SMILES, type, autocorrection, embed URL
- **Used**: ALWAYS (for all renderer engines)
- **Fallback**: If fails → Priority 1

### Priority 1: Direct PubChem API (FALLBACK)
- **URL**: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/NOMENCLATURE/property/CanonicalSMILES,IsomericSMILES/JSON`
- **Returns**: SMILES only (no autocorrection or type detection)
- **Used**: Only if MolView Search fails
- **Note**: Cannot detect proteins or minerals

### Priority 2: PubChem Autocomplete (FALLBACK)
- **URL**: `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/NOMENCLATURE/json?limit=5`
- **Returns**: Similar compound names
- **Used**: Only if Priority 1 fails
- **Example**: "sphingomyelin" → finds "Sphingomyelin 16:0"

---

## 🎛️ Renderer Engine Modes

### Client-Side Mode (`settings.rendererEngine === 'client-side'`)
1. Get SMILES from MolView Search API
2. Render SVG using `SmilesDrawer` library (browser-side)
3. **No server calls** for rendering (only for data lookup)
4. **Fastest** option (no network delay for rendering)

### Server-Side Modes (`moleculeviewer`, `mol2chemfig`, `pubchem`)
1. Get SMILES from MolView Search API
2. Send SMILES to rendering server:
   - MoleculeViewer: `http://localhost:5000/img/smiles?smiles=...`
   - mol2chemfig: `http://localhost:5001/m2cf/submit`
   - PubChem: `http://localhost:5002/pubchem/img/...`
3. Server returns SVG
4. Extension displays SVG

---

## 🔍 Code Locations

### MolView Search Server
- **File**: `Molview/molview/search-server.js`
- **Port**: 8001
- **Endpoint**: `/search?q=QUERY`
- **Fixed**: Lines 292-315 (SMILES property fetching)

### Extension - SMILES Bridge
- **File**: `chem-extension/content.js`
- **Function**: `smilesBridge(name, options)` (Lines 384-500)
- **Priority 0**: Lines 406-445 (MolView Search API)
- **Protein Detection**: Lines 415-424
- **Mineral Detection**: Lines 427-436

### Extension - Client-Side Rendering
- **File**: `chem-extension/content.js`
- **Function**: `renderClientSide(moleculeData, img)` (Lines 2208-2625)
- **MolView Search**: Lines 2257-2279 (Priority 0)
- **PubChem Fallback**: Lines 2281-2320 (Priority 1)

### Extension - Server-Side Rendering
- **File**: `chem-extension/content.js`
- **Function**: `loadMoleculeViewerImage(img)` (Lines 2662-2900)
- **Protein/Mineral Handling**: Lines 2757-2827
- **Uses**: `smilesBridge()` at line 2755

---

## 📋 Summary Table

| Input | MolView Type | SMILES? | Renderer | Output |
|-------|-------------|---------|----------|--------|
| `chem:histamine:` | compound | Yes | Client-Side | SVG (SmilesDrawer) |
| `chem:histamine:` | compound | Yes | MoleculeViewer | SVG (server) |
| `chem:rhinovirus:` | biomolecule | No | ANY | iframe (MolView 3D) |
| `chem:quartz:` | mineral | No | ANY | iframe (Crystal viewer) |
| `chem:C1=CC=CC=C1:` | N/A | Yes | Client-Side | SVG (SMILES input) |

---

## ✅ Key Benefits

1. **Single API Call**: One query to MolView Search gets:
   - SMILES data
   - Molecule type (compound/protein/mineral)
   - Autocorrection
   - Embed URL
   - SDF/PDB/CIF data

2. **Intelligent Routing**:
   - Compounds → 2D/3D structure rendering
   - Proteins → 3D protein viewer
   - Minerals → Crystal structure viewer

3. **Autocorrection**:
   - "histamine" → "Histamine" ✅
   - "cafeine" → "caffeine" ✅
   - Fuzzy matching for typos

4. **No More Wrong Results**:
   - **Before**: `chem:rhinovirus:` → random PubChem compound ❌
   - **After**: `chem:rhinovirus:` → 3D protein structure ✅

---

## 🚀 Testing

Open `test_molview_integration.html` in browser with extension enabled.

Test cases include:
- ✅ Regular compounds
- ✅ Autocorrection
- ✅ Proteins
- ✅ Minerals
- ✅ SMILES input
- ✅ Edge cases

---

**Date**: 2025-11-28
**Status**: ✅ Complete and Working
**Integration**: MolView Search API → ChemParser Extension → User
