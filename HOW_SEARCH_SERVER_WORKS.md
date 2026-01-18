# How the Search Server Works

## Architecture Overview

The search server (`search-server.js`) is a **hybrid search system** that combines local datasets with live API calls to provide comprehensive chemical data.

```
User Query
    ↓
[Local Search] + [API Lookups] → [Processed Result]
```

---

## Step 1: Local Dataset Search

### Local Datasets (Pre-loaded at Startup)

The server loads two JavaScript files at startup:

```javascript
// Line 151-160
const pdbCode = fs.readFileSync(path.join(DATASETS_DIR, 'PDBNames.js'), 'utf8');
const mineralCode = fs.readFileSync(path.join(DATASETS_DIR, 'MineralNames.js'), 'utf8');
```

**PDBNames.js** contains:
- ~180,000 PDB structures (proteins, biomolecules)
- Fields: `name`, `pdbids`, `label`
- Examples: "Hemoglobin", "Insulin", "DNA Polymerase"
- **No** actual 3D coordinates - just names and PDB IDs

**MineralNames.js** contains:
- ~5,000 mineral names from Crystallography Open Database (COD)
- Fields: `name`, `codid` (COD ID)
- Examples: "Calcite", "Quartz", "Pyrite"
- **No** actual crystallographic data - just names and COD IDs

### Local Search Process

```javascript
// Line 242-244
const macromolecules = new AutocompleteBuilder(PDBNames.macromolecules, "name");
const minerals = new AutocompleteBuilder(MineralNames.records, "name");

const MIN_SIM = 40;
const MAX_NUMBER = 10;

let mix = macromolecules.sort(query, MIN_SIM, MAX_NUMBER)
    .concat(minerals.sort(query, MIN_SIM, MAX_NUMBER));
```

This searches the local datasets using **similarity matching** (like autocomplete):
- Takes user query (e.g., "hemoglobin")
- Matches against 180k PDB names + 5k mineral names
- Returns top 10 matches with similarity scores

**Result: Fast, instant suggestions (no network call)**

---

## Step 2: Direct PubChem Lookup

Before using the best local match, the server checks if the **exact query** is a valid PubChem compound:

```javascript
// Line 256-266
const directLookupUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(query)}/cids/JSON`;
try {
    const directDataStr = await fetchUrl(directLookupUrl);
    if (!directDataStr.includes('PUGREST.NotFound')) {
        const directData = JSON.parse(directDataStr);
        if (directData.IdentifierList.CID.length > 0) {
            // Query is a valid compound - add it with HIGH priority
            console.log(`✅ Direct PubChem lookup found: "${query}" is a valid compound (CID: ${...})`);
```

**Example:**
- User types: "cyclohexane"
- Server queries: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/cyclohexane/cids/JSON`
- **Response:** `{ IdentifierList: { CID: [8450] } }`
- **Added to results** with highest priority (similarity = 300)

**Why this step?** PubChem autocomplete returns variants (like "cyclohexanemethanol") but not the exact name "cyclohexane". This ensures exact matches are always preferred.

---

## Step 3: PubChem Autocomplete Suggestions

```javascript
// Line 279-293
const autocpUrl = `https://pubchem.ncbi.nlm.nih.gov/pcautocp/pcautocp.cgi?dict=pc_compoundnames&n=10&q=${encodeURIComponent(query)}`;

try {
    const autocpDataStr = await fetchUrl(autocpUrl);
    const autocpData = JSON.parse(autocpDataStr);

    const pubchemResults = [];
    if (autocpData.autocp_array) {
        for (let name of autocpData.autocp_array) {
            pubchemResults.push({
                name: name,
                label: ucfirst(humanize(name)),
                PubChem_name: name,
                source_type: 'pubchem'
            });
        }
    }

    mix = mix.concat(pubchemResults);
```

This fetches autocomplete suggestions from PubChem, adds them to the local results.

**Example:**
- User types: "aspi"
- PubChem returns: ["aspirin", "aspidosperma", "aspidium"]
- **These are added** to the local search results
- All are re-sorted by similarity

---

## Step 4: Find Best Match & Route

```javascript
// Line 331-338
let bestMatch = sortedMix[0];  // Top result
const queryLower = query.toLowerCase();

for (let i = 0; i < sortedMix.length; i++) {
    if (sortedMix[i].name && sortedMix[i].name.toLowerCase() === queryLower) {
        bestMatch = sortedMix[i];
        console.log(`✅ Found exact match: ${bestMatch.name}`);
        break;
    }
}
```

After sorting, the server checks if it found an exact match and routes based on type:

```
IF bestMatch has `pdbids` field
    → Route to RCSB (Protein)
ELSE IF bestMatch has `codid` field
    → Route to COD (Mineral)
ELSE
    → Route to PubChem (Chemical Compound)
```

---

## Step 5: Fetch Data from APIs

### For Proteins (PDB):

```javascript
// Line 366-385
if (bestMatch.pdbids) {
    const pdbId = bestMatch.pdbids[0];
    result.primary_type = 'biomolecule';
    result.pdbid = pdbId;
    result.source_url = `https://www.rcsb.org/structure/${pdbId}`;
    result.pdb_url = `https://files.rcsb.org/view/${pdbId}.pdb`;
    result.embed_url = `https://molview.org/?pdbid=${pdbId}`;
    
    // Add RCSB image URLs
    result.image_url = `https://cdn.rcsb.org/images/structures/${pdbId.toLowerCase()}_model-1.jpeg`;
    result.assembly_image_url = `https://cdn.rcsb.org/images/structures/${pdbId.toLowerCase()}_assembly-1.jpeg`;
    
    // Fetch PDB file
    try {
        const pdbContent = await fetchUrl(result.pdb_url);
        result.sdf = pdbContent;
    }
}
```

**What it fetches from RCSB:**
- ✅ 3D coordinates (PDB file)
- ✅ Structure images (PNG from CDN)
- ✅ Metadata (protein name, resolution, etc.)

---

### For Minerals (COD):

```javascript
// Line 387-430
} else if (bestMatch.codid) {
    const codid = bestMatch.codid;
    result.primary_type = 'mineral';
    result.codid = codid;
    result.source_url = `http://www.crystallography.net/cod/${codid}.html`;
    result.cif_url = `https://www.crystallography.net/cod/${codid}.cif`;
    result.embed_url = `https://embed.molview.org/v1/?codid=${codid}`;
    
    // Fetch CIF file and parse it
    try {
        const cifContent = await fetchUrl(result.cif_url);
        if (!cifContent.includes('404')) {
            result.sdf = cifContent;  // CIF is treated as SDF
            
            // Parse CIF to extract chemical data
            const cifData = parseCIFData(cifContent);
            result.chemical_name = cifData.chemical_name;
            result.formula = cifData.formula;
            result.canonical_smiles = cifData.canonical_smiles;
        }
    }
}
```

**What it fetches from COD:**
- ✅ Crystallographic coordinates (CIF file)
- ✅ Chemical formula
- ✅ Crystal structure metadata
- ⚠️ SMILES (if available in CIF, otherwise fallback to PubChem)

---

### For Chemicals (PubChem):

```javascript
// Line 434-480
} else {
    result.primary_type = 'compound';
    let cid = null;
    
    // Get CID from compound name
    const cidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(bestMatch.name)}/cids/JSON`;
    try {
        const cidDataStr = await fetchUrl(cidUrl);
        const cidData = JSON.parse(cidDataStr);
        cid = cidData.IdentifierList.CID[0];
        
        // Get SMILES, formula, IUPAC name
        const propsUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/Title,SMILES,ConnectivitySMILES,MolecularFormula,IUPACName/JSON`;
        try {
            const propsDataStr = await fetchUrl(propsUrl);
            const props = propsData.PropertyTable.Properties[0];
            result.isomeric_smiles = props.SMILES;
            result.canonical_smiles = props.ConnectivitySMILES;
            result.formula = props.MolecularFormula;
        }
    }
    
    if (cid) {
        result.embed_url = `https://embed.molview.org/v1/?cid=${cid}`;
        result.sdf_url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF`;
        
        // Fetch 3D SDF
        const sdf3dUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`;
        try {
            const sdfContent = await fetchUrl(sdf3dUrl);
            result.sdf = sdfContent;  // 3D structure coordinates
        }
    }
}
```

**What it fetches from PubChem:**
- ✅ CID (Compound ID)
- ✅ SMILES (structure notation)
- ✅ Chemical formula
- ✅ 3D coordinates (SDF file)
- ✅ 2D coordinates (if 3D not available)

---

## Complete Data Flow Example

### User searches for "hemoglobin":

```
1. LOCAL SEARCH
   - Search local PDBNames.js
   - Find: { name: "Hemoglobin", pdbids: ["4hhb"] }

2. DIRECT PUBCHEM LOOKUP
   - Query: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/hemoglobin/cids/JSON
   - Response: { CID: [5460609] } ← Found as chemical too!

3. PUBCHEM AUTOCOMPLETE
   - Query: https://pubchem.ncbi.nlm.nih.gov/pcautocp/pcautocp.cgi?q=hemoglobin
   - Response: ["Hemoglobin", "Hemoglobin A1C", "Deoxygenated Hemoglobin"]

4. BEST MATCH
   - Mix contains: PDB "Hemoglobin" + PubChem results
   - Exact match found: "Hemoglobin" (from PDB with pdbids)
   - ✅ Route to RCSB (because it has pdbids field)

5. FETCH FROM RCSB
   - PDB ID: 4hhb
   - Fetch: https://files.rcsb.org/view/4hhb.pdb (3D coordinates)
   - Fetch: https://cdn.rcsb.org/images/structures/4hhb_model-1.jpeg (image)

6. RESULT
{
  name: "Hemoglobin",
  primary_type: "biomolecule",
  pdbid: "4hhb",
  source_url: "https://www.rcsb.org/structure/4hhb",
  embed_url: "https://molview.org/?pdbid=4hhb",
  image_url: "https://cdn.rcsb.org/images/structures/4hhb_model-1.jpeg",
  sdf: "<PDB file content with 3D coordinates>",
  _info: { format: "compact" }
}
```

---

### User searches for "calcite":

```
1. LOCAL SEARCH
   - Search local MineralNames.js
   - Find: { name: "Calcite", codid: "1010928" }

2. DIRECT PUBCHEM LOOKUP
   - Query: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/calcite/cids/JSON
   - Response: "PUGREST.NotFound" (Calcite is not a chemical compound)

3. BEST MATCH
   - Exact match found: "Calcite" (from local dataset with codid)
   - ✅ Route to COD (because it has codid field)

4. FETCH FROM COD
   - Fetch: https://www.crystallography.net/cod/1010928.cif (crystal structure)
   - Parse CIF for: chemical_name, formula, SMILES (if available)
   - If no SMILES in CIF:
     - Fallback: Query PubChem for "calcium carbonate"
     - Get SMILES from PubChem

5. RESULT
{
  name: "Calcite",
  formula: "CaCO3",
  primary_type: "mineral",
  codid: "1010928",
  source_url: "http://www.crystallography.net/cod/1010928.html",
  embed_url: "https://embed.molview.org/v1/?codid=1010928",
  cif_url: "https://www.crystallography.net/cod/1010928.cif",
  sdf: "<CIF file content with crystal coordinates>",
  canonical_smiles: "[Ca+2].[O-]C(=O)[O-]"
}
```

---

### User searches for "caffeine":

```
1. LOCAL SEARCH
   - Search PDBNames.js: no match
   - Search MineralNames.js: no match

2. DIRECT PUBCHEM LOOKUP
   - Query: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/caffeine/cids/JSON
   - Response: { CID: [2519] } ✅ FOUND

3. BEST MATCH
   - No local match, but direct PubChem found CID
   - ✅ Route to PubChem (no pdbids, no codid)

4. FETCH FROM PUBCHEM
   - Get compound properties: SMILES, formula, name
   - Fetch: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2519/SDF?record_type=3d
   - Result: 3D coordinates in SDF format

5. RESULT
{
  name: "Caffeine",
  chemical_name: "1,3,7-Trimethylxanthine",
  cid: 2519,
  formula: "C8H10N4O2",
  isomeric_smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  canonical_smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  embed_url: "https://embed.molview.org/v1/?cid=2519",
  sdf: "<3D structure coordinates>"
}
```

---

## Summary

| Data Type | Local Dataset | Live API Calls |
|-----------|---------------|----------------|
| **Proteins** | ✅ PDBNames.js (names + PDB IDs) | ✅ RCSB (3D coords, images) |
| **Minerals** | ✅ MineralNames.js (names + COD IDs) | ✅ COD (crystal structures) |
| **Chemicals** | ❌ None | ✅ PubChem (3D coords, SMILES, formula) |

### Network Calls Summary:

**Always made:**
1. `pubchem.ncbi.nlm.nih.gov` - Direct lookup + autocomplete
2. `rcsb.org` or `crystallography.net` - If protein or mineral found

**Sometimes made:**
3. `cdn.rcsb.org` - RCSB structure images
4. `eutils.ncbi.nlm.nih.gov` - Entrez fallback (only if PubChem fails)

### Advantages:

- ✅ **Fast:** Local datasets provide instant suggestions
- ✅ **Comprehensive:** API calls provide complete data
- ✅ **Intelligent routing:** Automatically chooses best data source
- ✅ **Fallback chain:** If one API fails, tries another
- ✅ **Hybrid:** Works offline for local searches, fetches live data on demand

---

## What the Extension Does With This Data:

When you search for a molecule in the extension:

1. **Extension** sends: `http://localhost:8001/search?q=hemoglobin`
2. **Search server** returns: Full JSON with 3D coords, images, embed URLs
3. **Extension** displays:
   - RCSB protein image
   - "View 3D" button → `https://molview.org/?pdbid=4hhb`
   - 2D/3D renderer based on coordinates

