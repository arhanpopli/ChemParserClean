# 3Dmol.js Implementation Plan: Remove MolView Dependency

## Current State

The `chem-extension/3dmol-viewer.js` currently:
1. Tries to fetch 3D SDF from PubChem
2. Falls back to MolView iframe when PubChem has no 3D data

## Goal

Implement the CIR (Chemical Identifier Resolver) fallback chain directly in 3dmol-viewer.js to handle molecules like cardiolipin without needing MolView.

## Implementation Steps

### Phase 1: CIR Integration for Compounds

Add CIR fallback to `loadMolecule()` function:

```javascript
// NEW: Try CIR when PubChem fails
async function fetch3DFromCIR(smiles) {
    const url = `https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(smiles)}/file?format=sdf&get3d=True`;
    const response = await fetch(url);
    if (response.ok) {
        const text = await response.text();
        // CIR returns HTML 404 page if not found
        if (!text.includes('<h1>Page not found')) {
            return text;
        }
    }
    return null;
}

async function getSMILESFromPubChem(cid) {
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/IsomericSMILES/JSON`;
    const response = await fetch(url);
    if (response.ok) {
        const data = await response.json();
        return data.PropertyTable?.Properties?.[0]?.IsomericSMILES;
    }
    return null;
}
```

### Phase 2: Update loadMolecule() Fallback Chain

```javascript
async function loadMolecule(cid) {
    let sdfData = null;

    // Step 1: Try PubChem 3D
    try {
        sdfData = await fetchPubChem3D(cid);
    } catch (e) { }

    // Step 2: CIR fallback using SMILES
    if (!sdfData) {
        console.log('PubChem 3D failed, trying CIR...');
        const smiles = await getSMILESFromPubChem(cid);
        if (smiles) {
            sdfData = await fetch3DFromCIR(smiles);
        }
    }

    // Step 3: 2D as 3D fallback (flat molecule)
    if (!sdfData) {
        console.log('CIR failed, using 2D structure...');
        sdfData = await fetchPubChem2D(cid);
    }

    // Step 4: Final fallback to MolView (only if all else fails)
    if (!sdfData) {
        fallbackToMolView(cid, true);
        return;
    }

    // Render with 3Dmol.js
    renderMolecule(sdfData);
}
```

### Phase 3: Biomolecule Support (RCSB)

The search server already returns PDB URLs. Update to handle PDB directly:

```javascript
async function loadBiomolecule(pdbid) {
    const url = `https://files.rcsb.org/view/${pdbid}.pdb`;
    const response = await fetch(url);
    if (response.ok) {
        const pdbData = await response.text();
        renderPDB(pdbData);
    }
}

function renderPDB(pdbData) {
    const viewer = $3Dmol.createViewer('viewer', { backgroundColor: bgColor });
    viewer.addModel(pdbData, 'pdb');
    viewer.setStyle({}, { cartoon: { color: 'spectrum' }});
    viewer.zoomTo();
    viewer.render();
}
```

### Phase 4: Mineral Support (COD)

Use the search server's COD integration:

```javascript
async function loadMineral(codid) {
    // Try direct COD fetch
    const url = `https://www.crystallography.net/cod/${codid}.cif`;

    // Need CORS proxy or server-side fetch
    // Option 1: Use search server to proxy
    // Option 2: Use CORS proxy service

    const proxyUrl = `http://localhost:8001/proxy?url=${encodeURIComponent(url)}`;
    const response = await fetch(proxyUrl);
    if (response.ok) {
        const cifData = await response.text();
        renderCIF(cifData);
    }
}

function renderCIF(cifData) {
    const viewer = $3Dmol.createViewer('viewer', { backgroundColor: bgColor });
    viewer.addModel(cifData, 'cif');
    viewer.setStyle({}, { sphere: { colorscheme: 'Jmol' }});
    viewer.zoomTo();
    viewer.render();
}
```

## API Endpoints Reference

| Purpose | URL |
|---------|-----|
| PubChem 3D SDF | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d` |
| PubChem 2D SDF | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=2d` |
| PubChem SMILES | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON` |
| CIR 3D from SMILES | `https://cactus.nci.nih.gov/chemical/structure/{smiles}/file?format=sdf&get3d=True` |
| RCSB PDB | `https://files.rcsb.org/view/{pdbid}.pdb` |
| COD CIF | `https://www.crystallography.net/cod/{codid}.cif` |
| Local Search | `http://localhost:8001/search?q={query}` |

## Server-Side Enhancements (Optional)

Add to search server (port 8001) for better reliability:

### 1. CIR Proxy Endpoint
```javascript
app.get('/cir/3d', async (req, res) => {
    const smiles = req.query.smiles;
    const url = `https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(smiles)}/file?format=sdf&get3d=True`;
    const response = await fetch(url);
    res.send(await response.text());
});
```

### 2. COD Proxy Endpoint
```javascript
app.get('/cod/cif', async (req, res) => {
    const codid = req.query.codid;
    const url = `https://www.crystallography.net/cod/${codid}.cif`;
    const response = await fetch(url);
    res.send(await response.text());
});
```

### 3. Enhanced Search Response
```json
{
    "query": "cardiolipin",
    "found": true,
    "category": "compound",
    "cid": 166177218,
    "smiles": "CCCCCC...",
    "sdf_3d_available": false,
    "fallback_method": "cir",
    "embed_url": "3dmol-viewer.html?cid=166177218"
}
```

## Testing Checklist

- [ ] Aspirin (CID: 2244) - Has PubChem 3D
- [ ] Cardiolipin (CID: 166177218) - No PubChem 3D, needs CIR fallback
- [ ] Insulin (PDB: 4INS) - Biomolecule/protein
- [ ] Calcite (CODID: 1010928) - Mineral
- [ ] Rhinovirus (PDB: 4RHV) - Large biomolecule

## File Changes Required

1. **chem-extension/3dmol-viewer.js**
   - Add CIR fetch function
   - Add SMILES fetch function
   - Update loadMolecule() fallback chain
   - Add PDB rendering support
   - Add CIF rendering support

2. **Molview/molview/search-server.js** (optional)
   - Add CIR proxy endpoint
   - Add COD proxy endpoint
   - Enhance search response with 3D availability info

## Priority Order

1. **High**: CIR fallback for compounds (fixes cardiolipin issue)
2. **Medium**: Direct PDB loading for biomolecules
3. **Low**: Direct CIF loading for minerals (CORS issues)

## Notes

- CIR service may be slow for complex molecules
- Consider caching CIR results in search server
- 2D as 3D is acceptable as final fallback (flat molecule better than no molecule)
- Keep MolView fallback as emergency option for edge cases
