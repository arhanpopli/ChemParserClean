# Search Server MolView Integration Update

## Changes Made

The search server (`Molview/molview/search-server.js`) has been updated to use official MolView embed URLs instead of localhost:8000.

### Updated URLs:

**1. PubChem Compounds (Chemicals):**
- **Old:** `http://localhost:8000/embed/v2/?cid={cid}`
- **New:** `https://embed.molview.org/v1/?cid={cid}`
- **Example:** `https://embed.molview.org/v1/?cid=16617`

**2. Crystallography/Minerals (COD):**
- **Old:** `http://localhost:8000/?codid={codid}`
- **New:** `https://embed.molview.org/v1/?codid={codid}`
- **Example:** `https://embed.molview.org/v1/?codid=5000121`
- **Also updated CIF URL from:** `http://localhost:8000/php/cif.php?codid={codid}`
- **To:** `https://www.crystallography.net/cod/{codid}.cif`

**3. Proteins/Biomolecules (RCSB PDB):**
- **Old:** `http://localhost:8000/embed/v2/?pdbid={pdbid}`
- **New:** `https://molview.org/?pdbid={pdbid}`
- **Example:** `https://molview.org/?pdbid=4rhv`

## Example Search Results

### Calcite (Mineral):
```json
{
    "query": "calcite",
    "corrected_query": "Calcite",
    "name": "Calcite",
    "chemical_name": null,
    "formula": null,
    "canonical_smiles": null,
    "isomeric_smiles": null,
    "sdf": {
        "available": false
    },
    "source_url": "http://www.crystallography.net/cod/1010928.html",
    "primary_type": "mineral",
    "sub_type": null,
    "embed_url": "https://embed.molview.org/v1/?codid=1010928",
    "codid": "1010928",
    "cif_url": "https://www.crystallography.net/cod/1010928.cif",
    "_info": {
        "format": "compact",
        "note": "To get full SDF/PDB/CIF content, add &format=full to the URL"
    }
}
```

### Chemical Compound:
```json
{
    "embed_url": "https://embed.molview.org/v1/?cid=16617",
    "cid": "16617"
}
```

### Protein/Biomolecule:
```json
{
    "embed_url": "https://molview.org/?pdbid=4rhv",
    "pdbid": "4rhv"
}
```

## Benefits

1. **No local server dependency**: The extension can now work without running the local MolView PHP server on port 8000
2. **Always up-to-date**: Uses the official MolView service with latest features and bug fixes
3. **Better reliability**: Official servers have better uptime than local development servers
4. **Simpler deployment**: One less server to manage (only need search-server on port 8001)

## Usage

Start the search server:
```powershell
cd C:\Users\Kapil\Personal\STUFF\Chemparser
node Molview\molview\search-server.js
```

Test queries:
```powershell
# Chemical compound
Invoke-RestMethod -Uri "http://localhost:8001/search?q=caffeine&format=compact"

# Mineral
Invoke-RestMethod -Uri "http://localhost:8001/search?q=calcite&format=compact"

# Protein
Invoke-RestMethod -Uri "http://localhost:8001/search?q=hemoglobin&format=compact"
```

## Files Modified

- `Molview/molview/search-server.js` (lines 410-440, 628)
  - Removed `localMolViewUrl` variable
  - Updated PDB embed_url to use `https://molview.org/`
  - Updated COD embed_url to use `https://embed.molview.org/v1/`
  - Updated COD cif_url to use `https://www.crystallography.net/cod/`
  - Updated PubChem embed_url to use `https://embed.molview.org/v1/`
