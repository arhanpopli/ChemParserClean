MolView's embed/v1 shows **FLAT 2D molecules** for complex molecules (e.g., cardiolipin) while V2 shows proper **3D bent structures**. This is because V1 lacks the Sketcher component needed to extract SMILES.

## ROOT CAUSE FOUND

**V1 shows flat molecules because:**
1. PubChem doesn't have 3D coordinates for cardiolipin
2. V1 tries to use `MolFile` class to extract SMILES but **MolFile is NOT in the v1 build bundle**
3. The try/catch silently fails and loads the flat 2D SDF

**V2 shows bent 3D molecules because:**
1. V2 loads 2D SDF into the **Sketcher** (MolPad 2D editor)
2. Uses `Sketcher.getSMILES()` to extract SMILES from the sketcher
3. Sends SMILES to CIR with `get3d=True`
4. **CIR generates 3D coordinates using molecular mechanics/force field**
5. The 3D structure has proper bond angles, torsions, and conformations

**The magic is `get3d=True` in the CIR API call** - this tells CIR to compute 3D geometry!

## Main App Fallback Chain
```javascript
// Full app uses Sketcher.getSMILES()
Request.PubChem.sdf(cid, false, // Get 3D
  success, // Got 3D from PubChem
  function() { // 3D failed
    var smiles = Sketcher.getSMILES(); // Extract from 2D editor
    Request.CIR.resolve(smiles, false, // get3d=True
      function(mol3d) { Model.loadMOL(mol3d); }, // CIR generated 3D!
      function() { Model.loadMOL(mol2d); } // Final fallback: flat 2D
    );
  });
```

**Embed V1:** Has CIR code but no Sketcher - can't extract SMILES from 2D
**Embed V2:** Full app (50KB), includes Sketcher, uses it to get SMILES
**Why V2 works:** Has Sketcher.getSMILES() which V1 doesn't have

**Key Services:**
- CIR: `https://cactus.nci.nih.gov/chemical/structure/{SMILES}/file?format=sdf&get3d=True`
- PubChem 2D: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sdf?record_type=2d`
- PubChem 3D: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sdf?record_type=3d`

**Data Source Comparison:**
- PubChem (small molecules): 3D often available, use CIR fallback
- PubChem (complex): 3D usually missing, requires CIR from SMILES
- RCSB (proteins): Always has 3D (PDB format)
- COD (minerals/crystals): Always has 3D (CIF format)

**3Dmol.js Implementation:**
1. Try PubChem 3D endpoint
2. Get SMILES from PubChem API (or extract from 2D SDF)
3. Call CIR to generate 3D
4. Fallback to 2D as 3D if all fail

**Required Endpoints:**
```
PubChem 3D: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sdf?record_type=3d
PubChem 2D: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sdf?record_type=2d
PubChem SMILES: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/json
CIR 3D: https://cactus.nci.nih.gov/chemical/structure/{SMILES}/file?format=sdf&get3d=True
RCSB PDB: https://files.rcsb.org/view/{pdbid}.pdb
COD CIF: https://www.crystallography.net/cod/{codid}.cif
```

**Quickstart Implementation:**
```javascript
async function load3DStructure(cid) {
    try {
        return await fetchPubChem3D(cid);
    } catch (e) {
        const smiles = await fetchPubChemSMILES(cid);
        try {
            return await fetchCIR3D(smiles);
        } catch (e) {
            return await fetchPubChem2D(cid); // 2D fallback
        }
    }
}
```

---

## IMPLEMENTATION COMPLETE (2024-12-05)

Added to `chem-extension/3dmol-viewer.js`:

### New Functions
- `getSMILESFromPubChem(cid)` - Fetches IsomericSMILES from PubChem
- `generate3DFromCIR(smiles)` - Generates 3D SDF via CIR with `get3d=True`
- `get2DFromPubChem(cid)` - Gets 2D SDF as last resort fallback

### Updated `loadMolecule()` Fallback Chain
```
Step 1: PubChem 3D conformer → if fails:
Step 2: PubChem computed 3D → if fails:
Step 3: CIR 3D from SMILES (get3d=True) → if fails:
Step 4: 2D as 3D (flat molecule) → if fails:
Step 5: MolView iframe (emergency)
```

### Test with Cardiolipin
```
chem-extension/3dmol-viewer.html?cid=166177218
```
Should now show proper 3D bent structure instead of flat molecule!
