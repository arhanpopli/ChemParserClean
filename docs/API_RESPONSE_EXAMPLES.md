# API Response Examples - Source Information

## API Returns Source Information in All Responses

### Endpoint: `/api/nomenclature-to-smiles`

**Request:**
```json
POST /api/nomenclature-to-smiles
{
    "nomenclature": "ethanol"
}
```

**Response (from OPSIN):**
```json
{
    "error": null,
    "smiles": "C(C)O",
    "nomenclature": "ethanol",
    "source": "OPSIN Parser"
}
```

**Response (from PubChem):**
```json
{
    "error": null,
    "smiles": "CCO",
    "nomenclature": "ethanol",
    "source": "PubChem (CID: 702)"
}
```

**Response (from ChemDoodle):**
```json
{
    "error": null,
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "nomenclature": "aspirin",
    "source": "ChemDoodle Database"
}
```

**Response (from Fallback):**
```json
{
    "error": null,
    "smiles": "c1ccccc1",
    "nomenclature": "benzene",
    "source": "Fallback Dictionary"
}
```

**Response (Failed lookup):**
```json
{
    "error": "Compound 'fakemolecule999' not found in PubChem database",
    "smiles": null,
    "nomenclature": "fakemolecule999",
    "source": null
}
```

---

## Endpoint: `/api/nomenclature-to-svg`

**Request:**
```json
POST /api/nomenclature-to-svg
{
    "nomenclature": "aspirin",
    "width": 400,
    "height": 400,
    "options": {
        "aromatic_circles": true
    }
}
```

**Response:**
```json
{
    "error": null,
    "svg": "<svg>...</svg>",
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "nomenclature": "aspirin",
    "source": "ChemDoodle Database",
    "info": {
        "molecular_weight": 180.16,
        "formula": "C9H8O4",
        "num_atoms": 13,
        "num_bonds": 12,
        "num_aromatic_rings": 1,
        "logp": 1.19,
        "hbd": 2,
        "hba": 4
    }
}
```

---

## Source Values

The `source` field tells you WHERE the SMILES came from:

| Source Value | Meaning | Speed |
|--------------|---------|-------|
| `ChemDoodle Database` | Found in pre-loaded ChemDoodle compounds | <1ms |
| `OPSIN Parser` | Parsed using OPSIN IUPAC parser | 2-5s |
| `Fallback Dictionary` | Found in Python fallback dictionary | <1ms |
| `PubChem (CID: xxx)` | Found in PubChem API (with Compound ID) | 500-1000ms |
| `null` | Compound not found in any source | N/A |

---

## Frontend Display

The frontend can show this information to the user:

```html
<div id="info-box">
    <p><strong>Compound:</strong> Aspirin</p>
    <p><strong>SMILES:</strong> CC(=O)Oc1ccccc1C(=O)O</p>
    <p><strong>Source:</strong> ChemDoodle Database</p>
    <p><strong>Molecular Weight:</strong> 180.16 g/mol</p>
</div>
```

---

## JavaScript Example

```javascript
async function fetchCompound(nomenclature) {
    const response = await fetch('/api/nomenclature-to-svg', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            nomenclature: nomenclature,
            options: {...}
        })
    });
    
    const data = await response.json();
    
    if (data.error) {
        console.error('Error:', data.error);
    } else {
        console.log('SMILES found from:', data.source);  // ← Shows source
        console.log('Molecular weight:', data.info.molecular_weight);
        // Render SVG...
    }
}
```

---

## Testing with curl

```bash
# Test OPSIN Parser source
curl -X POST http://localhost:5000/api/nomenclature-to-smiles \
  -H "Content-Type: application/json" \
  -d '{"nomenclature": "ethanol"}' | jq .

# Test ChemDoodle Database source
curl -X POST http://localhost:5000/api/nomenclature-to-smiles \
  -H "Content-Type: application/json" \
  -d '{"nomenclature": "aspirin"}' | jq .

# Test PubChem source (compound not in other sources)
curl -X POST http://localhost:5000/api/nomenclature-to-smiles \
  -H "Content-Type: application/json" \
  -d '{"nomenclature": "acetone"}' | jq .

# Test error handling
curl -X POST http://localhost:5000/api/nomenclature-to-smiles \
  -H "Content-Type: application/json" \
  -d '{"nomenclature": "fakemolecule999"}' | jq .
```

---

## Summary

✅ **Source information is already included** in all API responses
✅ **Tells you WHERE** the SMILES came from (OPSIN, PubChem, ChemDoodle, etc.)
✅ **Easy to display** in frontend UI
✅ **Helps users understand** reliability of the result
✅ **Useful for debugging** if lookups fail
