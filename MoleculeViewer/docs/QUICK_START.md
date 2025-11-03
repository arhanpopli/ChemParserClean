# MoleculeViewer - Quick Start Guide

## Installation

### 1. Install Dependencies

```bash
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer
pip install -r requirements.txt
```

### 2. Setup OPSIN (Optional but Recommended)

For IUPAC nomenclature support:

```bash
python setup_opsin.py
```

This downloads and extracts `opsin-cli.jar` for converting IUPAC names to SMILES.

---

## Running the Server

### Option 1: Direct Command
```bash
python start.py
```

### Option 2: Using Batch File (Windows)
```bash
moleculeviewer.bat
```

### Option 3: Using Flask Directly
```bash
flask run --host 0.0.0.0 --port 5000
```

The server will start on: **http://localhost:5000**

---

## Basic Usage

### 1. Web Interface
Open http://localhost:5000 in your browser for the interactive UI.

### 2. API Example: Get Benzene Structure
```bash
curl -X POST http://localhost:5000/api/nomenclature-to-svg \
  -H "Content-Type: application/json" \
  -d '{"nomenclature":"benzene"}'
```

### 3. Python Example
```python
import requests

response = requests.post('http://localhost:5000/api/nomenclature-to-svg', json={
    'nomenclature': 'aspirin',
    'width': 600,
    'height': 500
})

data = response.json()
if data['error'] is None:
    print(f"SMILES: {data['smiles']}")
    print(f"Molecular Weight: {data['info']['molecular_weight']} g/mol")
    # Save SVG
    with open('aspirin.svg', 'w') as f:
        f.write(data['svg'])
else:
    print(f"Error: {data['error']}")
```

---

## Supported Features

### Nomenclature Lookup
Supports ~17+ ChemDoodle pre-loaded compounds + IUPAC nomenclature via OPSIN.

### Compound Library
Built-in support for:
- Common drugs (aspirin, ibuprofen, acetaminophen, caffeine)
- Aromatic compounds (benzene, naphthalene, pyridine)
- Complex compounds (morphine, strychnine, cubane)

### Molecular Visualization
- SVG output (scalable, transparent background)
- Aromatic bond representation
- 2D coordinate generation

### Molecular Properties
- Molecular weight
- Chemical formula
- Atom/bond count
- LogP (lipophilicity)
- H-bond donors/acceptors
- Aromatic ring count

---

## Common Compounds

```
aspirin         → CC(=O)Oc1ccccc1C(=O)O
benzene         → c1ccccc1
caffeine        → CN1C=NC2=C1C(=O)N(C(=O)N2C)C
naphthalene     → c1cc2ccccc2cc1
phenol          → Oc1ccccc1
pyridine        → c1ccncc1
```

For complete list, see API_REFERENCE.md

---

## Troubleshooting

### Port Already in Use
If port 5000 is in use:
```bash
python start.py --port 5001
```

### OPSIN Not Working
- Verify `opsin-cli.jar` exists in the MoleculeViewer root
- Check Java is installed: `java -version`
- Run setup: `python setup_opsin.py`

### Import Errors
Reinstall dependencies:
```bash
pip install --upgrade rdkit flask flask-cors
```

### SVG Not Rendering
- Ensure SMILES is valid (try known compounds first)
- Check browser console for rendering errors
- Verify RDKit is properly installed

---

## Architecture

```
MoleculeViewer/
├── app/
│   ├── __init__.py           (Flask app initialization)
│   ├── api.py                (REST endpoints)
│   ├── chemistry.py          (SMILES/nomenclature logic)
│   └── chemdoodle_compounds.py (17 pre-extracted compounds)
├── templates/
│   └── index.html            (Web UI)
├── static/                   (CSS, JS assets)
├── start.py                  (Entry point)
├── test_complexes.py         (Integration tests)
├── requirements.txt          (Python dependencies)
└── opsin-cli.jar             (OPSIN binary, ~13.8 MB)
```

---

## Multi-Tier Lookup System

When you request a compound name:

```
User: "aspirin"
  ↓
[Tier 1] ChemDoodle Database? → FOUND ✓
  ↓
Return: "CC(=O)Oc1ccccc1C(=O)O"
```

If ChemDoodle doesn't have it:

```
User: "2-methylnaphthalene"
  ↓
[Tier 1] ChemDoodle Database? → NOT FOUND ✗
  ↓
[Tier 2] OPSIN Parser? → TRY PARSING ✓
  ↓
Return: SMILES from OPSIN
```

If neither work:

```
User: "random_compound"
  ↓
[Tier 1] ChemDoodle? → NOT FOUND ✗
  ↓
[Tier 2] OPSIN? → FAILED ✗
  ↓
[Tier 3] Fallback Dictionary? → NOT FOUND ✗
  ↓
[Tier 4] PubChem API? → TRY ONLINE ✓
  ↓
Return: SMILES from PubChem
```

---

## Performance

- **ChemDoodle Lookup**: ~0ms (instant)
- **OPSIN Parsing**: ~500-2000ms (Java startup + parsing)
- **SVG Generation**: ~100-500ms (depends on molecule complexity)
- **Total for compound**: Typically 1-3 seconds

---

## Next Steps

- See **API_REFERENCE.md** for complete endpoint documentation
- See **ARCHITECTURE.md** for technical details
- Check **test_complexes.py** for more examples
