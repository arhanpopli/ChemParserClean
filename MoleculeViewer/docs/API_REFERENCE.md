# MoleculeViewer API Reference

## Overview
MoleculeViewer is a Flask-based web service for converting chemical nomenclature and SMILES strings to SVG visualizations. It uses a multi-tier nomenclature lookup system: ChemDoodle → OPSIN → Fallback Dictionary → PubChem.

## Base URL
```
http://localhost:5000
```

## Endpoints

### 1. Health Check
**GET** `/health`

Check if the service is running.

**Response:**
```json
{
  "status": "ok"
}
```

---

### 2. Nomenclature to SMILES
**POST** `/api/nomenclature-to-smiles`

Convert a chemical name to SMILES string.

**Request:**
```json
{
  "nomenclature": "aspirin"
}
```

**Response (Success):**
```json
{
  "error": null,
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "nomenclature": "aspirin"
}
```

**Response (Error):**
```json
{
  "error": "Could not convert 'aspirin' to SMILES",
  "smiles": null,
  "nomenclature": "aspirin"
}
```

---

### 3. SMILES to SVG
**POST** `/api/smiles-to-svg`

Convert SMILES string to SVG molecule visualization.

**Request:**
```json
{
  "smiles": "c1ccccc1",
  "width": 600,
  "height": 500
}
```

**Response (Success):**
```json
{
  "error": null,
  "svg": "<svg>...</svg>",
  "smiles": "c1ccccc1",
  "info": {
    "molecular_weight": 78.11,
    "formula": "C6H6",
    "num_atoms": 6,
    "num_bonds": 6,
    "num_aromatic_rings": 1,
    "logp": 2.13,
    "hbd": 0,
    "hba": 0
  }
}
```

---

### 4. Nomenclature to SVG
**POST** `/api/nomenclature-to-svg`

Convert chemical name directly to SVG (combines nomenclature and SMILES lookup).

**Request:**
```json
{
  "nomenclature": "caffeine",
  "width": 600,
  "height": 500
}
```

**Response (Success):**
```json
{
  "error": null,
  "svg": "<svg>...</svg>",
  "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  "nomenclature": "caffeine",
  "info": {
    "molecular_weight": 194.19,
    "formula": "C8H10N4O2",
    "num_atoms": 24,
    "num_bonds": 24,
    "num_aromatic_rings": 2,
    "logp": 0.16,
    "hbd": 0,
    "hba": 4
  }
}
```

---

### 5. Molecule Information
**POST** `/api/molecule-info`

Get molecular properties for a given SMILES.

**Request:**
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O"
}
```

**Response:**
```json
{
  "error": null,
  "info": {
    "molecular_weight": 180.16,
    "formula": "C9H8O4",
    "num_atoms": 21,
    "num_bonds": 20,
    "num_aromatic_rings": 1,
    "logp": 1.19,
    "hbd": 2,
    "hba": 4
  }
}
```

---

## Supported Compounds (ChemDoodle Tier-1)

These compounds are available for instant lookup:

- **aspirin** → CC(=O)Oc1ccccc1C(=O)O
- **benzene** → c1ccccc1
- **caffeine** → CN1C=NC2=C1C(=O)N(C(=O)N2C)C
- **cubane** → C12C3C4C1C5C1C2C2C3C3C4C4C5C5C1C1C2C2C3C3C4C4C5C1C12C34
- **cyclobutadiene** → C1=CC=C1
- **cyclohexane** → C1CCCCC1
- **cyclopentadiene** → C1=CCC=C1
- **furan** → O1C=CC=C1
- **hexane** → CCCCCC
- **mannitol** → C(C(C(C(C(CO)O)O)O)O)O
- **morphine** → CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)OC3C(C4)O
- **napthalene** → c1cc2ccccc2cc1
- **phenol** → Oc1ccccc1
- **pyridine** → c1ccncc1
- **pyrrole** → c1cc[nH]c1
- **strychnine** → CN1C[C@H]2[C@H]3[C@H]4C=C[C@H]([C@@]45CCN3[C@@H]2C1)Oc1c5ccc2ccccc12
- **thiophene** → c1ccsc1

---

## Error Handling

All endpoints return appropriate HTTP status codes:

- **200** - Success
- **400** - Bad request (missing required fields, invalid SMILES)
- **404** - Compound not found
- **500** - Server error

Error responses include an `error` field with descriptive message:
```json
{
  "error": "Invalid SMILES string",
  "smiles": null
}
```

---

## Multi-Tier Nomenclature System

When converting nomenclature to SMILES, the system tries in order:

1. **Tier 1 - ChemDoodle Database**: Pre-extracted compounds for instant lookup (~0ms)
2. **Tier 2 - OPSIN Parser**: Java-based IUPAC nomenclature parser (requires opsin-cli.jar)
3. **Tier 3 - Fallback Dictionary**: 20+ common compounds and drugs
4. **Tier 4 - PubChem API**: Online database lookup (requires internet)

---

## Examples

### Example 1: Get benzene structure
```bash
curl -X POST http://localhost:5000/api/nomenclature-to-svg \
  -H "Content-Type: application/json" \
  -d '{"nomenclature":"benzene"}'
```

### Example 2: Analyze aspirin
```bash
curl -X POST http://localhost:5000/api/molecule-info \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}'
```

### Example 3: Convert SMILES to image
```bash
curl -X POST http://localhost:5000/api/smiles-to-svg \
  -H "Content-Type: application/json" \
  -d '{"smiles":"c1ccccc1","width":600,"height":500}'
```
