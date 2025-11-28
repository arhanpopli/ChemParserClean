# MolView API Implementation Guide

## Overview
The browser extension now queries your **local MolView instance** at `http://localhost:5003` to intelligently determine how to render molecules.

## Required API Endpoint

You need to implement this endpoint in your local MolView server:

### `GET /api/search?q={molecule_name}`

**Purpose:** Search for a molecule and return its category and data.

### Response Format

```json
{
  "found": boolean,
  "category": string,  // "compound", "mineral", "protein", "biomolecule", "pdb", "crystal"
  "smiles": string,    // Optional: for compounds/minerals
  "sdf": string,       // Optional: for compounds/minerals (preferred)
  "pdb_id": string     // Optional: for proteins
}
```

### Example Responses

#### 1. Small Molecule (Compound)
```json
{
  "found": true,
  "category": "compound",
  "smiles": "CCO",
  "sdf": "... SDF data ..."
}
```

#### 2. Protein/Biomolecule
```json
{
  "found": true,
  "category": "protein",
  "pdb_id": "1RHV"
}
```

#### 3. Not Found
```json
{
  "found": false
}
```

## How the Extension Uses This

### Flow Diagram

```
User types: chem:rhinovirus:+3d:
    ↓
Extension queries: http://localhost:5003/api/search?q=rhinovirus
    ↓
MolView API responds with category
    ↓
┌─────────────────────────────────────┐
│ Category: "compound" or "mineral"   │
│ → Use SMILES/SDF                    │
│ → Render locally with 3Dmol.js      │
└─────────────────────────────────────┘
┌─────────────────────────────────────┐
│ Category: "protein" or "biomolecule"│
│ → Embed local MolView iframe        │
│ → URL: localhost:5003/?q=rhinovirus │
└─────────────────────────────────────┘
```

### Fallback Chain

1. **Try local MolView API** (`/api/search`)
2. If API fails → **Try PubChem** (for small molecules)
3. If PubChem fails → **Embed local MolView** (`localhost:5003/?q=...`)

## Testing

### Test Cases

1. **Ethanol** (small molecule): `chem:ethanol:+3d:`
2. **Rhinovirus** (protein): `chem:rhinovirus:+3d:`
3. **Cardiolipin** (complex molecule): `chem:cardiolipin:+3d:`

## Current Status

✅ Extension updated to use `localhost:5003`
✅ API endpoint structure defined
⏳ **You need to implement:** `/api/search` endpoint in your MolView server
