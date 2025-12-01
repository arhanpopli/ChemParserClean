# ChemParser Unified Search System

## Overview

The ChemParser extension now uses a **unified search API** running on port 8001 that handles ALL chemical queries with autocorrect, name normalization, and intelligent data retrieval from multiple sources (PubChem, PDB, COD minerals database).

## Architecture

```
User types: chem:rhinovirus:
    ↓
Extension sends query to: http://localhost:8001/search?q=rhinovirus
    ↓
Search API:
  1. Searches local databases (PDB proteins, COD minerals)
  2. Queries PubChem autocomplete
  3. Autocorrects typos using similarity matching
  4. Returns best match with:
     - Corrected name
     - Canonical SMILES
     - Isomeric SMILES (for stereochemistry)
     - Compound type (compound/biomolecule/mineral)
     - Source URLs
    ↓
Extension receives:
  {
    "query": "rhinovirus",
    "corrected_query": null,  // or corrected name if typo
    "name": "Rhinovirus",
    "canonical_smiles": "...",
    "isomeric_smiles": "...",
    "primary_type": "biomolecule",
    "pdbid": "1RHV",
    "source_url": "https://www.rcsb.org/structure/1RHV"
  }
    ↓
Extension renders molecule using:
  - MoleculeViewer (default)
  - OR selected rendering engine
```

## Starting the System

### Quick Start (All Servers)
```bash
1-start-all.bat
```

This starts:
- **Port 5001**: Mol2ChemFig Server (Flask)
- **Port 5000**: MoleculeViewer Backend (Node.js)
- **Port 5002**: PubChem Server (Node.js)
- **Port 8000**: MolView PHP Server (main viewer)
- **Port 8001**: **MolView Search API** (unified search)

### Individual Servers
```bash
# MolView Search API (REQUIRED)
cd Molview\molview
node search-server.js

# MolView PHP Server (needed for embeds)
cd Molview\molview
php -S localhost:8000

# MoleculeViewer (for rendering)
cd MoleculeViewer
node server.js
```

## Using the Extension

### 1. Select "MolView Search" Engine

In the extension popup:
1. Open the ChemParser extension popup
2. Select **"MolView Search"** rendering engine
3. Reload the page

### 2. Type Chemical Queries

**Correct names:**
- `chem:ethanol:`
- `chem:aspirin:`
- `chem:rhinovirus:`
- `chem:quartz:`

**Typos (will be autocorrected):**
- `chem:booxite:` → Autocorrects to **Brookite**
- `chem:aspriin:` → Autocorrects to **Aspirin**
- `chem:etanol:` → Autocorrects to **Ethanol**

**SMILES (also supported):**
- `chem:CCO:`
- `chem:CC(=O)OC1=CC=CC=C1C(=O)O:`

### 3. Autocorrect Display

When a typo is detected, the extension shows:

```
✓ Autocorrected: booxite → Brookite
[Molecular structure displayed below]
```

## Features

### ✅ Autocorrect System
- Fuzzy string matching using `similar_text` algorithm
- Prefix matching bonus (words starting with query get higher rank)
- Works for compounds, proteins, and minerals

### ✅ Multi-Database Search
1. **PubChem**: Small molecules, drugs, compounds
2. **PDB (Protein Data Bank)**: Proteins, biomolecules (e.g., rhinovirus)
3. **COD (Crystallography Open Database)**: Minerals, crystals (e.g., brookite, quartz)

### ✅ Stereochemistry Support
- Returns **canonical SMILES** (default)
- Returns **isomeric SMILES** (if "Use 3D SMILES" is enabled)
- Extension automatically uses correct SMILES based on settings

### ✅ Unified Data
- No need to check multiple APIs
- Single endpoint for all queries
- Consistent response format
- Automatic fallback between sources

## API Response Format

### Compact Mode (Default)
```json
{
  "query": "booxite",
  "corrected_query": "Brookite",
  "name": "Brookite",
  "canonical_smiles": null,
  "isomeric_smiles": null,
  "sdf": {
    "available": true,
    "size_bytes": 230,
    "lines": 8,
    "format": "CIF"
  },
  "source_url": "http://www.crystallography.net/cod/9004137.html",
  "primary_type": "mineral",
  "codid": "9004137",
  "embed_url": "http://localhost:8000/?codid=9004137"
}
```

### Full Mode (with complete structure data)
```bash
http://localhost:8001/search?q=aspirin&format=full
```

Returns full SDF/PDB/CIF content in the `sdf` field.

## Benefits

### For Users
- ✅ **Typo-tolerant**: No need to type perfect chemical names
- ✅ **Smart**: Automatically finds the right compound/protein/mineral
- ✅ **Fast**: Single API call gets all data
- ✅ **Comprehensive**: Works for compounds, proteins, and minerals

### For Developers
- ✅ **Simple integration**: One API endpoint for everything
- ✅ **Consistent format**: All responses follow same structure
- ✅ **Extensible**: Easy to add more databases
- ✅ **Cached**: Results cached by search API for performance

## Examples

### Example 1: Compound with Typo
```
Input: chem:aspriin:
API Query: http://localhost:8001/search?q=aspriin
Result:
  - Corrected to "Aspirin"
  - SMILES: CC(=O)OC1=CC=CC=C1C(=O)O
  - Rendered as 2D structure
  - Shows: "✓ Autocorrected: aspriin → Aspirin"
```

### Example 2: Protein
```
Input: chem:rhinovirus:
API Query: http://localhost:8001/search?q=rhinovirus
Result:
  - Found in PDB database
  - PDB ID: 1RHV (or similar)
  - Rendered as 3D protein structure
  - No autocorrect needed
```

### Example 3: Mineral with Typo
```
Input: chem:booxite:
API Query: http://localhost:8001/search?q=booxite
Result:
  - Corrected to "Brookite"
  - COD ID: 9004137
  - Rendered as crystal structure
  - Shows: "✓ Autocorrected: booxite → Brookite"
```

## Troubleshooting

### "localhost refused to connect"
**Problem**: Search API not running on port 8001

**Solution**:
```bash
cd Molview\molview
node search-server.js
```

### "No results found"
**Problem**: Query doesn't match any known compound/protein/mineral

**Solution**:
- Try alternative spellings
- Use SMILES instead: `chem:CCO:`
- Check if compound exists in PubChem

### Servers not starting
**Problem**: Port conflicts or missing dependencies

**Solution**:
```bash
# Kill all processes
util-stop-all.bat

# Restart
1-start-all.bat
```

## Configuration

### Enable 3D SMILES (Stereochemistry)
1. Open extension popup
2. Enable **"Use 3D SMILES for stereochemistry"**
3. Search API will return isomeric SMILES with stereochemistry markers

### Change Rendering Engine
Even with MolView Search, you can choose how molecules are rendered:
- **MoleculeViewer** (default, SVG diagrams)
- **mol2chemfig** (LaTeX-style chemistry)
- **PubChem** (direct PubChem images)
- **Client-Side** (offline rendering with SmilesDrawer)

The search API autocorrect works with ALL rendering engines!

## Technical Details

### Search Algorithm
1. Query local databases (PDB, COD minerals)
2. Query PubChem autocomplete API
3. Merge results and remove duplicates
4. Sort by similarity score:
   - Similar_text algorithm (longest common substring)
   - Prefix bonus (+100 if query matches start of name)
5. Return best match

### Minimum Similarity Threshold
- Configured in `search-server.js`: `MIN_SIM = 40`
- Adjustable based on desired strictness

### Supported Compound Types
- **compound**: Small molecules from PubChem
- **biomolecule**: Proteins, nucleic acids from PDB
- **mineral**: Crystal structures from COD

## API Endpoints

### Search Endpoint
```
GET http://localhost:8001/search?q={query}
GET http://localhost:8001/embed/v1/search?q={query}
```

**Query Parameters:**
- `q` or `text`: Search query (required)
- `format`: `compact` (default) or `full`

**Response:** JSON with compound data

### Health Check
```bash
# Check if API is running
curl http://localhost:8001/search?q=aspirin
```

Should return JSON with aspirin data.

## Future Enhancements

Potential improvements:
- [ ] Add more databases (ChEBI, DrugBank, etc.)
- [ ] Implement caching layer (Redis)
- [ ] Add synonym support
- [ ] Implement query suggestions
- [ ] Add batch search endpoint
- [ ] Support for chemical reactions

## Credits

- **MolView**: Herman Bergwerf (https://molview.org)
- **PubChem**: National Library of Medicine
- **PDB**: RCSB Protein Data Bank
- **COD**: Crystallography Open Database
- **ChemParser**: Unified search implementation

---

**Last Updated**: 2025-01-28
**Version**: 3.0 (Unified Search)
