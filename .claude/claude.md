# ChemParser Project Context

## Project Overview
ChemParser is a chemical structure visualization system with Chrome extension integration. It consists of three main servers:

1. **MoleculeViewer Server** (Port 5000) - Node.js server for SMILES/nomenclature to SVG conversion
2. **Mol2ChemFig Server** (Port 5001) - Flask server wrapping mol2chemfig Docker backend for superior rendering
3. **PubChem Server** (Port 5002) - Flask server for fetching images and 3D models from PubChem

## Current File Structure

### Servers
- `MoleculeViewer/server.js` - Node.js server on port 5000
- `mol2chemfig_server.py` - Flask server on port 5001
- `pubchem_server.py` - Flask server on port 5002

### Cache Directories
- `MoleculeViewer/cache/moleculeviewer/` - MoleculeViewer cache
- `cache/mol2chemfig/` - Mol2ChemFig cache
- `pubchem-cache/` - PubChem cache

### Extension
- `chem-extension/` - Chrome extension files
  - `content.js` - Content script for text replacement
  - `popup.js` - Extension popup
  - `popup.html` - Extension popup UI

## Key Technical Details

### Cache System
Both servers use hash-based caching:
- **MoleculeViewer**: MD5 hash of `type:value:widthxheight`
- **Mol2ChemFig**: SHA256 hash of `smiles:sorted_options`
- **PubChem**: MD5 hash of `img_cid_size_type`

### Docker Backend
Mol2ChemFig uses Docker container at `http://localhost:8000`:
- `/m2cf/submit` - Simple SMILES to ChemFig
- `/m2cf/apply` - SMILES with options
- `/m2cf/layers` - Layered SVG generation
- `/m2cf/search` - Name to SMILES lookup

### Extension Integration
Extension replaces `chem:chemicalname:` patterns with images:
- Tries MoleculeViewer first
- Falls back to Mol2ChemFig if needed
- Can use PubChem for direct image links

## Active Issues to Fix

1. **Cache Deduplication** - Same SMILES stored with different cache keys
2. **Options Persistence** - ChemFig options don't persist across searches
3. **Cache Link Display** - Links not shown after applying options
4. **MoleculeViewer Generation** - Not generating new SVGs for some compounds
5. **3D SMILES Support** - Need to integrate OPSIN for 3D stereochemistry

## Testing Approach
- Each server should be testable independently
- Test files available: `test_m2cf_full.html`, `test_pubchem.html`
- Extension can be tested by loading unpacked extension

## Important Notes
- All servers already have separate cache folders (DONE)
- PubChem server already exists (DONE)
- Focus on fixing bugs and optimizations
- Maintain backward compatibility with existing cache
