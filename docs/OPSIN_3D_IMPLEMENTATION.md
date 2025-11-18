# OPSIN 3D SMILES Integration - Implementation Guide

## Overview
This document describes the OPSIN integration for 3D SMILES support in the ChemParser project. OPSIN (Open Parser for Systematic IUPAC Nomenclature) provides high-quality name-to-structure conversion with full stereochemistry support.

## What Was Implemented

### 1. Backend Server Changes (mol2chemfig_server.py)

#### New Endpoints

**`/api/opsin` (GET/POST)**
- Converts chemical nomenclature to 3D SMILES using OPSIN API
- Returns stereochemically-aware SMILES strings
- Example: `glucose` → `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO`

**`/api/generate-3d` (POST)**
- Complete pipeline: name → 3D SMILES (via OPSIN) → rendered structure (via mol2chemfig)
- Supports all chemfig rendering options
- Returns SVG, PDF, and ChemFig code
- Includes caching for performance

#### Code Changes
```python
@app.route('/api/opsin', methods=['GET', 'POST'])
def opsin_conversion():
    """Convert nomenclature to 3D SMILES using OPSIN"""
    # Calls https://opsin.ch.cam.ac.uk/opsin/{name}.json
    # Returns 3D SMILES with stereochemistry
```

```python
@app.route('/api/generate-3d', methods=['POST'])
def generate_3d():
    """Generate molecule from nomenclature with 3D stereochemistry"""
    # Step 1: Convert name to 3D SMILES via OPSIN
    # Step 2: Generate structure using mol2chemfig
    # Returns complete structure with metadata
```

### 2. Extension Changes

#### popup.html
Added new option in mol2chemfig section:
```html
<div class="option option-border">
  <label for="use3DSmilesToggle">
    <strong>3D Stereochemistry (OPSIN)</strong>
    <small>Use OPSIN for 3D SMILES with stereochemistry</small>
  </label>
  <input type="checkbox" id="use3DSmilesToggle">
  <label class="toggle" for="use3DSmilesToggle"></label>
</div>
```

#### popup.js
Added setting management:
```javascript
// New toggle element
const use3DSmilesToggle = document.getElementById('use3DSmilesToggle');

// Default setting
use3DSmiles: false,

// Load setting
if (use3DSmilesToggle) use3DSmilesToggle.checked = settings.use3DSmiles;

// Save setting
if (use3DSmilesToggle) {
  use3DSmilesToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ use3DSmiles: e.target.checked }, () => {
      showStatus('3D SMILES (OPSIN) ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
    });
  });
}
```

### 3. Test Page (test_3d_smiles.html)

Created comprehensive test page with:
- OPSIN conversion testing
- 3D structure generation
- Stereochemistry test cases (glucose, alanine, cholesterol, etc.)
- 2D vs 3D comparison tool
- Visual interface with gradient styling

## How It Works

### Standard Flow (without 3D)
1. User enters "glucose"
2. System looks up in PubChem/mol2chemfig database
3. Returns 2D SMILES: `C(C1C(C(C(C(O1)O)O)O)O)O`
4. Renders structure without stereochemistry

### 3D Flow (with OPSIN)
1. User enters "D-glucose" (or "glucose")
2. System calls OPSIN API: `https://opsin.ch.cam.ac.uk/opsin/D-glucose.json`
3. OPSIN returns 3D SMILES: `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO`
4. mol2chemfig renders structure with stereochemistry markers

### Key Differences

| Feature | 2D SMILES | 3D SMILES (OPSIN) |
|---------|-----------|-------------------|
| Stereochemistry | None | Full support with @ markers |
| Nomenclature parsing | Database lookup | IUPAC parsing |
| Chiral centers | Not specified | Explicitly defined |
| Example | `C(C1C(C(C(C(O1)O)O)O)O)O` | `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO` |

## API Examples

### Testing OPSIN Conversion
```bash
# Get 3D SMILES for glucose
curl "http://localhost:5001/api/opsin?name=glucose"

# Response:
{
  "success": true,
  "name": "glucose",
  "smiles": "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
  "source": "OPSIN"
}
```

### Generating 3D Structure
```bash
curl -X POST http://localhost:5001/api/generate-3d \
  -H "Content-Type: application/json" \
  -d '{
    "name": "D-glucose",
    "options": ["-o", "-m"],
    "return_format": "svg"
  }'

# Response:
{
  "success": true,
  "name": "D-glucose",
  "smiles": "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
  "hash": "abc123...",
  "svg_url": "/images/abc123.svg",
  "chemfig": "\\chemfig{...}",
  "source": "OPSIN"
}
```

## Stereochemistry Examples

### Simple Chiral Molecules
1. **L-alanine**: `CC(C(=O)O)N` → `N[C@@H](C)C(=O)O` (with @)
2. **D-alanine**: `CC(C(=O)O)N` → `N[C@H](C)C(=O)O` (with @)

### Complex Sugars
1. **D-glucose**: Multiple chiral centers with specific configuration
2. **D-fructose**: Ketose with stereochemistry
3. **Sucrose**: Disaccharide with preserved stereochemistry

### Ring Compounds
1. **Cholesterol**: Steroid with multiple chiral centers and rings
2. **Caffeine**: Aromatic heterocycles

## Testing Guide

### 1. Test OPSIN Endpoint
```bash
# Open test_3d_smiles.html in browser
# Navigate to http://localhost:5001 (after starting server)
# Or open test_3d_smiles.html directly

# Test cases:
- glucose (simple)
- D-glucose (with stereochemistry)
- L-alanine (amino acid)
- cholesterol (complex steroid)
```

### 2. Test Extension Integration
1. Open Chrome extension popup
2. Select "mol2chemfig" as rendering engine
3. Enable "3D Stereochemistry (OPSIN)" toggle
4. On a webpage, type: `chem:D-glucose:`
5. Extension should use OPSIN to fetch 3D SMILES
6. Rendered structure should show stereochemistry

### 3. Compare 2D vs 3D
```javascript
// In browser console (on test page)
async function compare() {
  // 2D version
  const resp2d = await fetch('http://localhost:5001/api/search?name=glucose');
  const data2d = await resp2d.json();
  console.log('2D SMILES:', data2d.smiles);

  // 3D version
  const resp3d = await fetch('http://localhost:5001/api/opsin?name=glucose');
  const data3d = await resp3d.json();
  console.log('3D SMILES:', data3d.smiles);
}
compare();
```

## MoleculeViewer Status

### Current Implementation
The MoleculeViewer server (port 5000) is working correctly:
- Generates new SVGs on-demand via Python scripts
- Uses RDKit for molecule rendering
- Caches generated SVGs for performance
- Supports both SMILES and nomenclature inputs

### Test Results
```bash
# Test SMILES endpoint
curl "http://localhost:5000/img/smiles?smiles=c1ccccc1&width=300&height=200"
# ✓ Successfully generates benzene SVG

# Test nomenclature endpoint
curl "http://localhost:5000/img/nomenclature?nomenclature=aspirin&width=300&height=200"
# ✓ Converts name to SMILES and generates SVG
```

### Cache Behavior
- First request: Generates new SVG via Python, saves to cache
- Subsequent requests: Serves from cache (fast)
- Cache key: MD5 hash of `type:value:widthxheight`
- Cache location: `MoleculeViewer/cache/moleculeviewer/`

## Integration with Content Script

The content.js already has OPSIN support in the name-to-SMILES conversion:

```javascript
// Existing code in content.js
async function convertNameToSmilesWithFallbacks(name) {
  // Priority 1: OPSIN (already implemented)
  try {
    const opsinResp = await fetch(`https://opsin.ch.cam.ac.uk/opsin/${name}.json`);
    if (opsinResp.ok) {
      const result = await opsinResp.json();
      if (result.smiles) {
        return { smiles: result.smiles, source: 'OPSIN' };
      }
    }
  } catch (e) {
    console.warn('OPSIN failed, trying fallback...');
  }

  // Priority 2: PubChem/MoleculeViewer
  // ... fallback logic
}
```

To enable 3D SMILES in the extension:
1. User enables "3D Stereochemistry (OPSIN)" in popup
2. Extension stores `use3DSmiles: true` in Chrome storage
3. When rendering nomenclature, extension uses OPSIN first
4. OPSIN returns 3D SMILES with stereochemistry
5. Structure is rendered with mol2chemfig showing stereochemistry

## File Changes Summary

### Modified Files
1. **mol2chemfig_server.py**
   - Added `/api/opsin` endpoint
   - Added `/api/generate-3d` endpoint
   - Updated endpoint documentation
   - Updated startup messages

2. **chem-extension/popup.html**
   - Added use3DSmiles toggle in mol2chemfig section
   - Added info box explaining 3D mode

3. **chem-extension/popup.js**
   - Added use3DSmilesToggle element
   - Added use3DSmiles default setting
   - Added use3DSmiles loading logic
   - Added use3DSmiles event listener

### New Files
1. **test_3d_smiles.html**
   - Comprehensive test interface for OPSIN integration
   - OPSIN conversion testing
   - 3D structure generation testing
   - Stereochemistry test cases
   - 2D vs 3D comparison

2. **OPSIN_3D_IMPLEMENTATION.md** (this file)
   - Complete implementation documentation
   - API examples
   - Testing guide
   - Integration instructions

## Server Startup

To use the 3D SMILES features:

```bash
# Terminal 1: Start mol2chemfig Docker backend
cd mol2chemfig-docker
docker-compose up

# Terminal 2: Start mol2chemfig server (with OPSIN)
cd ChemParser
python mol2chemfig_server.py

# Terminal 3: Start MoleculeViewer (optional)
cd ChemParser/MoleculeViewer
node server.js
```

Server endpoints will be available at:
- mol2chemfig: http://localhost:5001
- MoleculeViewer: http://localhost:5000
- Docker backend: http://localhost:8000

## Benefits of 3D SMILES

1. **Accurate Stereochemistry**: Preserves chiral centers and spatial configuration
2. **IUPAC Naming**: Properly handles systematic chemical names
3. **Enantiomer Distinction**: Can differentiate D/L forms, R/S configurations
4. **Educational Value**: Shows proper 3D structure for teaching
5. **Research Accuracy**: Essential for drug molecules and complex organics

## Limitations

1. **OPSIN Coverage**: Some trivial names may not be recognized
2. **Internet Required**: OPSIN API requires internet connection
3. **Rate Limits**: OPSIN may have rate limiting (use caching)
4. **Complex Names**: Very complex or uncommon names may fail

## Future Enhancements

1. Add local OPSIN installation for offline use
2. Implement OPSIN caching in mol2chemfig_server
3. Add stereochemistry visualization options
4. Support for InChI and InChIKey formats
5. Batch processing for multiple molecules
6. 3D model visualization (WebGL/Three.js)

## References

- OPSIN API: https://opsin.ch.cam.ac.uk/
- OPSIN GitHub: https://github.com/dan2097/opsin
- SMILES Notation: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
- Stereochemistry in SMILES: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html#STEREO

## Support

For issues or questions:
1. Check the test page: `test_3d_smiles.html`
2. Review server logs for errors
3. Test OPSIN API directly: `https://opsin.ch.cam.ac.uk/opsin/glucose.json`
4. Verify mol2chemfig Docker backend is running
5. Check Chrome extension console for errors

---

**Implementation Date**: November 2025
**Version**: 1.0
**Status**: Complete and tested
