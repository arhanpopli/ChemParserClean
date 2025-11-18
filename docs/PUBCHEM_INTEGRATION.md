# PubChem Integration Documentation

## Overview

PubChem integration adds a third rendering engine to the Chemistry Extension, allowing users to fetch molecular structure images and 3D models directly from the PubChem database.

## Features

### 1. Direct PubChem Image Fetching
- Fetches 2D structure images from PubChem API
- Supports chemical nomenclature (e.g., "histamine", "caffeine")
- Supports SMILES notation (e.g., "CCO", "c1ccccc1")
- Configurable image sizes: small (100x100), large (300x300), or custom

### 2. 3D Model Support
- Optional "3D View" button on each molecule image
- Opens PubChem's 3D viewer in a new window
- Displays 3D conformer models with rotation capabilities
- Downloads 3D structure data in SDF format

### 3. Caching System
- Local caching of fetched images
- Separate cache for 2D images and 3D models
- Reduces API calls and improves performance

## Installation

### 1. Install Python Dependencies

```bash
pip install -r requirements_pubchem.txt
```

Required packages:
- Flask >= 2.3.0
- flask-cors >= 4.0.0
- requests >= 2.31.0

### 2. Start the PubChem Server

**Windows:**
```bash
start_pubchem.bat
```

**Linux/Mac:**
```bash
python pubchem_server.py
```

The server will start on `http://localhost:5002`

### 3. Configure Chrome Extension

1. Open the Chrome extension settings (click the extension icon)
2. Under "Rendering Engine", select **🌐 PubChem**
3. Configure PubChem options:
   - **Image Size**: Choose small, large, or extra large
   - **Enable 3D Models**: Toggle to show/hide 3D view buttons
   - **Default View Type**: 2D structure or 3D model
4. Reload the page to apply changes

## Usage

### Basic Usage

Use the same syntax as other rendering engines:

```
chem:histamine:
chem:caffeine:
chem:CCO:
chem:c1ccccc1:
```

### 3D Viewer

When 3D models are enabled:
1. A "🔮 3D" button appears in the top-right corner of each molecule image
2. Click the button to open PubChem's 3D viewer in a new window
3. The 3D viewer allows rotation, zooming, and viewing different conformers

## API Endpoints

### PubChem Server Endpoints

#### 1. Direct Image (PNG)
```
GET /pubchem/img/{name}?size={size}&type={type}
```

**Parameters:**
- `name`: Chemical name or SMILES (required)
- `size`: Image size - `small`, `large`, or custom like `500x500` (default: `large`)
- `type`: Record type - `2d` or `3d` (default: `2d`)

**Example:**
```
http://localhost:5002/pubchem/img/histamine?size=large&type=2d
```

#### 2. Image Info (JSON)
```
GET /pubchem/image?name={name}&size={size}&type={type}
```

Returns JSON with image URL, CID, and metadata.

#### 3. 3D Model (SDF)
```
GET /pubchem/3d-model?name={name}&format={format}&conformer={index}
```

**Parameters:**
- `name`: Chemical name or SMILES (required)
- `format`: `sdf` or `json` (default: `sdf`)
- `conformer`: Conformer index, 0 for default (default: `0`)

**Example:**
```
http://localhost:5002/pubchem/3d-model?name=histamine&format=sdf
```

#### 4. 3D Viewer
```
GET /pubchem/3d-viewer?name={name}&embed={true/false}
```

**Parameters:**
- `name`: Chemical name or SMILES (required)
- `embed`: `true` for HTML page, `false` for redirect (default: `false`)

**Example:**
```
http://localhost:5002/pubchem/3d-viewer?name=histamine&embed=true
```

#### 5. Compound Info
```
GET /pubchem/info?name={name}
```

Returns comprehensive information about the compound including CID, URLs, and available conformers.

#### 6. Cache Management
```
GET /cache-info
DELETE /clear-cache
```

### PubChem REST API (Direct)

The server uses PubChem's PUG REST API:

**Image:**
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/PNG?image_size=large
```

**3D Structure (SDF):**
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF?record_type=3d
```

**Conformers:**
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/conformers/TXT
```

## Architecture

### Server Structure (`pubchem_server.py`)

```
pubchem_server.py
├── Helper Functions
│   ├── get_compound_cid()      # Convert name/SMILES to PubChem CID
│   ├── fetch_pubchem_image()   # Fetch PNG from PubChem
│   ├── fetch_pubchem_sdf()     # Fetch 3D model (SDF format)
│   └── get_conformers()        # Get list of 3D conformers
│
├── Routes
│   ├── /                       # API info page
│   ├── /health                 # Health check
│   ├── /pubchem/image          # Image info (JSON)
│   ├── /pubchem/img/<name>     # Direct PNG image
│   ├── /pubchem/3d-model       # 3D model (SDF/JSON)
│   ├── /pubchem/3d-viewer      # 3D viewer page
│   ├── /pubchem/info           # Compound info
│   ├── /cache-info             # Cache statistics
│   └── /clear-cache            # Clear cache
│
└── Cache System
    ├── pubchem-cache/images/   # Cached PNG images
    └── pubchem-cache/sdf/      # Cached 3D models
```

### Extension Integration

1. **Content Script (`content.js`)**
   - Added `PUBCHEM_API` constant
   - Added `loadPubChemImage()` function
   - Added `add3DViewButton()` function
   - Modified `loadMoleculeImage()` to support PubChem mode

2. **Popup UI (`popup.html` & `popup.js`)**
   - Added PubChem radio button option
   - Added PubChem options panel
   - Added settings for image size, 3D enable, and record type

## Testing

### Test File

A comprehensive test file is provided: `test_pubchem.html`

**To run tests:**
1. Start the PubChem server: `start_pubchem.bat`
2. Open `test_pubchem.html` in Chrome
3. Enable the extension and select PubChem mode
4. Reload the page
5. Verify all molecules render correctly
6. Test 3D viewer buttons

### Test Coverage

The test file includes:
- Chemical nomenclature (histamine, caffeine, aspirin, etc.)
- SMILES notation (CCO, c1ccccc1, etc.)
- Complex molecules (adrenaline, penicillin, morphine)
- Neurotransmitters (serotonin, GABA, dopamine)

## Advantages

### Compared to Local Rendering (MoleculeViewer)
- ✅ No need for RDKit installation
- ✅ Access to PubChem's extensive database
- ✅ 3D model support out of the box
- ✅ Professional quality images
- ❌ Requires internet connection
- ❌ Limited customization options

### Compared to mol2chemfig
- ✅ Simpler setup (no Docker required)
- ✅ 3D model support
- ✅ Faster image loading
- ❌ No LaTeX/ChemFig output
- ❌ No custom rendering options

## Troubleshooting

### Server Won't Start
- Check if port 5002 is available
- Verify Python dependencies are installed
- Check `pip install -r requirements_pubchem.txt`

### Images Not Loading
- Verify PubChem server is running on port 5002
- Check browser console for error messages
- Test direct API access: `http://localhost:5002/pubchem/img/histamine`
- Verify internet connection (server needs to reach PubChem API)

### 3D Button Not Appearing
- Check extension settings: "Enable 3D Models" should be ON
- Reload the page after changing settings
- Verify PubChem mode is selected

### Compound Not Found
- PubChem may not have the compound in its database
- Try using SMILES notation instead of nomenclature
- Check spelling of chemical names
- Some compounds may only have 2D structures (no 3D models)

## Cache System

### Cache Location
```
pubchem-cache/
├── images/     # PNG images (2D and 3D renders)
└── sdf/        # 3D structure files
```

### Cache Key Format
- Images: MD5 hash of `img_{cid}_{size}_{type}`
- SDF files: MD5 hash of `sdf_{cid}_{conformer_index}`

### Clear Cache
```bash
curl -X DELETE http://localhost:5002/clear-cache
```

Or use the `/cache-info` endpoint to view cache statistics.

## Future Enhancements

Potential improvements for future versions:

1. **Multiple Conformer Support**
   - Display multiple 3D conformers
   - Allow user to select preferred conformer

2. **Embedded 3D Viewer**
   - Integrate 3D viewer directly in page
   - Use WebGL-based molecule viewer (e.g., 3Dmol.js, NGL Viewer)

3. **Batch Loading**
   - Preload multiple molecules at once
   - Reduce API call overhead

4. **Property Display**
   - Show molecular weight, formula
   - Display chemical properties from PubChem

5. **Download Options**
   - Download as PNG, SVG, or SDF
   - Export to different formats

6. **Similarity Search**
   - Find similar compounds
   - Display alternative structures

## References

- [PubChem PUG REST API Documentation](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [PubChem 3D Viewer](https://pubchem.ncbi.nlm.nih.gov/docs/3d-structure-viewer)
- [PubChem Database](https://pubchem.ncbi.nlm.nih.gov/)

## License

Part of the Chemparser project. See main project README for license information.

## Support

For issues, questions, or contributions, refer to the main Chemparser project documentation.
