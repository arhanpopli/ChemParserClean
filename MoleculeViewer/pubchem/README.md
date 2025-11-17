# 🌐 PubChem Integration Server

A lightweight Node.js server that provides direct access to PubChem molecular structures and 3D viewers.

## 🚀 Features

- **2D Structure Images**: Direct image links from PubChem (no local hosting required)
- **3D Interactive Viewers**: Beautiful 3D molecular viewers with controls
- **Smart Caching**: CID (Compound ID) lookups are cached for performance
- **CORS Enabled**: Works seamlessly with browser extensions
- **Multiple Formats**: Support for PNG images, SDF files, and embedded viewers

## 📦 Installation

```bash
# Navigate to the pubchem directory
cd MoleculeViewer/pubchem

# Install dependencies
npm install
```

## ▶️ Starting the Server

### Windows
```bash
start.bat
```

### macOS/Linux
```bash
npm start
```

The server will start on **http://localhost:5002**

## 🔗 API Endpoints

### 1. Get 2D Image (Direct PNG)
Returns the actual PNG image data from PubChem.

```
GET http://localhost:5002/pubchem/img/{name}
```

**Parameters:**
- `name` - Compound name or SMILES (required)
- `size` - Image size: `small`, `large`, or custom like `500x500` (optional, default: `large`)
- `type` - Record type: `2d` or `3d` (optional, default: `2d`)

**Example:**
```html
<img src="http://localhost:5002/pubchem/img/histamine">
<img src="http://localhost:5002/pubchem/img/caffeine?size=500x500">
```

### 2. Get Image Info
Returns JSON with image URLs and metadata.

```
GET http://localhost:5002/pubchem/image?name={name}
```

**Example Response:**
```json
{
  "cid": 774,
  "name": "histamine",
  "image_url": "http://localhost:5002/pubchem/img/histamine?size=large&type=2d",
  "direct_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/774/PNG?image_size=large&record_type=2d",
  "cached": true
}
```

### 3. Get 3D Model (SDF File)
Returns the 3D structure data file.

```
GET http://localhost:5002/pubchem/3d-model?name={name}
```

**Parameters:**
- `name` or `cid` - Compound identifier (required)
- `format` - `sdf` or `json` (optional, default: `sdf`)
- `conformer` - Conformer index (optional)

### 4. Open 3D Viewer
Opens an interactive 3D molecular viewer.

```
GET http://localhost:5002/pubchem/3d-viewer?name={name}&embed=true
```

**Parameters:**
- `name` or `cid` - Compound identifier (required)
- `embed` - `true` for embedded viewer, `false` to redirect to PubChem (optional, default: `false`)

**Example:**
```javascript
// Open in new window
window.open('http://localhost:5002/pubchem/3d-viewer?name=histamine&embed=true', 
            '_blank', 'width=1000,height=700');
```

**3D Viewer Features:**
- 🎨 Multiple rendering styles (Ball & Stick, Stick, Space Filling, Wireframe)
- 💫 Auto-rotation option
- 🔬 Show/hide hydrogens
- 🔗 Direct link to PubChem
- 💾 Download SDF file

### 5. Get Compound Info
Returns comprehensive information about a compound.

```
GET http://localhost:5002/pubchem/info?name={name}
```

**Example Response:**
```json
{
  "cid": 774,
  "name": "histamine",
  "pubchem_url": "https://pubchem.ncbi.nlm.nih.gov/compound/774",
  "3d_viewer_url": "https://pubchem.ncbi.nlm.nih.gov/compound/774#section=3D-Conformer",
  "image_url_2d": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/774/PNG?image_size=large",
  "image_url_3d": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/774/PNG?image_size=large&record_type=3d",
  "sdf_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/774/record/SDF?record_type=3d",
  "conformers_count": 10,
  "local_endpoints": {
    "image": "http://localhost:5002/pubchem/img/histamine",
    "3d_model": "http://localhost:5002/pubchem/3d-model?cid=774",
    "3d_viewer": "http://localhost:5002/pubchem/3d-viewer?cid=774"
  }
}
```

## 🧪 Testing

Open `test_pubchem_3d.html` in your browser (located in the project root) to test all functionality:

```bash
# From project root
# Simply open test_pubchem_3d.html in your browser
```

## 🎨 Chrome Extension Integration

The extension automatically uses this server when the **PubChem** renderer is selected.

### Enable 3D Viewer in Extension:

1. Open extension popup
2. Go to **Developer Options** tab
3. Under **3D Viewer Mode**, enable **Enable 3D Viewer**
4. Reload any page with chemistry formulas
5. Look for the **🔮 3D** button on molecule images

## 📁 Directory Structure

```
pubchem/
├── server.js              # Main server file
├── package.json           # Dependencies
├── start.bat              # Windows startup script
├── static/                # Static files
│   └── structure-3d-webgl.min.js   # PubChem 3D library
└── README.md             # This file
```

## 💾 Caching

- **CID Lookups**: Cached in `../cache/pubchem/cid_cache.json`
- **Images**: Direct from PubChem (no local caching)
- **3D Models**: Not cached (fetched on demand)

This minimizes storage requirements while maintaining fast lookups.

## 🔧 Configuration

The server runs on **port 5002** by default. To change:

Edit `server.js`:
```javascript
const PORT = 5002; // Change this to your preferred port
```

## 🌟 Example Usage in Extension

When you write `chem:histamine:` on a webpage with the extension active:

1. Extension detects the formula
2. If **PubChem** renderer is selected:
   - Fetches 2D image from `http://localhost:5002/pubchem/img/histamine`
   - Displays the image inline
3. If **Enable 3D Viewer** is ON:
   - Shows a **🔮 3D** button
   - Clicking opens interactive 3D viewer

## 📚 Supported Compounds

Any compound available on PubChem can be used:
- **By Name**: `histamine`, `caffeine`, `aspirin`, `dopamine`, `glucose`, etc.
- **By SMILES**: `CCO` (ethanol), `c1ccccc1` (benzene), etc.

## 🐛 Troubleshooting

### Server won't start
```bash
# Make sure Node.js is installed
node --version

# Install dependencies
npm install

# Check if port 5002 is already in use
netstat -ano | findstr :5002
```

### Compound not found
- Check spelling of compound name
- Try using SMILES notation instead
- Verify compound exists on PubChem.ncbi.nlm.nih.gov

### 3D viewer not loading
- Check browser console for errors
- Ensure `embed=true` parameter is included
- Verify internet connection (viewer loads PubChem data)

## 📝 License

MIT

## 🤝 Contributing

This is part of the ChemParser project. Contributions welcome!

---

**Server Status:** ✅ Running on http://localhost:5002

Made with ❤️ for chemistry visualization
