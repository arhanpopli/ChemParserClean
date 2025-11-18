# MoleculeViewer - Deployment Package

## 🚀 Quick Start (30 seconds)

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. (Optional) Setup OPSIN parser for IUPAC names
# Windows:
setup_opsin.bat
# Linux/Mac:
python setup_opsin.py

# 3. Start server
python run_server.py

# 4. Open browser
http://localhost:5000
```

Server will be running at **http://localhost:5000**

---

## 📋 What This Is

MoleculeViewer is a web-based chemical molecule visualization tool that:

- **Converts SMILES to SVG** - Renders chemical structures as interactive SVG
- **IUPAC Nomenclature** - Converts chemical names to SMILES (aspirin → structure)
- **Carbon Labels** - Shows CH₃, CH₂, CH, C labels for methyl groups
- **Aromatic Circles** - Displays aromatic rings with indicator circles (benzene)
- **Multiple Parsers** - ChemDoodle, OPSIN, PubChem, Fallback dictionary

---

## 🔧 Installation

### Requirements

- **Python 3.8+** 
- **pip** (Python package manager)
- **Java** (optional, for OPSIN IUPAC parser)

### Step 1: Install Python Dependencies

```bash
pip install -r requirements.txt
```

**What gets installed:**
- `Flask` - Web server
- `RDKit` - Molecular structure rendering
- `flask-cors` - Cross-origin requests

### Step 2: (Optional) Setup OPSIN Parser

For converting IUPAC names to SMILES (systematic nomenclature):

**Windows:**
```bash
setup_opsin.bat
```

**Linux/Mac:**
```bash
python setup_opsin.py
```

This downloads and configures the OPSIN IUPAC name parser.

### Step 3: Verify Installation

```bash
python -c "from app.chemistry import smiles_to_svg; print('✓ Installed successfully')"
```

Should print: `✓ Installed successfully`

---

## 🎯 Quick Usage

### Start Server

```bash
python run_server.py
```

Output:
```
==================================================
  MoleculeViewer Flask Server
==================================================
Server available at:
  http://localhost:5000
  http://127.0.0.1:5000
Press CTRL+C to stop
```

### Test It

Open browser: **http://localhost:5000**

Try entering:
- **SMILES**: `CC(C)C` (isobutane)
- **Chemical Name**: `aspirin`

### Stop Server

Press **Ctrl+C** in terminal

---

## ⚙️ Configuration

### Change Port

Edit `run_server.py`:
```python
app.run(host='0.0.0.0', port=8000, debug=False)  # Change 5000 to 8000
```

### Change Parser

Edit `app/config.py`:
```python
NOMENCLATURE_PARSER = 'auto'  # 'chemdoodle', 'opsin', 'pubchem', 'fallback'
```

### Adjust CH₃ Label Size

Edit `app/config.py`:
```python
CARBON_LABEL_FONT_SIZE = 32  # Increase or decrease (24-48 recommended)
CARBON_LABEL_SCALING = 'auto'  # 'auto' or 'fixed'
```

---

## 📁 Project Structure

```
MoleculeViewer/
├── app/                      # Main application code
│   ├── __init__.py
│   ├── api.py               # REST API endpoints
│   ├── chemistry.py         # Molecule processing
│   ├── config.py            # Configuration settings
│   └── chemdoodle_compounds.py  # Compound database
├── static/                  # Frontend (CSS, JS)
├── templates/               # HTML templates
├── requirements.txt         # Python dependencies
├── run_server.py           # Start server (USE THIS)
├── opsin-cli.jar           # OPSIN parser (optional)
└── app/config.py           # Settings
```

---

## 🌐 API Endpoints

### Convert SMILES to SVG

```javascript
POST http://localhost:5000/api/smiles-to-svg
Content-Type: application/json

{
  "smiles": "CC(C)C",
  "width": 400,
  "height": 300,
  "options": {
    "show_methyls": true,
    "aromatic_circles": true
  }
}
```

### Convert Nomenclature to SMILES

```javascript
POST http://localhost:5000/api/nomenclature-to-smiles
Content-Type: application/json

{
  "nomenclature": "aspirin"
}
```

Response:
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "source": "ChemDoodle Database"
}
```

---

## 🔍 Features

### Show Methyls (CH₃ Labels)
Option: `show_methyls: true`
Shows `CH₃` labels for terminal methyl groups

### Show All Carbons
Option: `show_carbons: true`  
Labels: `CH₃`, `CH₂`, `CH`, `C`

### Aromatic Circles
Option: `aromatic_circles: true`
Shows dashed circles in aromatic rings (benzene visualization)

### Explicit Hydrogens
Option: `explicit_hydrogens: true/false/default`
- `true` - Show all H atoms
- `false` - Hide all H atoms  
- `default` - RDKit default

---

## 🚨 Troubleshooting

### "ModuleNotFoundError: No module named 'rdkit'"

**Solution**: Install dependencies
```bash
pip install -r requirements.txt
```

### "Port 5000 already in use"

**Solution**: Use different port in `run_server.py`:
```python
app.run(port=8001)  # Use 8001 instead
```

### "OPSIN parser not working"

**Solution**: 
1. Verify Java is installed: `java -version`
2. Run setup: `python setup_opsin.py`
3. Check `opsin-cli.jar` exists

### "Nomenclature lookup fails"

**Solution**: Check parser in `app/config.py`:
```python
NOMENCLATURE_PARSER = 'auto'  # Try all parsers
ENABLE_PUBCHEM = True  # Enable internet lookup
```

---

## 📝 Configuration Files

### `app/config.py` - Main Settings

**Carbon Labels:**
- `CARBON_LABEL_FONT_SIZE` - Label size (default: 32)
- `CARBON_LABEL_SCALING` - 'auto' or 'fixed'
- `CARBON_LABEL_SCALE_FACTOR` - Size multiplier (0.55)

**Nomenclature Parsers:**
- `NOMENCLATURE_PARSER` - Which parser to use
- `ENABLE_CHEMDOODLE` - Use ChemDoodle database
- `ENABLE_OPSIN` - Use OPSIN IUPAC parser
- `ENABLE_PUBCHEM` - Use PubChem API
- `PARSER_TIMEOUT` - Seconds to wait

**Aromatic Circles:**
- `AROMATIC_CIRCLE_RADIUS_6` - 6-membered ring radius
- `AROMATIC_CIRCLE_RADIUS_5` - 5-membered ring radius

---

## 🧪 Testing

### Test Nomenclature Parser

```bash
python test_parser_config.py
```

Tests various compound names with your parser configuration.

### Test SMILES Rendering

```python
python -c "
import sys
sys.path.insert(0, '.')
from app.chemistry import smiles_to_svg

error, svg = smiles_to_svg('CC(C)C', 400, 300, {'show_methyls': True})
print('✓ Works!' if not error else f'✗ Error: {error}')
"
```

---

## 🌍 Deployment Scenarios

### Local Development
```bash
python run_server.py
# Access: http://localhost:5000
```

### LAN Network
```bash
python run_server.py
# Access: http://<your-ip>:5000
```

### Production (gunicorn)
```bash
pip install gunicorn
gunicorn --workers 4 --bind 0.0.0.0:5000 app.api:app
```

### Docker
```bash
docker build -t moleculeviewer .
docker run -p 5000:5000 moleculeviewer
```

---

## 📚 Additional Documentation

- **`PARSER_CONFIG_GUIDE.md`** - Detailed parser configuration
- **`CONFIG_GUIDE.md`** - Carbon label sizing guide
- **`QUICKSTART.md`** - Features and examples
- **`START_HERE.md`** - Feature overview

---

## 📞 Support

### Check Dependencies
```bash
python -c "import flask, rdkit, rdkit.Chem; print('All dependencies installed')"
```

### View Server Logs
The server prints all requests to console while running.

### Test Connectivity
```bash
curl http://localhost:5000/
```

### Enable Debug Mode (Development Only)
Edit `run_server.py`:
```python
app.run(debug=True)  # Enables hot reload
```

---

## ✅ Checklist for Deployment

- [ ] Python 3.8+ installed
- [ ] Dependencies installed: `pip install -r requirements.txt`
- [ ] (Optional) OPSIN setup completed
- [ ] Server starts: `python run_server.py`
- [ ] Accessible: `http://localhost:5000`
- [ ] Basic test works: Enter SMILES or compound name
- [ ] Configuration reviewed: `app/config.py`

---

**Ready to deploy!** 🚀

For detailed configuration options, see `app/config.py` comments.
