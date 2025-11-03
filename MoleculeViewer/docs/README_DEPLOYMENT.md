# MoleculeViewer - Deployment Ready Package

**Status: ✅ READY FOR SERVER DEPLOYMENT**

This package contains everything needed to run MoleculeViewer as a standalone web application on a server.

---

## 📦 What's in the Package

### Core Application
- **`app/`** - Complete Flask web application
  - `__init__.py` - Flask app initialization
  - `api.py` - REST API endpoints
  - `chemistry.py` - Molecule rendering and parsing logic
  - `config.py` - Centralized configuration
  - `chemdoodle_compounds.py` - Compound database

- **`static/`** - Frontend assets
  - `index.html` - Web interface
  - CSS stylesheets
  - JavaScript files

- **`templates/`** - HTML templates

### Entry Points
- **`run_server.py`** - Main server startup script ⭐ USE THIS

### Dependencies
- **`requirements.txt`** - All Python packages with exact versions
  - Flask==2.3.0
  - rdkit==2024.9.1
  - Werkzeug==2.3.0
  - flask-cors==4.0.0

### Setup Tools
- **`setup_opsin.bat`** - Windows OPSIN installer
- **`setup_opsin.py`** - Linux/Mac/Windows OPSIN installer
- **`opsin-cli.jar`** - OPSIN parser (optional, improves IUPAC nomenclature support)

### Documentation
- **`DEPLOYMENT.md`** - Complete setup and configuration guide ⭐ START HERE
- **`PACKAGE_CHECKLIST.md`** - Verification checklist
- **`CONFIG_GUIDE.md`** - Carbon label customization
- **`PARSER_CONFIG_GUIDE.md`** - Parser configuration options
- **`README.md`** - Project overview

---

## 🚀 30-Second Quick Start

### Step 1: Install Dependencies
```bash
pip install -r requirements.txt
```

### Step 2: Start Server
```bash
python run_server.py
```

### Step 3: Open Browser
```
http://localhost:5000
```

**Done!** The web interface is now running.

---

## 📋 Complete Setup (with OPSIN for IUPAC names)

### Step 1: Install Python Dependencies
```bash
pip install -r requirements.txt
```

### Step 2: Install Java (Required only if using OPSIN)
- Windows: https://java.com/en/download/
- Linux: `sudo apt install default-jre`
- Mac: `brew install openjdk`

### Step 3: Setup OPSIN (Optional but Recommended)
```bash
# Windows:
python setup_opsin.bat

# Linux/Mac:
python setup_opsin.py
```

### Step 4: Start Server
```bash
python run_server.py
```

### Step 5: Access Web Interface
```
http://localhost:5000
```

---

## ⚙️ Configuration

All settings are in **`app/config.py`**. Key options:

### Carbon Label Sizing
```python
CARBON_LABEL_FONT_SIZE = 32        # Font size in pixels
CARBON_LABEL_SCALING = 'auto'      # 'auto' or 'fixed'
CARBON_LABEL_SCALE_FACTOR = 0.55   # Scale multiplier for auto mode
```

### Nomenclature Parsers
```python
NOMENCLATURE_PARSER = 'auto'  # Options: 'auto', 'chemdoodle', 'opsin', 'pubchem', 'fallback'
ENABLE_CHEMDOODLE = True      # Fast, basic names
ENABLE_OPSIN = True           # IUPAC names (requires Java)
ENABLE_FALLBACK = True        # Simple fallback
ENABLE_PUBCHEM = True         # Comprehensive database
```

### Server Settings
```python
DEBUG = False                 # Set to True for development
HOST = '0.0.0.0'             # Accessible from network
PORT = 5000                  # Change if needed
```

See `app/config.py` for all available options.

---

## 🌐 API Reference

### GET `/` 
Returns the web interface (HTML page)

### POST `/api/render`
Render a molecule structure
```json
{
  "smiles": "CCO"
}
```

### POST `/api/render-name`
Render from compound name
```json
{
  "name": "ethanol"
}
```

### POST `/api/name-to-smiles`
Convert compound name to SMILES
```json
{
  "name": "ethanol"
}
```

See `DEPLOYMENT.md` for complete API documentation.

---

## 🔧 Deployment Scenarios

### Scenario 1: Development Machine
```bash
python run_server.py
```
Access at: `http://localhost:5000`

### Scenario 2: Production Server
```bash
# Using Gunicorn (better for production)
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 app:app

# Or basic Flask (for testing)
python run_server.py
```

### Scenario 3: Docker Container
```bash
docker-compose up
```

### Scenario 4: Different Port
Edit `app/config.py`:
```python
PORT = 8080  # Change from 5000 to 8080
```

---

## ✅ Verification Checklist

Before deployment, verify:

- [ ] Python 3.8+ installed: `python --version`
- [ ] Dependencies installed: `pip list | grep -i flask`
- [ ] Server starts: `python run_server.py` (should show "Running on http://localhost:5000")
- [ ] Web interface loads: Open `http://localhost:5000` in browser
- [ ] Basic test works: Enter "ethanol" in search box
- [ ] (Optional) Java installed: `java -version` (if using OPSIN)
- [ ] (Optional) OPSIN works: `python setup_opsin.py` runs successfully

---

## 🐛 Troubleshooting

### "ModuleNotFoundError: No module named 'flask'"
**Solution:** Run `pip install -r requirements.txt`

### Port 5000 Already in Use
**Solution:** 
1. Use different port in `app/config.py`: `PORT = 8080`
2. Or kill existing process: `lsof -ti:5000 | xargs kill` (Linux/Mac)

### OPSIN Not Working
**Solution:**
1. Verify Java: `java -version`
2. Re-run setup: `python setup_opsin.py`
3. Or disable: Set `ENABLE_OPSIN = False` in `app/config.py`

### Molecule Visualization Issues
**Solution:** Check `app/config.py` settings:
- Adjust `CARBON_LABEL_FONT_SIZE` (try 28-36)
- Change `CARBON_LABEL_SCALING` ('auto' vs 'fixed')

See `DEPLOYMENT.md` for more troubleshooting help.

---

## 📚 Documentation Structure

| Document | Purpose | Audience |
|----------|---------|----------|
| **DEPLOYMENT.md** | Complete setup guide | Everyone (start here) |
| **CONFIG_GUIDE.md** | Carbon label customization | Customization needed |
| **PARSER_CONFIG_GUIDE.md** | Parser configuration | Advanced users |
| **PACKAGE_CHECKLIST.md** | Verification checklist | Before deployment |
| **README.md** | Project overview | New to system |

---

## 🎯 Next Steps

1. **First Time?**
   - Read `DEPLOYMENT.md` for complete setup
   - Follow Quick Start above

2. **Want to Customize?**
   - See `CONFIG_GUIDE.md` for label sizing
   - See `PARSER_CONFIG_GUIDE.md` for parser options

3. **Deploying to Server?**
   - Follow "Complete Setup" above
   - See "Deployment Scenarios" in `DEPLOYMENT.md`

4. **Something Not Working?**
   - Check "Troubleshooting" section below
   - See `DEPLOYMENT.md` Troubleshooting section

---

## 📞 Support

If something doesn't work:

1. **Check logs:** Server output shows error messages
2. **Verify config:** Check `app/config.py` settings
3. **Try defaults:** Reset settings to defaults in `app/config.py`
4. **Check documentation:** `DEPLOYMENT.md` has extensive troubleshooting

---

## 📄 Package Contents Summary

```
MoleculeViewer/
├── app/                           # Flask application
│   ├── __init__.py
│   ├── api.py                     # REST endpoints
│   ├── chemistry.py               # Rendering logic
│   ├── config.py                  # ⭐ All settings here
│   └── chemdoodle_compounds.py    # Compound database
│
├── static/                        # Frontend assets
│   └── [HTML, CSS, JS files]
│
├── templates/                     # Templates
│   └── [HTML templates]
│
├── run_server.py                  # ⭐ Start with this
├── requirements.txt               # Dependencies
├── opsin-cli.jar                  # OPSIN parser
├── setup_opsin.bat                # Windows OPSIN setup
├── setup_opsin.py                 # Linux/Mac/Windows OPSIN setup
│
└── [Documentation files]
    ├── DEPLOYMENT.md              # ⭐ Complete guide
    ├── PACKAGE_CHECKLIST.md
    ├── CONFIG_GUIDE.md
    ├── PARSER_CONFIG_GUIDE.md
    └── README.md
```

---

## ✨ Features Included

✅ **Molecule Visualization**
- SVG-based rendering
- Interactive structure display
- Support for benzene rings and functional groups

✅ **Multiple Input Methods**
- SMILES notation
- IUPAC nomenclature (requires OPSIN)
- Common compound names
- ChemDoodle database lookup

✅ **Configurable Rendering**
- Carbon label sizing (fixed or auto-scaling)
- Aromatic circle positioning
- Font customization

✅ **REST API**
- JSON endpoints
- CORS support for cross-origin requests
- Extensible architecture

✅ **Multi-Parser Support**
- ChemDoodle (fast)
- OPSIN (IUPAC names)
- PubChem (comprehensive)
- Fallback (basic)

---

## 🎉 You're Ready!

This package is **production-ready**. Everything needed is included:

- ✅ All Python code
- ✅ All dependencies specified
- ✅ Configuration system ready
- ✅ Setup scripts included
- ✅ Comprehensive documentation
- ✅ API reference included
- ✅ Troubleshooting guide included

**Next step:** Run `python run_server.py` and enjoy! 🚀

---

**Last Updated:** 2024
**Status:** Ready for Deployment ✅
**Version:** 2.0 (Production Ready)
