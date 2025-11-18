# MoleculeViewer - Standalone Web Application

Ready-to-deploy molecular structure viewer with support for SMILES, compound names, and IUPAC nomenclature.

## Quick Start (30 seconds)

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Start server
python run_server.py

# 3. Open browser
# http://localhost:5000
```

## Documentation

- **`README_DEPLOYMENT.md`** - Start here! Complete setup and configuration guide
- **`DEPLOYMENT.md`** - Detailed deployment guide with troubleshooting
- **`CONFIG_GUIDE.md`** - Carbon label sizing options
- **`PARSER_CONFIG_GUIDE.md`** - Nomenclature parser configuration
- **`PACKAGE_CHECKLIST.md`** - Deployment verification checklist

## What's Included

✅ Flask web server  
✅ Molecule visualization with RDKit  
✅ Multiple nomenclature parsers (ChemDoodle, OPSIN, PubChem, Fallback)  
✅ Auto-scaling labels  
✅ REST API endpoints  
✅ Configuration system  
✅ Setup scripts for optional OPSIN  

## Requirements

- Python 3.8+
- Dependencies in `requirements.txt` (Flask, RDKit, CORS)
- (Optional) Java 8+ for OPSIN IUPAC parser support

## Setup with OPSIN (5 minutes)

```bash
# Install dependencies
pip install -r requirements.txt

# Install Java (if needed)
# https://java.com/en/download/

# Setup OPSIN
python setup_opsin.py

# Start server
python run_server.py
```

## Configuration

All settings in `app/config.py`:
- Carbon label sizing
- Nomenclature parser selection
- Server host/port
- Aromatic circle styling

See `CONFIG_GUIDE.md` for details.

## API Endpoints

- `GET /` - Web interface
- `POST /api/render` - Render SMILES structure
- `POST /api/render-name` - Render from compound name
- `POST /api/name-to-smiles` - Convert name to SMILES

See `README_DEPLOYMENT.md` for full API documentation.

## Deployment

### Development
```bash
python run_server.py
```

### Production
```bash
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 app:app
```

### Docker
```bash
docker-compose up
```

## Troubleshooting

**"ModuleNotFoundError"**
```bash
pip install -r requirements.txt
```

**Port already in use**
- Change PORT in `app/config.py` or kill existing process

**OPSIN not working**
- Verify Java: `java -version`
- Re-run setup: `python setup_opsin.py`
- Or disable: `ENABLE_OPSIN = False` in `app/config.py`

See `DEPLOYMENT.md` for more help.

## Features

- **Multiple Input Formats:** SMILES, IUPAC names, common names
- **Smart Rendering:** Auto-scaling labels, proper arene ring display
- **Configurable:** Customize appearance and behavior
- **REST API:** JSON endpoints for integration
- **Multi-Parser:** Different parsers for different name formats

## File Structure

```
MoleculeViewer/
├── app/                 # Flask application
├── static/              # Frontend (HTML, CSS, JS)
├── templates/           # HTML templates
├── run_server.py        # Start server ← USE THIS
├── requirements.txt     # Dependencies
├── setup_opsin.py       # Optional OPSIN setup
└── [documentation]      # Guides and references
```

## Next Steps

1. **New to this?** Read `README_DEPLOYMENT.md`
2. **Want to customize?** Read `CONFIG_GUIDE.md`
3. **Need API docs?** See `README_DEPLOYMENT.md` API section
4. **Deploying to server?** See `DEPLOYMENT.md`

---

**Status:** ✅ Ready for deployment  
**All dependencies included:** ✅  
**Documentation:** ✅ Complete  
**Setup scripts:** ✅ Included
