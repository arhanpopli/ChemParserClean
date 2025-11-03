# 📑 MoleculeViewer Package - Complete File Index

## 🚀 Start Here

| File | Purpose | Read When |
|------|---------|-----------|
| **DELIVERY_SUMMARY.md** | Overview of everything included | First (this explains what you got) |
| **START_HERE_README.md** | Quick start guide | Second (before running anything) |
| **README_DEPLOYMENT.md** | Complete deployment guide | Before deploying to server |

## ▶️ To Run the Application

```bash
1. pip install -r requirements.txt          # Install dependencies
2. python setup_opsin.py                    # Optional: improve IUPAC parser
3. python run_server.py                     # Start the server
4. Open: http://localhost:5000              # View in browser
```

## 📚 Documentation (Read in Order of Need)

### Essential
| File | Content | For Whom |
|------|---------|----------|
| **START_HERE_README.md** | Quick overview + links | Everyone first |
| **README_DEPLOYMENT.md** | Detailed guide with API | Complete reference |
| **DEPLOYMENT.md** | Extended scenarios + troubleshooting | Server deployment |

### Customization
| File | Content | For Whom |
|------|---------|----------|
| **CONFIG_GUIDE.md** | Carbon label sizing options | UI customization |
| **PARSER_CONFIG_GUIDE.md** | Nomenclature parser setup | Advanced configuration |
| **PACKAGE_CHECKLIST.md** | Verification steps | Before deployment |

## 💻 Application Files

### Core Application
```
app/
├── __init__.py              # Flask app initialization
├── api.py                   # REST API endpoints
├── chemistry.py             # Molecule rendering engine ⭐
├── config.py                # Configuration hub ⭐⭐⭐ EDIT THIS
└── chemdoodle_compounds.py  # Compound database
```

### Frontend
```
static/                       # HTML, CSS, JavaScript files
templates/                    # Jinja2 HTML templates
run_server.py                 # Main entry point ⭐⭐ RUN THIS
```

## 📦 Dependencies & Setup

```
requirements.txt              # All Python packages (exact versions)
setup_opsin.py                # Cross-platform OPSIN installer (Python)
setup_opsin.bat               # Windows OPSIN installer (Batch)
opsin-cli.jar                 # OPSIN parser JAR (optional)
```

**Key Versions:**
- Flask==2.3.0
- rdkit==2024.9.1
- Werkzeug==2.3.0
- flask-cors==4.0.0

## 🔧 Configuration

**Main config file:**
- `app/config.py` ← Edit this to customize everything

**Key settings:**
```python
# Server
PORT = 5000

# Carbon labels
CARBON_LABEL_FONT_SIZE = 32
CARBON_LABEL_SCALING = 'auto'

# Parsers
NOMENCLATURE_PARSER = 'auto'
ENABLE_OPSIN = True
ENABLE_CHEMDOODLE = True
```

## 🎯 Common Tasks

### "Just want to run it"
1. `pip install -r requirements.txt`
2. `python run_server.py`
3. Open `http://localhost:5000`
4. Done! ✅

### "Want to customize labels"
1. Read: `CONFIG_GUIDE.md`
2. Edit: `app/config.py`
3. Change: `CARBON_LABEL_FONT_SIZE`, `CARBON_LABEL_SCALING`

### "Want to change parser"
1. Read: `PARSER_CONFIG_GUIDE.md`
2. Edit: `app/config.py`
3. Change: `NOMENCLATURE_PARSER` (options: 'auto', 'chemdoodle', 'opsin', 'pubchem', 'fallback')

### "Need to deploy to server"
1. Read: `DEPLOYMENT.md` → "Deployment Scenarios"
2. Or: `README_DEPLOYMENT.md` → "Deployment Scenarios"
3. Copy folder to server
4. Run same commands as local

### "Something isn't working"
1. Check: `DEPLOYMENT.md` → "Troubleshooting"
2. Or: `README_DEPLOYMENT.md` → "Troubleshooting"
3. Check server console for error messages

## 📊 Package Contents Summary

```
✅ Flask web server (complete)
✅ Molecule rendering (RDKit, complete)
✅ Multiple parsers (ChemDoodle, OPSIN, PubChem, Fallback - all included)
✅ Configuration system (fully implemented)
✅ REST API (ready to use)
✅ Web interface (HTML/CSS/JS included)
✅ All dependencies (requirements.txt with exact versions)
✅ Setup scripts (Windows and cross-platform)
✅ Documentation (focused, complete, non-excessive)
```

## ✨ Key Features

✅ **Multiple Input Formats**
- SMILES notation
- IUPAC names (with OPSIN)
- Common compound names
- Cheminformatics databases

✅ **Smart Rendering**
- Auto-scaling labels
- Proper aromatic ring display
- Functional group support
- Stereochemistry handling

✅ **Highly Configurable**
- Label appearance
- Parser selection
- Server settings
- Output formatting

✅ **REST API**
- `/api/render` - Render SMILES
- `/api/render-name` - Render from name
- `/api/name-to-smiles` - Convert name to SMILES
- See README_DEPLOYMENT.md for full API

✅ **Production Ready**
- Error handling
- Configuration validation
- Multiple parser fallback
- CORS support

## 🎓 Quick Reference

| Need | Read | Action |
|------|------|--------|
| Quick start | START_HERE_README.md | Follow 3 steps |
| Full guide | README_DEPLOYMENT.md | Follow detailed instructions |
| Customize appearance | CONFIG_GUIDE.md | Edit app/config.py |
| Change parser | PARSER_CONFIG_GUIDE.md | Edit app/config.py |
| Deploy to server | DEPLOYMENT.md | Follow scenario instructions |
| API docs | README_DEPLOYMENT.md | See API Reference section |
| Troubleshoot | DEPLOYMENT.md | See Troubleshooting section |
| Verify setup | PACKAGE_CHECKLIST.md | Run checklist items |

## 📋 Files at a Glance

### To Run
- ⭐⭐⭐ `run_server.py` - Start here
- ⭐⭐⭐ `requirements.txt` - Install dependencies first

### To Configure
- ⭐⭐⭐ `app/config.py` - All settings in one file
- CONFIG_GUIDE.md - How to customize labels
- PARSER_CONFIG_GUIDE.md - How to configure parsers

### To Understand
- START_HERE_README.md - Quick overview
- README_DEPLOYMENT.md - Complete reference
- DELIVERY_SUMMARY.md - What you got

### To Deploy
- DEPLOYMENT.md - Production scenarios
- setup_opsin.py - Optional OPSIN setup
- PACKAGE_CHECKLIST.md - Verification

## 🚀 Deployment Quick Links

**Local Development:**
```bash
pip install -r requirements.txt
python run_server.py
```

**Production (Gunicorn):**
```bash
pip install gunicorn
gunicorn -w 4 -b 0.0.0.0:5000 app:app
```

**Docker:**
```bash
docker-compose up
```

**See README_DEPLOYMENT.md for more scenarios**

## ✅ Pre-Flight Checklist

Before you start:

- [ ] Python 3.8+ installed
- [ ] Read START_HERE_README.md
- [ ] Run: `pip install -r requirements.txt`
- [ ] Run: `python run_server.py`
- [ ] Open: `http://localhost:5000`
- [ ] Try: Enter "ethanol" in search

If all ✅, you're ready!

## 💡 Pro Tips

1. **The only file you need to edit is `app/config.py`**
   - All settings are there
   - Well-commented
   - Easy to understand

2. **If something breaks:**
   - Check server console output
   - Read DEPLOYMENT.md troubleshooting
   - Reset config to defaults

3. **To use different port:**
   - Edit `app/config.py`
   - Change `PORT = 5000` to `PORT = 8080`

4. **To improve IUPAC parsing:**
   - Run `python setup_opsin.py`
   - Make sure Java is installed
   - Restart server

5. **To integrate with other apps:**
   - Use REST API
   - See README_DEPLOYMENT.md API section
   - All endpoints return JSON

## 📞 Common Issue Resolution

| Issue | Solution |
|-------|----------|
| "ModuleNotFoundError: No module named 'flask'" | Run: `pip install -r requirements.txt` |
| Port 5000 already in use | Edit `app/config.py`: `PORT = 8080` |
| OPSIN not working | Ensure Java 8+ installed, run: `python setup_opsin.py` |
| Labels too small/large | See CONFIG_GUIDE.md or edit `CARBON_LABEL_FONT_SIZE` in app/config.py |
| Parser not recognizing names | Edit `NOMENCLATURE_PARSER` in app/config.py or read PARSER_CONFIG_GUIDE.md |

**For more issues:** See DEPLOYMENT.md Troubleshooting

---

## 🎉 You're All Set!

This package contains **everything** you need:

✅ Complete application code  
✅ All dependencies  
✅ Configuration system  
✅ Setup automation  
✅ Comprehensive documentation  
✅ Ready for production  

**Next Step:**
1. Read `START_HERE_README.md`
2. Run `python run_server.py`
3. Enjoy! 🚀

---

**Package Status:** ✅ COMPLETE AND READY FOR DEPLOYMENT
