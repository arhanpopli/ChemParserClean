# ChemParser - Complete Setup Guide

## 📋 What's in This Repository

This is a **lightweight source-code repository** optimized for Claude AI code editing. It contains all the code Claude needs to modify, with heavy dependencies excluded (they're installed via npm/pip).

### 🎯 What's Included:

**Backend Services:**
- `mol2chemfig_server.py` - LaTeX chemical rendering (Flask)
- `pubchem_server.py` - PubChem integration (Node.js)
- `launcher-server.js` - Web interface hub

**Chemistry Engines:**
- `MoleculeViewer/` - SVG-based SMILES renderer (Node.js + Python)
- `mol2chemfig/` - Mol2ChemFig source code
- `backend_source_*` - Backend processing modules
- `chemistry/` - Chemistry utilities

**Frontend:**
- `frontend_source/` - Vue.js frontend source
- `chem-extension/` - Chrome extension

**Interfaces:**
- `unified-interface.html` - Web hub for all tools
- `mol2chemfig-full-interface.html` - Complete Mol2ChemFig with 8 ChemFig options
- Test HTML files in `tests/` and `frontend_static/`

**Configuration:**
- `docker-compose.yml` - Docker setup
- `requirements*.txt` - Python dependencies
- `package.json` - Node.js dependencies
- `1-start-all.bat` - Main startup script
- `dev-start-*.bat` - Individual server scripts

**Documentation:**
- `docs/` folder - 118 markdown files
- `README_LIGHTWEIGHT.md` - This repo's philosophy
- `MoleculeViewer/README_LIGHTWEIGHT.md` - MoleculeViewer setup

---

## 🚀 Quick Start (After Cloning)

### 1. Install Dependencies

**Python packages:**
```bash
pip install -r requirements_server.txt
pip install -r requirements_pubchem.txt
```

**Node packages:**
```bash
npm install
cd MoleculeViewer && npm install && cd ..
```

### 2. Optional: Install OPSIN (For 3D SMILES)

```bash
cd MoleculeViewer
python setup_opsin.py
cd ..
```

### 3. Start All Servers

```bash
1-start-all.bat
```

This starts:
- **Mol2ChemFig** (Flask) on port 5001
- **MoleculeViewer** (Node.js) on port 5000
- **PubChem** (Node.js) on port 5002
- Opens `unified-interface.html` in browser

### 4. Access Web Interfaces

After running `1-start-all.bat`:
- **Unified Hub:** http://localhost:5000/unified-interface.html
- **Mol2ChemFig Complete:** http://localhost:5000/mol2chemfig-full-interface.html
- **MoleculeViewer Direct:** http://localhost:5000

---

## 📁 Folder Structure

```
ChemParserClean/
├── backend_source_*        # Backend processing code
├── chemistry/              # Chemistry utilities
├── chem-extension/         # Chrome extension source
├── docs/                   # 118 documentation files
├── frontend_source/        # Vue.js frontend
├── frontend_static/        # Static assets
├── MoleculeViewer/         # Chemistry rendering engine
│   ├── app/                # Flask backend
│   ├── templates/          # HTML templates
│   └── server.js           # Node.js server
├── mol2chemfig/            # Mol2ChemFig source
├── tests/                  # Test files
├── *.py                    # Backend servers
├── *.html                  # Web interfaces
├── *.bat                   # Startup scripts
├── *.txt                   # Configuration
└── docker-compose.yml      # Docker setup
```

---

## 🔧 Individual Server Startup

If you want to run servers separately:

**Mol2ChemFig (Flask on 5001):**
```bash
dev-start-mol2chemfig.bat
```

**MoleculeViewer (Node.js on 5000):**
```bash
dev-start-moleculeviewer.bat
```

**PubChem (Node.js on 5002):**
```bash
dev-start-pubchem.bat
```

---

## 🛑 Stop All Servers

```bash
util-stop-all.bat
```

---

## ✅ Server Status

Check if all servers are running:
```bash
util-status.bat
```

---

## 📝 Environment Variables

Create a `.env` file in the root with:

```
# Server Configuration
FLASK_ENV=development
FLASK_DEBUG=1

# Ports (default if not set)
MOL2CHEMFIG_PORT=5001
MOLECULEVIEWER_PORT=5000
PUBCHEM_PORT=5002

# Optional: API Keys
PUBCHEM_API_KEY=your_key_here

# Optional: Docker settings
MOL2CHEMFIG_BACKEND=http://localhost:8000
```

See `.env.example` for more options.

---

## 🐳 Docker Setup (Optional)

If Docker is installed:

```bash
docker-compose up -d
```

This starts the mol2chemfig backend on port 8000 (optional for enhanced features).

---

## 🧪 Testing

Run test files to verify setup:

```bash
# Python tests
python test_runner.py

# Node.js tests (if available)
npm test
```

---

## 🔄 Git Workflow with Claude

### After Claude Makes Changes:

1. **Pull changes:**
   ```bash
   git pull origin main
   ```

2. **See what changed:**
   ```bash
   git log --oneline -5
   ```

3. **Reinstall if dependencies changed:**
   ```bash
   npm install
   pip install -r requirements_server.txt
   ```

4. **Test locally:**
   ```bash
   1-start-all.bat
   ```

---

## ❓ Troubleshooting

### Port Already in Use:
```bash
util-stop-all.bat
```

### Module Not Found:
```bash
pip install -r requirements_server.txt
npm install
```

### Docker Issues:
```bash
docker ps
docker-compose down
docker-compose up -d
```

---

## 📚 Full Documentation

See `docs/` folder for:
- `00_START_HERE.md` - Quick start
- `BATCH_FILES_GUIDE.md` - Batch file explanations
- `INTERFACES_GUIDE.md` - Interface documentation
- `FIXES_SUMMARY.md` - Recent fixes
- Architecture and deployment docs

---

## 🎯 For Claude AI Editing:

Claude has full access to:
- ✅ All source code (.py, .js, .vue, .html, .css)
- ✅ Configuration files
- ✅ Documentation
- ✅ Test files

Claude should NOT modify:
- ❌ node_modules/ (install via npm)
- ❌ venv/ (install via pip)
- ❌ cache/ folders
- ❌ .git/ (handled via GitHub)

---

## 📞 Key Files for Claude:

| File | Purpose |
|------|---------|
| `mol2chemfig_server.py` | Main Mol2ChemFig Flask server |
| `pubchem_server.py` | PubChem integration |
| `MoleculeViewer/server.js` | MoleculeViewer Node server |
| `unified-interface.html` | Web interface hub |
| `mol2chemfig-full-interface.html` | Complete Mol2ChemFig UI |
| `1-start-all.bat` | Master startup script |
| `requirements_*.txt` | Python dependencies |
| `package.json` | Node.js dependencies |

---

**Happy coding!** 🚀
