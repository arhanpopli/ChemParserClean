# ChemParser - Lightweight Repository

This repository contains **ONLY the source code** that Claude needs to edit. Heavy libraries and runtime data are excluded.

## 📦 What's Included (3-5 MB):
✅ All Python backend files (.py)
✅ All HTML/CSS/JS interfaces
✅ Batch scripts for server management
✅ Configuration files (package.json, docker-compose.yml, requirements.txt)
✅ Documentation (118 markdown files in docs/)
✅ Source code folders:
  - `mol2chemfig/` - Mol2ChemFig source
  - `backend_source_mol2chemfig/` - Backend processing
  - `backend_source_chemistry/` - Chemistry utilities
  - `chemistry/` - Chemistry core
  - `frontend_source/` - Frontend source files
  - `chem-extension/` - Chrome extension

## 🚫 What's Excluded (150 MB - in .gitignore):
❌ `MoleculeViewer/` - 159 MB library (run `npm install` locally)
❌ `ChemDoodleWeb-11.0.0/` - 16 MB library (download separately)
❌ `node_modules/` - 2.5 MB (run `npm install`)
❌ `cache/` - Runtime SVG cache
❌ `pubchem-cache/` - PubChem data cache
❌ `backend_data/` - Runtime backend storage
❌ `mol2chemfig_storage/` - Processing temp files
❌ `__pycache__/` - Python bytecode
❌ Log files (*.log)

## 🚀 Setup After Clone:

1. **Install Node dependencies:**
   ```bash
   npm install
   ```

2. **Install Python dependencies:**
   ```bash
   pip install -r requirements_server.txt
   pip install -r requirements_pubchem.txt
   ```

3. **Download MoleculeViewer** (if needed):
   - Place in `MoleculeViewer/` folder
   - Or clone from original source

4. **Download ChemDoodleWeb** (if needed):
   - Download from https://web.chemdoodle.com/
   - Extract to `ChemDoodleWeb-11.0.0/`

5. **Start servers:**
   ```bash
   1-start-all.bat
   ```

## 💡 Why This Structure?

Claude AI only needs to **edit source code**. The heavy libraries:
- Are static/unchanging
- Can be installed via npm/pip
- Would make the repo 200MB+ (unnecessary)
- Can't be edited by Claude anyway

This keeps the repo **lightweight and fast** for AI code editing! ⚡

## 📚 Full Documentation:

See `docs/` folder for:
- Setup guides
- API references  
- Architecture documentation
- Deployment instructions
