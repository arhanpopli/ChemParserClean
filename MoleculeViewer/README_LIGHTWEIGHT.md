# MoleculeViewer - Source Code Only

This folder contains the **source code** for MoleculeViewer. Heavy dependencies are excluded.

## ✅ What's Included:
- `app/` - Flask application source code
- `templates/` - HTML templates
- `server.js` - Node.js server
- `*.py` - Python backend scripts
- `docs/` - Documentation
- All configuration files

## ❌ What's Excluded (in .gitignore):
- `venv/` - 133 MB Python virtual environment
- `opsin-cli.jar` - 13 MB OPSIN library
- `pubchem/` - 6.7 MB PubChem data cache
- `node_modules/` - 4.5 MB NPM packages
- `svg-cache/` - Runtime SVG cache
- `cache/` - Runtime cache

## 🚀 Setup After Clone:

1. **Install Python dependencies:**
   ```bash
   cd MoleculeViewer
   pip install -r requirements.txt
   ```

2. **Install Node dependencies:**
   ```bash
   npm install
   ```

3. **Download OPSIN (optional):**
   ```bash
   python setup_opsin.py
   ```
   Or manually download `opsin-cli.jar` from https://github.com/dan2097/opsin

4. **Start server:**
   ```bash
   node server.js
   # or
   python server.py
   ```

## 💡 Why This Structure?

MoleculeViewer is a separate chemistry rendering engine that:
- Converts SMILES to SVG
- Converts chemical names to structures (via OPSIN)
- Provides an alternative to Mol2ChemFig when Docker isn't available

Claude needs to edit the **source code** (app/, templates/, *.py, *.js) but doesn't need the heavy runtime dependencies (venv, node_modules, jars).
