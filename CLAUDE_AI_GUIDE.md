# 🤖 Claude AI - Quick Reference Guide

This guide is specifically for Claude AI to understand the codebase structure and what to edit.

## 📝 Key Files Claude Will Edit

### Backend Servers (Python/Flask/Node.js)

| File | Purpose | Edit For |
|------|---------|----------|
| `mol2chemfig_server.py` | Main Mol2ChemFig Flask server | Fix bugs, add features, handle requests |
| `pubchem_server.py` | PubChem integration server | Add PubChem features, fix 3D viewer |
| `MoleculeViewer/server.js` | Node.js chemistry renderer | SVG generation, SMILES conversion |
| `MoleculeViewer/app/*.py` | Flask chemistry backend | Chemistry calculations, SVG rendering |

### Chemistry Processing

| File | Purpose | Edit For |
|------|---------|----------|
| `mol2chemfig/processor.py` | Core Mol2ChemFig logic | Chemical rendering, LaTeX generation |
| `mol2chemfig/molecule.py` | Molecule representation | Structure handling, bond/atom processing |
| `backend_source_mol2chemfig/*` | Mol2ChemFig modules | Option handling, coordinate calculations |
| `chemistry/utils.py` | Chemistry utilities | SMILES parsing, nomenclature conversion |

### Frontend & Interfaces

| File | Purpose | Edit For |
|------|---------|----------|
| `unified-interface.html` | Web hub for all tools | Tab switching, interface layout |
| `mol2chemfig-full-interface.html` | Complete Mol2ChemFig UI | ChemFig options, drawing tools |
| `MoleculeViewer/templates/*.html` | MoleculeViewer pages | UI layout, user input forms |
| `frontend_source/**/*.vue` | Vue.js components | Advanced frontend work |

### Configuration & Scripts

| File | Purpose | Edit For |
|------|---------|----------|
| `1-start-all.bat` | Master startup script | Server startup order, environment setup |
| `docker-compose.yml` | Docker orchestration | Container setup, port mapping |
| `requirements_server.txt` | Python dependencies | Version management, new packages |
| `package.json` | Node.js dependencies | Script management, new packages |

---

## 🎯 Common Tasks for Claude

### Task 1: Fix a Server Bug

1. **Find the error file** - Usually in `mol2chemfig_server.py` or `pubchem_server.py`
2. **Check the endpoint** - Look for `@app.route()` decorator
3. **Read error details** - Check error message in terminal/logs
4. **Edit the code** - Make the fix
5. **Test locally** - Run `1-start-all.bat` and test

### Task 2: Add a New ChemFig Option

1. **Edit** `mol2chemfig/options.py` - Define the new option
2. **Edit** `mol2chemfig/molecule.py` - Implement the logic
3. **Edit** `mol2chemfig-full-interface.html` - Add UI control (line 780-792 has OPTIONS array)
4. **Test** - Open interface and try the new option

### Task 3: Fix SVG Generation

1. **Edit** `MoleculeViewer/app/chemistry.py` - SVG generation logic
2. **Edit** `MoleculeViewer/server.js` - Server endpoint
3. **Edit** `MoleculeViewer/templates/index.html` - If UI change needed
4. **Test** - Try rendering a molecule

### Task 4: Add PubChem Feature

1. **Edit** `pubchem_server.py` - Add new endpoint
2. **Edit** `unified-interface.html` - Add UI for feature
3. **Test** - Check browser console for errors

---

## 🔍 Code Patterns to Know

### Flask Endpoint Pattern (Python)
```python
@app.route('/api/endpoint', methods=['GET', 'POST'])
def endpoint():
    # Get request data
    data = request.json
    # Process
    result = process(data)
    # Return JSON
    return jsonify(result)
```

### HTML/JavaScript Pattern (Interfaces)
```html
<button onclick="callBackend()">Click Me</button>
<script>
async function callBackend() {
    const response = await fetch('http://localhost:5001/api/endpoint', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ data: 'value' })
    });
    const result = await response.json();
    console.log(result);
}
</script>
```

### ChemFig Options Pattern
ChemFig options are defined in `mol2chemfig/options.py`:
```python
OPTIONS = {
    'aromatic': { 'flag': '-o', 'tooltip': 'Show aromatic circles' },
    'fancy_bonds': { 'flag': '-f', 'tooltip': 'Fancy bond rendering' },
    # Add more options here
}
```

---

## 🧪 Testing Commands

```bash
# Start all servers
1-start-all.bat

# Check server status
util-status.bat

# Stop all servers
util-stop-all.bat

# View logs
type server.log

# Test a specific endpoint
curl http://localhost:5001/api/test

# Check port 5001 (Mol2ChemFig)
netstat -ano | findstr :5001
```

---

## 📚 Important Files to READ First

Before editing, Claude should understand:

1. **`docs/ARCHITECTURE.md`** - System overview
2. **`README_LIGHTWEIGHT.md`** - Why this repo structure
3. **`COMPLETE_SETUP_GUIDE.md`** - Full setup instructions
4. **`INTERFACES_GUIDE.md`** - Interface documentation
5. **`mol2chemfig/README.md`** - If exists, Mol2ChemFig docs
6. **`MoleculeViewer/README.md`** - MoleculeViewer docs

---

## 🚨 DO NOT EDIT

These files should NOT be edited by Claude (they'll be overwritten):
- ❌ `node_modules/**` - Use `package.json` instead
- ❌ `venv/**` - Use `requirements.txt` instead
- ❌ `.git/**` - Use git commands
- ❌ Cache folders - Generated at runtime
- ❌ `.env` - Use `.env.example` as template

---

## 🔄 Git Workflow for Claude

When Claude makes changes:

1. **Claude edits code** - E.g., `mol2chemfig_server.py`
2. **Claude tests locally** - Runs `1-start-all.bat`
3. **Claude commits:**
   ```bash
   git add .
   git commit -m "Fix: Describe the fix"
   git push origin main
   ```
4. **User pulls changes:**
   ```bash
   git pull origin main
   ```

---

## 💡 API Endpoints Reference

### Mol2ChemFig Server (Port 5001)
- `POST /api/mol2chemfig` - Convert molecule to Mol2ChemFig
- `GET /api/options` - Get available ChemFig options
- `POST /api/opsin` - Convert name to 3D SMILES

### MoleculeViewer Server (Port 5000)
- `GET /render?smiles=...` - Render SVG from SMILES
- `GET /mol2chemfig-full-interface.html` - Full interface

### PubChem Server (Port 5002)
- `GET /search?name=...` - Search PubChem
- `GET /structure?cid=...` - Get 3D structure

---

## 🎓 Learning Resources

- **Mol2ChemFig Docs:** `docs/MOL2CHEMFIG_API.md`
- **SMILES Format:** `docs/OPSIN_3D_IMPLEMENTATION.md`
- **ChemFig Guide:** `docs/CHEMFIG_OPTIONS_PERSISTENCE.md`
- **PubChem Integration:** `docs/PUBCHEM_3D_GUIDE.md`

---

## ❓ Questions Claude Should Ask

Before diving into code, Claude should consider:

1. **Where does this error come from?** Check `server.log` first
2. **Which server handles this?** Mol2ChemFig, MoleculeViewer, or PubChem?
3. **What's the user asking for?** New feature or bug fix?
4. **Which interface shows this?** unified-interface.html or mol2chemfig-full-interface.html?
5. **Do I need to update tests?** Check `tests/` folder

---

## 🎯 Quick Wins for Claude

Easy tasks to get started:
1. ✅ Update documentation in `docs/`
2. ✅ Add print statements for debugging
3. ✅ Fix typos in error messages
4. ✅ Improve comments in code
5. ✅ Add new HTML buttons/forms
6. ✅ Update requirements.txt with new packages

---

**Made specifically for Claude AI to understand and contribute to ChemParser!** 🤖

Good luck with the code! 🚀
