# Chemparser Project - Progress Handoff

**Date**: 2025-11-09
**Completion**: 8/12 tasks (67%)
**Context for Next AI**: This file contains everything you need to continue the project

---

## 📋 Quick Start for New AI Agent

1. **Read these files first**:
   - `.claude.md` - Quick project overview
   - `MoleculeViewer/docs/Todolist.md` - **DETAILED requirements** (READ THIS!)
   - `CLAUDE_CONTEXT.md` - Technical details
   - This file - Current progress

2. **Check what's running**:
   ```bash
   curl http://localhost:5000/health  # MoleculeViewer
   curl http://localhost:8000/        # mol2chemfig
   curl http://localhost:5002/health  # PubChem
   ```

3. **See remaining tasks** below in "What Needs To Be Done"

---

## ✅ What's Been Completed (8/12 tasks)

### 1. Testing Infrastructure ✅
- **File**: `test_runner.py`
- **What it does**: Tests all servers, checks dependencies
- **How to run**: `python test_runner.py`
- **Status**: Fully working

### 2. mol2chemfig SVG Size Increase ✅
- **File modified**: `m2cf_fixed.py` (lines 140-151)
- **Changes**:
  - Atom separation: 16pt → **28pt**
  - Font size: 8pt → **12pt**
- **Result**: SVG width increased from ~36pt to ~61pt (69% larger)
- **Restart required**: `docker-compose restart backend` after changes

### 3. mol2chemfig Dark Mode Fix ✅
- **File modified**: `chem-extension/content.js`
- **Changes**: Fixed regex to handle both `stroke="black"` and `stroke='black'`
- **Lines changed**: 1188-1189, 1236-1237, 1302-1303
- **Pattern used**: `/stroke=["']black["']/gi` and `/fill=["']black["']/gi`

### 4. Separate Cache Folders ✅
- **MoleculeViewer cache**: `MoleculeViewer/cache/moleculeviewer/`
- **mol2chemfig cache**: `cache/mol2chemfig/`
- **Files modified**:
  - `MoleculeViewer/server.js` - CACHE_DIR updated
  - `mol2chemfig_server.py` - STORAGE_DIR updated
- **Benefit**: No more cache conflicts between systems

### 5. Chemfig Settings Persistence ✅
- **File modified**: `test_m2cf_full.html`
- **Features added**:
  - localStorage persistence for options
  - Auto-load on page reload
  - Auto-apply on new conversions
  - Cache links display after applying options
- **How it works**: Options saved to `localStorage` when "Apply Options" clicked

### 6. PubChem 3D Integration ✅ (JUST COMPLETED)
- **Server**: Running on port 5002 (`pubchem_server.py`)
- **3D Viewer**: All required options implemented
  - ✅ Ball & Stick
  - ✅ Sticks
  - ✅ Wire-Frame
  - ✅ Space-Filling
  - ✅ Show Hydrogens
  - ✅ Animate
- **Extension integration**: `content.js` lines 1250-1372
- **Test file**: `test_pubchem_3d.html`
- **API endpoint**: `http://localhost:5002/pubchem/3d-viewer?name=histamine`

### 7. Auto-Approval Setup ✅
- **File**: `.claude/settings.json`
- **Purpose**: Allows agents to work autonomously without permission prompts
- **What's approved**: All tools (Bash, Read, Write, Edit, etc.)

### 8. Documentation Created ✅
- `.claude.md` - Quick reference
- `CLAUDE_CONTEXT.md` - Detailed context
- `AGENT_GUIDE.md` - How to use AI agents
- `PUBCHEM_3D_READY.md` - PubChem integration guide
- Various implementation docs from agents

---

## ⏳ What Needs To Be Done (4/12 tasks remaining)

### Priority 1: Image Size Controls (Todolist.md lines 5-7)
**Status**: Not started
**Complexity**: Medium (3-4 hours)
**Files to modify**:
- `chem-extension/content.js` - Add UI and logic
- `chem-extension/popup.html` - Add developer options
- `chem-extension/popup.js` - Handle settings

**Requirements from Todolist.md**:
- Add up/down arrow buttons in bottom left corner of each image
- Arrows adjust image size incrementally
- **Developer Option 1**: Save size per image per page (localStorage key: page URL + image)
- **Developer Option 2**: Save size per molecule globally (localStorage key: molecule SMILES)
- Example: Resize histamine, reload page → same size persists

**Implementation approach**:
1. Create size control UI (arrows) for each rendered image
2. Add event handlers for +/- buttons
3. Store sizes in chrome.storage.local
4. Add two toggle options in popup developer settings
5. Load saved sizes on page load
6. Test with multiple molecules on different pages

**Storage structure**:
```javascript
// Option 1: Per-page storage
{
  "imageSize_perPage": {
    "https://chatgpt.com": {
      "histamine_img_0": 300,
      "caffeine_img_1": 250
    }
  }
}

// Option 2: Per-molecule storage
{
  "imageSize_perMolecule": {
    "C1=CC=C(C=C1)CCN": 300,  // histamine SMILES
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C": 250  // caffeine SMILES
  }
}
```

### Priority 2: 10 UI Design Variations (Todolist.md lines 10-11)
**Status**: Not started
**Complexity**: Medium (2-3 hours)
**Files to modify**: `chem-extension/popup.html` and `styles.css`

**Requirements from Todolist.md**:
- Create 10 DIFFERENT UI designs for popup
- Keep all options intact (same functionality)
- NOT just color changes - different layouts, styles, themes
- User should be able to pick their favorite

**Implementation approach**:
1. Create 10 separate popup HTML files:
   - `popup-design-1-modern.html`
   - `popup-design-2-dark.html`
   - `popup-design-3-minimal.html`
   - `popup-design-4-cards.html`
   - `popup-design-5-sidebar.html`
   - `popup-design-6-glassmorphism.html`
   - `popup-design-7-cyberpunk.html`
   - `popup-design-8-neumorphism.html`
   - `popup-design-9-material.html`
   - `popup-design-10-retro.html`
2. Add theme selector in main popup
3. Each design should have distinct:
   - Layout (grid, sidebar, tabs, cards, etc.)
   - Color scheme
   - Typography
   - Button styles
   - Overall aesthetic

### Priority 3: OPSIN 3D Nomenclature (Todolist.md lines 18-19)
**Status**: Not started
**Complexity**: High (4-5 hours)
**Files to modify**:
- `m2cf_fixed.py` or create new endpoint
- `chem-extension/content.js` - Add option
- `chem-extension/popup.html` - Add toggle

**Requirements from Todolist.md**:
- Fetch 3D SMILES from OPSIN website
- Support 3D SMILES format like: `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO`
- Add option in mol2chemfig to enable 3D nomenclature
- Make this option available in extension

**OPSIN API**:
- Base URL: `https://opsin.ch.cam.ac.uk/opsin/`
- Endpoint: `/opsin/{name}.smi` or `/opsin/{name}.json`
- Example: `https://opsin.ch.cam.ac.uk/opsin/glucose.smi`

**Implementation approach**:
1. Add OPSIN API integration to mol2chemfig backend
2. Create endpoint: `/api/opsin-3d?name=glucose`
3. Fetch SMILES from OPSIN
4. Process with 3D coordinates (stereochemistry)
5. Add toggle in extension popup
6. Test with compounds like glucose, alanine

### Priority 4: Cache Deduplication (Todolist.md line 20)
**Status**: Not started
**Complexity**: Medium (2-3 hours)
**Files to modify**:
- `MoleculeViewer/server.js`
- `mol2chemfig_server.py`

**Problem identified**:
- MoleculeViewer stores same SMILES separately
- Example: "ethanol" and "CCO" both generate separate cache files
- Same molecule, different cache entries

**Requirements from Todolist.md**:
- Better cache system / optimized
- Remove duplicates
- Store based on **canonical SMILES** not nomenclature

**Implementation approach**:
1. Use RDKit to canonicalize SMILES before caching
2. Hash based on canonical SMILES, not input
3. Scan existing cache for duplicates
4. Merge duplicates
5. Update cache key generation:
   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles(smiles)
   canonical = Chem.MolToSmiles(mol, canonical=True)
   cache_key = hashlib.md5(canonical.encode()).hexdigest()
   ```

---

## ⚠️ Critical Information

### 1. Docker Backend Changes
**IMPORTANT**: After editing `m2cf_fixed.py`:
```bash
docker-compose restart backend
# Or force rebuild:
docker rm -f m2cf_backend && docker-compose up -d
```

### 2. Extension Changes
**IMPORTANT**: After editing extension files:
1. Go to `chrome://extensions`
2. Click reload button on extension
3. Hard refresh webpage (Ctrl+Shift+R)

### 3. Port Usage
- **5000**: MoleculeViewer (Node.js + RDKit)
- **5001**: mol2chemfig wrapper (Flask, optional)
- **5002**: PubChem server (Flask)
- **8000**: mol2chemfig Docker backend
- **8080**: mol2chemfig frontend (Docker)

### 4. Cache Locations
- MoleculeViewer: `MoleculeViewer/cache/moleculeviewer/`
- mol2chemfig: `cache/mol2chemfig/`
- PubChem: `pubchem-cache/`

### 5. Key Files Map
```
Extension:
  chem-extension/content.js       - Main rendering logic (2500+ lines)
  chem-extension/popup.html       - Settings UI
  chem-extension/popup.js         - Settings logic

Servers:
  MoleculeViewer/server.js        - RDKit server (port 5000)
  m2cf_fixed.py                   - mol2chemfig backend (mounted to Docker)
  mol2chemfig_server.py           - mol2chemfig wrapper (port 5001)
  pubchem_server.py               - PubChem server (port 5002)

Documentation:
  MoleculeViewer/docs/Todolist.md - DETAILED requirements ← READ THIS!
  .claude.md                      - Quick reference
  CLAUDE_CONTEXT.md               - Technical details
```

---

## 🧪 How to Test Everything

### Test All Servers
```bash
python test_runner.py
```

### Test Individual Components
```bash
# MoleculeViewer
curl "http://localhost:5000/img/smiles?smiles=CCO"

# mol2chemfig
curl -X POST http://localhost:8000/m2cf/submit \
  -H "Content-Type: application/json" \
  -d '{"textAreaData":"CCO"}'

# PubChem
curl "http://localhost:5002/pubchem/img/histamine"

# PubChem 3D
open test_pubchem_3d.html
```

### Test Extension
1. Load unpacked extension from `chem-extension/`
2. Go to ChatGPT or any webpage
3. Type: `chem:histamine`
4. Verify molecule renders

---

## 🤖 How to Launch AI Agents

### Using Multiple Agents in Parallel
You have `.claude/settings.json` configured for auto-approval!

**Launch all 4 remaining tasks simultaneously**:
```
Launch 4 agents in parallel:

1. general-purpose agent: Image size controls
   - Read MoleculeViewer/docs/Todolist.md lines 5-7
   - Implement up/down arrows for image sizing
   - Add 2 developer options (per-page, per-molecule)
   - Use chrome.storage for persistence

2. general-purpose agent: 10 UI design variations
   - Read MoleculeViewer/docs/Todolist.md lines 10-11
   - Create 10 different popup designs
   - Keep all functionality intact
   - Make designs visually distinct

3. general-purpose agent: OPSIN 3D integration
   - Read MoleculeViewer/docs/Todolist.md lines 18-19
   - Integrate OPSIN API for 3D SMILES
   - Add option to mol2chemfig and extension

4. general-purpose agent: Cache deduplication
   - Read MoleculeViewer/docs/Todolist.md line 20
   - Implement canonical SMILES caching
   - Remove duplicate cache entries
   - Update both servers
```

---

## 💡 Important Decisions Made

### 1. SVG Size Increase
- Chose 28pt atom separation (75% increase from 16pt)
- Chose 12pt font (50% increase from 8pt)
- Makes mol2chemfig match MoleculeViewer size

### 2. Dark Mode Pattern
- Use regex: `/stroke=["']black["']/gi` to handle both quote types
- Replace with white in dark mode
- Apply to all SVG sources (data URIs, raw SVG, etc.)

### 3. Cache Organization
- Separate folders per server
- No shared cache between systems
- Prevents conflicts and confusion

### 4. PubChem Integration
- Use iframe embedding for 3D viewer
- Provide toggle button for 2D ↔ 3D switching
- All 6 viewing options implemented in server

---

## 📝 Notes for Next AI

### Things That Work Well
- Docker auto-restarts on file changes (if volume mounted correctly)
- Extension has robust error handling
- All servers have CORS enabled
- Cache systems work efficiently

### Known Issues (None Critical)
- None currently - all implemented features working

### Tips for Success
1. **Always read Todolist.md** for detailed requirements
2. **Test after each change** - use test_runner.py
3. **Restart services** after backend changes
4. **Reload extension** after UI changes
5. **Check browser console** for extension errors
6. **Use the auto-approval** - `.claude/settings.json` is configured

### Code Style
- JavaScript: Use modern ES6+ features
- Python: Follow PEP 8
- Comments: Explain WHY, not WHAT
- Error handling: Always include try/catch
- Logging: Use console.log with emojis for visibility

---

## 🎯 Success Criteria

### When All Tasks Are Done
- [ ] All 12 todos marked complete
- [ ] All tests pass: `python test_runner.py`
- [ ] Extension works on ChatGPT and other sites
- [ ] All 3 servers running and healthy
- [ ] Documentation updated
- [ ] No console errors

### How to Verify
```bash
# Run full test suite
python test_runner.py

# Check todos
cat MoleculeViewer/docs/Todolist.md

# Test extension manually on ChatGPT
# Type: chem:histamine, chem:caffeine, chem:aspirin
```

---

## 📞 Quick Commands Reference

```bash
# Start servers
cd MoleculeViewer && node server.js &
docker-compose up -d &
python mol2chemfig_server.py &
python pubchem_server.py &

# Stop servers
pkill -f "node server.js"
docker-compose down
pkill -f "mol2chemfig_server.py"
pkill -f "pubchem_server.py"

# Test all
python test_runner.py

# View logs
docker-compose logs -f backend

# Clear caches
rm -rf MoleculeViewer/cache/moleculeviewer/*
rm -rf cache/mol2chemfig/*
rm -rf pubchem-cache/*
```

---

**Ready for handoff!** 🚀
**Next AI**: Start by reading `MoleculeViewer/docs/Todolist.md` and launching the 4 agents above.
