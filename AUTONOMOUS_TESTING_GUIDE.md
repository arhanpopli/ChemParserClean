# 🤖 Autonomous Testing & Development Guide

## How I Can Work Through Your Todolist Autonomously

### What I Need From You (Setup Requirements)

1. **Keep Flask server running:**
   ```powershell
   cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"
   python run_server.py
   # Keep this window open
   ```

2. **Keep Docker mol2chemfig running** (if needed):
   ```powershell
   docker run -it -p 8000:8000 mol2chemfig:latest
   # Keep this window open
   ```

3. **That's it!** You don't need to do anything else.

---

## How I'll Work

### For Each Task in Todolist:

1. **Read the requirement** from `Todolist.md`
2. **Make code changes** in the relevant files
3. **Run automated tests** using:
   - Python test script: `python tests/test_api.py`
   - Visual test suite: Open `tests/test_visual.html` in browser
   - Terminal commands: Verify file changes, syntax errors, etc.
4. **Report results:**
   - ✅ PASS: Show evidence (test output, screenshot)
   - ❌ FAIL: Debug, fix code, retry
5. **Move to next task**

---

## Available Testing Tools I Have

### 1. **API Testing** (`tests/test_api.py`)
```bash
python tests/test_api.py
```
Tests:
- Server connection
- SMILES endpoint with options
- Cache filename encoding
- Multiple options together
- Nomenclature endpoint

**Example output:**
```
============================================================
TESTING: SMILES Endpoint with Options
============================================================
✅ PASS: Basic SMILES request
✅ PASS: Response has SVG
✅ PASS: Aromatic circles option accepted
✅ PASS: Cache URL contains 'aromatic_circles'
```

### 2. **Visual Testing** (`tests/test_visual.html`)
- Open in browser at `file:///C:/Users/Kapil/Personal/PROJECTS/Mol2chemfig/tests/test_visual.html`
- Runs 5 test cases
- Shows images and cache URLs
- Visual pass/fail indicators
- Summary statistics

### 3. **Browser Testing** (Simple Browser)
- I can open `http://localhost:5000` directly
- Take screenshots
- Verify UI changes
- Test interactive features

### 4. **Extension Testing** (Terminal + File verification)
- I can verify popup.html UI changes
- Check content.js logic
- Validate chrome.storage interactions
- Run Python code to simulate extension calls

---

## Workflow Example: Adding Image Size Controls

**Task:** "Add up/down arrows to control image size in extension"

**My workflow:**
```
1. Read requirement from Todolist.md
2. Open chem-extension/popup.html and read current structure
3. Add HTML for size controls:
   <button onclick="decreaseSize()">↓</button>
   <span id="sizeDisplay">100%</span>
   <button onclick="increaseSize()">↑</button>
4. Add JavaScript in popup.js to handle clicks and save to chrome.storage
5. Add CSS styling
6. Create test HTML file to verify UI renders correctly
7. Open test page with simple browser and verify visually
8. If visual looks good: run terminal test to verify storage is working
9. If all pass: Mark task as COMPLETE ✅
10. Move to next task
```

---

## What Tests You Should Expect Me to Run

### For Backend Changes:
```bash
python tests/test_api.py              # Automated API test
curl http://localhost:5000/img/smiles?smiles=c1ccccc1&show_carbons=true  # Manual verification
```

### For Extension Changes:
```bash
# Verify syntax
python -m py_compile chem-extension/content.js  # Won't work (JS), so I'll read file visually

# Check for specific patterns
grep -n "aromatic_circles" chem-extension/content.js
grep -n "chrome.storage.sync" chem-extension/popup.js
```

### For UI Changes:
```
1. Create test HTML file with new UI
2. Open in simple browser
3. Take screenshot
4. Compare with requirement
```

---

## Priority Order (My Recommendation)

Start with **easiest** tasks first to build momentum:

1. ✅ **Already Done:** Options flowing through to backend ✓
2. 📌 **NEXT (Easy):** Make mol2chemfig SVGs bigger by default
3. 📌 **Next (Medium):** Add image size controls (↑/↓ arrows)
4. 📌 **Next (Medium):** Add dev mode for saving size config
5. 📌 **Next (Hard):** Add separate cache folders for each renderer

---

## How to Give Me Work

### Option A: Update the Todolist.md file:
Edit the file with clear, numbered tasks:
```markdown
# PRIORITY 1 (Do This First)
- [ ] Make mol2chemfig SVGs bigger by default
  - File: chem-extension/content.js (line ~900)
  - Find: container width/height defaults
  - Change: Increase from 250px to 350px
  - Test: Open test_visual.html and compare sizes visually

# PRIORITY 2 (Do After 1 Passes)
- [ ] Add size control buttons
  - Location: bottom-left of each image
  - Buttons: ↑ (increase) and ↓ (decrease)
  - Save to: chrome.storage.sync
  - Test: Click buttons, reload page, verify size persists
```

### Option B: Just Tell Me:
"Work on feature X" and I'll:
1. Read Todolist.md for context
2. Break it into subtasks
3. Work through them autonomously
4. Report progress with evidence

---

## What I'll Produce for Each Task

✅ **When Task Passes:**
```
TASK: Make mol2chemfig SVGs bigger
STATUS: ✅ COMPLETE

CHANGES MADE:
- Modified: chem-extension/content.js (lines 900-950)
- Changed: SVG container width from 250px → 350px
- Changed: height from 200px → 300px

TEST RESULTS:
✅ test_visual.html: SVG displays 40% larger
✅ Console logs show new dimensions
✅ Persists across page reloads

EVIDENCE:
[Screenshot showing larger molecule diagram]
```

❌ **When Task Fails:**
```
TASK: Make mol2chemfig SVGs bigger
STATUS: ⚠️ IN PROGRESS

PROBLEM FOUND:
- SVG container has max-width CSS constraint
- Size change not reflected in rendered output

DEBUG INFO:
- Found max-width: 250px in popup.css line 350
- Solution: Need to update both JS and CSS

NEXT STEP:
- Updating CSS constraints
- Will re-test after
```

---

## Tools I Can Use

| Tool | What It Does | Example |
|------|-------------|---------|
| `run_in_terminal` | Run Python, PowerShell, git commands | `python tests/test_api.py` |
| `mcp_pylance_mcp_s_pylanceRunCodeSnippet` | Run Python code directly | Test option encoding logic |
| `open_simple_browser` | Open URL in VS Code browser | View http://localhost:5000 |
| `create_file` | Create new test files | Test HTML, test scripts |
| `replace_string_in_file` | Edit code files | Add size control buttons |
| `read_file` | Read code to understand | Study current implementation |
| `grep_search` | Find patterns in code | Find all chrome.storage calls |
| `manage_todo_list` | Track task progress | Mark tasks complete |

---

## What I CAN'T Do (Limitations)

❌ Test Chrome extension UI interactively (can't click buttons in extension)
   → *Workaround:* Create test HTML that simulates extension behavior

❌ View real-time browser rendering of extension popup
   → *Workaround:* Create standalone HTML version of UI to verify

❌ Install/run full Chrome extension automatically
   → *Workaround:* Verify code is correct, then you load extension manually for final test

---

## How to Start

### Right Now:

1. **Create test directory** (I'll do this):
   ```bash
   mkdir -p C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\tests
   ```

2. **Copy test files** (Already done ✅):
   - `tests/test_api.py` - Python API test suite
   - `tests/test_visual.html` - Visual test UI

3. **Update Todolist.md** with clear, actionable tasks

4. **Tell me:** "Start working on [FEATURE]"

5. **I'll:**
   - Read the requirement
   - Make changes
   - Run tests autonomously
   - Report progress

---

## Ready to Start?

**Yes!** I'm ready to:
- Read your `Todolist.md`
- Pick highest priority task
- Work through it autonomously
- Use the test suite to verify
- Report results with evidence
- Move to next task

**Just keep the Flask server running** and I'll handle everything else!

---

## Questions?

If you want me to start on a specific task, just say:
- "Fix mol2chemfig SVG sizing"
- "Add image size controls"
- "Add dark mode support for mol2chemfig"
- "Add separate cache folders"
- etc.

I'll take it from there! 🚀
