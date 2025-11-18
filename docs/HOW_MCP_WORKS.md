# How the MCP Server Works - Simple Explanation

## What is MCP?

MCP = Model Context Protocol

It's just a bridge between me (Claude) and your tools.

**Without MCP:**
- I tell you to run a command
- You run it manually
- You tell me results

**With MCP:**
- I call the MCP server directly
- MCP runs the command
- I get results immediately
- I can make decisions based on results

---

## How It Works (Step by Step)

### You Tell Me a Task
```
User: "Make mol2chemfig SVGs bigger"
```

### I Use MCP to Check System
```python
# I call: mcp_server.py check_status
Result:
{
  "flask_running": true,
  "tests_available": true
}
```

### I Make Code Changes
```
I edit: chem-extension/content.js
Change: width 250px → 350px
```

### I Use MCP to Run Tests
```python
# I call: mcp_server.py run_tests
Result:
[PASS] 13/13 tests pass ✅
```

### I Report Results
```
Done! All tests pass.
Changed: chem-extension/content.js lines 900-910
Tests: 13/13 PASS
```

---

## API Keys

**None needed!**

MCP server runs locally on your computer.
- No external services
- No API keys
- No authentication
- Just local Python scripts

---

## How It Runs Tests

### Test Script: `tests/test_api.py`

```python
# This script:
1. Connects to http://localhost:5000 (Flask)
2. Makes test requests to API endpoints
3. Checks if options are included in URLs
4. Verifies cache filenames have option names
5. Confirms performance (2nd request faster)
6. Reports: PASS or FAIL for each
```

**Example test:**
```
Test: Aromatic circles option
→ Send: http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true
← Expect: cache_url contains "aromatic_circles"
✅ PASS
```

All 13 tests run in ~10 seconds.

---

## How I Divide Work

### Single Agent (Me)

I work alone because:
- ✅ One context window (what I remember)
- ✅ Consistent decisions
- ✅ No coordination overhead
- ✅ Faster execution

**Workflow:**
```
Task 1 (20 min)
  └─ Make change
  └─ Run tests (10 sec)
  └─ Report done

Task 2 (20 min)
  └─ Make change
  └─ Run tests (10 sec)
  └─ Report done

(repeat)
```

### Why Not Multiple Agents?

**Multiple agents would:**
- Require coordination overhead
- Need to share state (complex)
- Have communication delays
- Be harder to debug if something breaks
- Waste context switching

**One agent is better for:**
- Sequential tasks (one after another)
- Interconnected changes (depend on each other)
- Quality consistency
- Debugging

---

## Where MCP Connects

### Connection Flow

```
You (tell me task)
    ↓
Me (Claude in VS Code)
    ↓
MCP Server (mcp_server.py)
    ↓
Local Tools:
  - python test_api.py (tests)
  - Flask (development server)
  - File system (edit code)
    ↓
Results come back to me
    ↓
I report to you
```

### No Network Communication

Everything runs **locally**:
- MCP server on your computer
- Flask server on your computer (http://localhost:5000)
- Tests on your computer
- Files on your computer

---

## When I'll Use MCP

### Scenario 1: Start of Day
```
I call: mcp_server.py check_status
→ "Is Flask running? Are tests available?"
→ If Flask stopped, I start it
```

### Scenario 2: Make Changes
```
I edit: chem-extension/content.js
I call: mcp_server.py run_tests
→ "Do tests pass?"
→ If fail: debug and fix
→ If pass: report done
```

### Scenario 3: Verify System
```
I call: mcp_server.py run_tests
→ Ensure 13/13 pass
→ Before reporting feature complete
```

---

## When You'll Use MCP

### Manual Commands (Optional)
```bash
# If you want to check things yourself:
python mcp_server.py check_status
python mcp_server.py run_tests
python mcp_server.py start_flask
```

But you don't need to. I'll handle this.

---

## The Actual Workflow

### Task: "Make mol2chemfig SVGs bigger"

**What Happens:**

```
Step 1: I read requirement (2 min)
Step 2: I find file (content.js)
Step 3: I edit code (3 min)
        width: 250px → 350px
Step 4: I check status (MCP)
        mcp_server.py check_status
        → Flask running ✅
Step 5: I run tests (MCP)
        mcp_server.py run_tests
        → [PASS] 13/13 ✅
Step 6: I create visual test (5 min)
Step 7: I report: "Done! Tests pass. SVGs 40% bigger"

Total: ~20 minutes
```

---

## No Lower LLM Needed

**Why?**
- One agent is efficient
- Tests are simple (just run .py)
- Tasks are sequential (not parallel)
- No need for multiple LLMs

**Just me working through tasks one by one.**

---

## How Files Transfer

### Code Flow

```
You have files:
chem-extension/content.js
MoleculeViewer/app/api.py
tests/test_api.py

I can:
✅ Read files (understand code)
✅ Edit files (make changes)
✅ Run scripts (via MCP)
✅ Report results (show you output)

All changes stay in your folder.
```

---

## Summary

### MCP Does:
- Runs test script (`python test_api.py`)
- Starts/stops Flask
- Checks system status
- Reads task list

### I Do:
- Make code decisions
- Edit files
- Interpret test results
- Report to you

### You Do:
- Keep Flask running (background)
- Tell me tasks
- Review results

### No One Needs:
- ✅ API keys
- ✅ External services
- ✅ Multiple agents
- ✅ Network communication

**Everything local, everything simple.**

---

## Ready?

Just tell me:
```
"Work on: Make mol2chemfig SVGs bigger"
```

I'll:
1. Use MCP to verify setup
2. Make changes
3. Run tests via MCP
4. Report when done

No other explanation needed!
