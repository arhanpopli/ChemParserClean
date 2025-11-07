# MCP Server Setup - 3 Steps

## What is this?
MCP = Model Context Protocol
Lets me (Claude) work on your code autonomously

## Setup (Pick One)

### Option 1: Automatic Setup
```bash
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
setup_mcp.bat
```

### Option 2: Manual Setup
```bash
# Go to project folder
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig

# Test MCP server
python mcp_server.py check_status

# Start Flask
python mcp_server.py start_flask

# Run tests
python mcp_server.py run_tests
```

## What the MCP server does

```
mcp_server.py run_tests      → Runs 13 automated tests
mcp_server.py start_flask    → Starts Flask on port 5000
mcp_server.py stop_flask     → Stops Flask
mcp_server.py check_status   → Shows what's running
mcp_server.py get_task       → Gets next task from Todolist.md
```

## That's it!

Now I can:
✅ Run tests automatically
✅ Start/stop Flask
✅ Check system status
✅ Read your Todolist
✅ Work autonomously

## Just tell me:
"Start working on: Make mol2chemfig SVGs bigger"

And I'll:
1. Use MCP to check Flask is running
2. Make code changes
3. Run tests via MCP
4. Report results

Done! ✅
