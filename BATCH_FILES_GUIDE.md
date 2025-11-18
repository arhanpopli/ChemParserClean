# ChemParser Batch Files - Quick Reference

## 🎯 Overview
All 18 redundant batch files have been consolidated into **7 organized scripts**:
- **1 main startup** (starts everything)
- **3 dev servers** (for individual testing)
- **2 utilities** (stop all, check status)
- **1 setup file** (MCP configuration)

---

## 📋 Your Batch Files

### 🚀 MAIN STARTUP (START HERE)
**`1-start-all.bat`** - Start the entire system
```
What it does:
  ✓ Kills any existing processes (avoids port conflicts)
  ✓ Starts Mol2ChemFig (Flask, Port 5001)
  ✓ Starts MoleculeViewer (Node.js, Port 5000)
  ✓ Starts PubChem 3D (Node.js, Port 5002)
  ✓ Opens unified interface in your browser

When to use: You want to run everything
Result: 3 terminal windows + browser opens
```

---

### 🔧 DEVELOPMENT SERVERS (Individual testing)

**`dev-start-mol2chemfig.bat`** - Mol2ChemFig only
```
Port: 5001
Type: Flask (Python)
Use when: Testing Mol2ChemFig rendering without other servers
Access: http://localhost:5001
```

**`dev-start-moleculeviewer.bat`** - MoleculeViewer only
```
Port: 5000
Type: Node.js
Use when: Testing SVG generation and nomenclature conversion
Access: http://localhost:5000
```

**`dev-start-pubchem.bat`** - PubChem 3D only
```
Port: 5002
Type: Node.js
Use when: Testing 3D molecule viewer and image fetching
Access: http://localhost:5002
```

---

### ⚙️ UTILITIES

**`util-stop-all.bat`** - Stop all running servers
```
What it does:
  ✓ Kills all Node.js processes
  ✓ Kills all Python processes
  
Safe to run: Yes, won't error if nothing is running
Use when: You want to shut everything down
```

**`util-status.bat`** - Check server status
```
What it does:
  ✓ Checks if each port (5000, 5001, 5002) is listening
  ✓ Shows which servers are ONLINE or OFFLINE

Use when: You want to verify servers are running
```

**`setup_mcp.bat`** - MCP server setup
```
Special use: Model Context Protocol configuration
(Not needed for regular operation)
```

---

## 🎬 Common Workflows

### "I want to start everything"
```
Double-click: 1-start-all.bat
```

### "I want to test just one server"
```
For Mol2ChemFig:   dev-start-mol2chemfig.bat
For MoleculeViewer: dev-start-moleculeviewer.bat
For PubChem 3D:     dev-start-pubchem.bat
```

### "Servers are acting weird"
```
1. Run: util-stop-all.bat
2. Wait 2 seconds
3. Run: 1-start-all.bat
```

### "I want to check if servers are running"
```
Double-click: util-status.bat
```

### "I want to stop everything"
```
Double-click: util-stop-all.bat
```

---

## 📊 File Organization Summary

| Category | Files | Purpose |
|----------|-------|---------|
| **Main** | 1-start-all.bat | Start entire system |
| **Dev** | dev-start-*.bat (3) | Individual server testing |
| **Utils** | util-stop-all.bat, util-status.bat | Control & monitoring |
| **Setup** | setup_mcp.bat | MCP configuration |

**Total: 7 files** (down from 18 redundant files)

---

## 🔗 Quick Links

- **Unified Interface**: http://localhost:5000/unified-interface.html
- **Mol2ChemFig**: http://localhost:5001
- **MoleculeViewer**: http://localhost:5000
- **PubChem**: http://localhost:5002

---

## ❌ Deleted Files (No longer needed)

The following 17 redundant files were removed:
```
quick_start.bat
start_all.bat
START_ALL_SERVERS.bat
SIMPLE_START_ALL.bat
start_backend.bat
START_BACKEND_SIMPLE.bat
start_docker.bat
start-docker.bat
START_HERE.bat
START_LAUNCHER.bat
start_mol2chemfig.bat
start_moleculeviewer.bat
start_pubchem.bat
START_PUBCHEM_3D.bat
status.bat
stop_all.bat
test.bat
```

These were causing confusion due to similar naming. All functionality is now in the 7 organized scripts above.

---

## 💡 Pro Tips

1. **Always use `1-start-all.bat` first** - It's the main entry point
2. **Keep browser tab open** - Servers will stay running even if you close the batch windows
3. **Ports persist** - If a script won't start, run `util-stop-all.bat` first
4. **Check status anytime** - Run `util-status.bat` to verify servers are responding

