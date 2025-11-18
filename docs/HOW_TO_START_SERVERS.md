# How to Start ChemParser Servers - TWO OPTIONS

I understand the control panel wasn't working for you. Here are **TWO OPTIONS** - one simple, one advanced:

---

## ⭐ OPTION 1: SIMPLE (RECOMMENDED FOR YOU)

### Just double-click: `SIMPLE_START_ALL.bat`

**What it does:**
- Directly starts all 3 servers in separate windows
- No control panel needed
- No Node.js backend
- Simple and reliable
- Each server runs in its own window

**What you'll see:**
- 3 new command windows pop up:
  1. "MoleculeViewer Server - Port 5000"
  2. "mol2chemfig Docker - Port 8000"
  3. "PubChem Server - Port 5002"
- Health checks confirm they're running
- Done!

**To stop:**
- Close each window, OR
- Run `stop_all.bat`

---

## 🎮 OPTION 2: CONTROL PANEL (If you want to try again)

### Double-click: `START_LAUNCHER.bat`

**What it does:**
- Starts a Node.js control server (launcher-server.js)
- Opens web control panel in browser
- Lets you start/stop servers with buttons

**Why it might not have worked before:**
- Port 3000 was already in use
- Node.js modules (express, cors) not installed
- Browser opened before server was ready
- Working directory was wrong

**I've fixed it now with:**
- ✅ Automatic server health checks (retries 10 times)
- ✅ Correct working directory set
- ✅ Clear error messages if it fails
- ✅ Only opens browser when server is confirmed ready

**If it still doesn't work:**
1. Check if port 3000 is in use: `netstat -ano | findstr :3000`
2. Look at the "ChemParser Launcher Control" window for errors
3. Make sure you have express installed: `npm install express cors`
4. Try OPTION 1 instead (it's simpler!)

---

## Comparison

| Feature | SIMPLE_START_ALL.bat | START_LAUNCHER.bat |
|---------|----------------------|---------------------|
| **Complexity** | Very Simple | Medium |
| **Dependencies** | Node.js, Python, Docker | + express, cors |
| **Control** | Manual (close windows) | Web UI buttons |
| **Reliability** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Speed** | Instant | 3-5 seconds |
| **Best for** | Quick start/stop | Frequent testing |

---

## Which Should You Use?

**Use SIMPLE_START_ALL.bat if:**
- ✅ You just want servers running
- ✅ You don't mind closing windows manually
- ✅ Control panel keeps failing
- ✅ You want it to "just work"

**Use START_LAUNCHER.bat if:**
- ✅ You like web UIs
- ✅ You start/stop servers frequently
- ✅ You want status monitoring
- ✅ The control panel actually works now

---

## What Each Server Does

### 1. MoleculeViewer (Port 5000)
- Uses RDKit to generate SVG images
- Converts SMILES and nomenclature to structures
- Fast and reliable

**Test it:**
```
http://localhost:5000/img/smiles?smiles=CCO
```

### 2. mol2chemfig (Port 8000)
- Docker-based ChemFig generator
- High-quality chemistry structures
- Supports advanced options (aromatic circles, etc.)

**Test it:**
```
http://localhost:8000/
```

### 3. PubChem (Port 5002)
- Fetches structures from PubChem database
- 3D models and 2D images
- Fallback for unknown compounds

**Test it:**
```
http://localhost:5002/health
```

---

## Troubleshooting

### Problem: "Port already in use"

**Solution:**
```bat
REM Find what's using the port
netstat -ano | findstr :5000

REM Kill the process (replace PID with actual number)
taskkill /F /PID <PID>
```

### Problem: "Docker not responding"

**Solution:**
1. Make sure Docker Desktop is running
2. Check Docker status: `docker ps`
3. Restart Docker: `docker-compose down && docker-compose up -d`

### Problem: "Node.js not found"

**Solution:**
1. Install Node.js from: https://nodejs.org/
2. Restart command prompt after install
3. Verify: `node --version`

### Problem: "Python not found"

**Solution:**
1. Install Python from: https://www.python.org/
2. Check "Add Python to PATH" during install
3. Verify: `python --version`

---

## Quick Commands Reference

```bat
REM Start everything (simple)
SIMPLE_START_ALL.bat

REM Start everything (control panel)
START_LAUNCHER.bat

REM Stop everything
stop_all.bat

REM Check status
status.bat

REM Check specific port
curl http://localhost:5000/health
curl http://localhost:8000/
curl http://localhost:5002/health
```

---

## My Recommendation for You

Based on your experience, I recommend:

### ⭐ Use `SIMPLE_START_ALL.bat`

It's straightforward, doesn't rely on the control panel, and just works. You can always try the control panel later if you want.

**Steps:**
1. Double-click `SIMPLE_START_ALL.bat`
2. Wait for 3 windows to open
3. Check the health checks (should all be OK)
4. Start using your extension!
5. When done, run `stop_all.bat`

That's it - no complicated control panel, no Node.js backend issues, just direct server startup.

---

## If You REALLY Want the Control Panel to Work

1. Make sure port 3000 is free: `netstat -ano | findstr :3000`
2. Install dependencies: `npm install express cors`
3. Run `START_LAUNCHER.bat`
4. Watch the output - it now shows detailed status
5. Should see "[SUCCESS] Launcher control server is running!"
6. Browser opens to http://localhost:3000/launcher.html
7. Click buttons - they should work now!

The main fix was making sure the working directory is set correctly and waiting for the server to be ready before opening the browser.

---

**TL;DR**: Use `SIMPLE_START_ALL.bat` - it works without any control panel complexity.
