# ChemParser Control Panel - Complete Guide

## The Problem You Had

When you opened `launcher.html` directly in your browser, it showed popups saying "Run this command" but didn't actually start anything. **That's because HTML/JavaScript in a browser can't execute system commands** for security reasons.

## The Solution: Two-Part System

The control panel is actually **TWO parts working together**:

1. **launcher-server.js** - Node.js backend that CAN execute commands
2. **launcher.html** - Web UI that sends requests to the backend

Think of it like a TV remote:
- The remote (launcher.html) = sends signals
- The TV (launcher-server.js) = actually does the work

## How to Use It (Simple 3 Steps)

### Step 1: Start the Control Server

**Double-click**: `START_LAUNCHER.bat`

This will:
- Start the launcher control server (Node.js backend)
- Automatically open the control panel in your browser
- Keep running in the background

**You'll see a new command window titled "ChemParser Launcher Control" - DON'T CLOSE IT!**

### Step 2: Use the Web Control Panel

The control panel will open automatically at: `http://localhost:3000/launcher.html`

You'll see:
- **Status badges** for each server (red = stopped, green = running)
- **Start/Stop buttons** for each server
- **Start All** and **Stop All** buttons
- **Real-time health checks** (auto-refreshes every 5 seconds)

### Step 3: Click Buttons to Control Servers

Now you can actually click buttons and servers will start/stop!

Example:
1. Click **"Start"** next to MoleculeViewer
2. Watch the status badge turn green
3. See confirmation in the log viewer
4. Server is now running on port 5000!

## What Each Server Does

| Server | Port | Purpose |
|--------|------|---------|
| **MoleculeViewer** | 5000 | RDKit-based SVG generation |
| **mol2chemfig** | 8000 | ChemFig LaTeX structures (Docker) |
| **PubChem** | 5002 | PubChem API integration |

## The Control Panel Features

### ✅ Real Control Panel Features:

1. **Actual Start/Stop** - Buttons execute commands to start/stop servers
2. **Status Monitoring** - Health checks every 5 seconds
3. **Visual Indicators** - Green = running, Red = stopped, Yellow = starting
4. **Log Viewer** - See what's happening in real-time
5. **Quick Links** - Jump to test pages and endpoints
6. **Stop All** - Emergency stop button for all servers

### ❌ What It's NOT:

- NOT just informational popups
- NOT fake buttons
- NOT a static HTML page

## How It Actually Works

When you click "Start MoleculeViewer":

1. Browser sends HTTP POST to: `http://localhost:3000/api/start/moleculeviewer`
2. launcher-server.js receives request
3. Backend executes: `node MoleculeViewer/server.js`
4. Server starts in a new command window
5. Health check confirms it's running
6. Status badge turns green
7. You see success message!

**The magic**: The Node.js backend (launcher-server.js) has permission to execute system commands, while the HTML in browser doesn't.

## Troubleshooting

### Problem: "Failed to fetch" or "Network error"

**Solution**: The launcher control server isn't running.

1. Close the browser
2. Double-click `START_LAUNCHER.bat` again
3. Make sure you see "ChemParser Launcher Control" window
4. Wait 3 seconds for server to start
5. Browser should open automatically

### Problem: Buttons don't do anything

**Solution**: Check the browser console (F12)

- If you see "Failed to fetch" → launcher-server.js not running
- If you see "success: false" → Check the error message
- If you see CORS errors → Use `http://localhost:3000/launcher.html` not `file://`

### Problem: Server shows "stopped" even after clicking start

**Solution**:

1. Check if the service actually started (look for new command windows)
2. Wait 5 seconds for health check
3. Click "Refresh Status" button
4. Check if port is already in use: `netstat -ano | findstr :5000`

### Problem: Docker (mol2chemfig) won't start

**Solution**:

1. Make sure Docker Desktop is running
2. Check if `.env` file exists (launcher creates it automatically)
3. Manually test: `docker-compose up -d`
4. Look at the Docker command window for errors

## Architecture Diagram

```
┌─────────────────────────────────────────────────┐
│  Your Browser (Chrome/Firefox)                  │
│  http://localhost:3000/launcher.html            │
│                                                  │
│  [Start MoleculeViewer Button] ← You click here │
└────────────────┬────────────────────────────────┘
                 │ HTTP POST /api/start/moleculeviewer
                 ↓
┌─────────────────────────────────────────────────┐
│  Launcher Control Server (launcher-server.js)   │
│  Node.js running on port 3000                   │
│                                                  │
│  Receives request → Executes command             │
└────────────────┬────────────────────────────────┘
                 │ spawn('node', ['server.js'])
                 ↓
┌─────────────────────────────────────────────────┐
│  MoleculeViewer Server                          │
│  Starts in new command window on port 5000      │
└─────────────────────────────────────────────────┘
```

## Manual Control (Alternative)

If you prefer not to use the control panel:

### Start Everything Manually:
```bat
REM Terminal 1
cd MoleculeViewer
node server.js

REM Terminal 2
docker-compose up -d

REM Terminal 3
python pubchem_server.py
```

### Stop Everything Manually:
```bat
REM Stop Node
taskkill /F /IM node.exe

REM Stop Docker
docker-compose down

REM Stop Python
taskkill /F /IM python.exe
```

## API Endpoints (For Developers)

The launcher-server.js provides these endpoints:

```
GET  /api/status                  - Get status of all servers
POST /api/start/:server           - Start specific server
POST /api/stop/:server            - Stop specific server
POST /api/start/all               - Start all servers
POST /api/stop/all                - Stop all servers
GET  /health                      - Launcher health check
```

### Example API Usage:

```bash
# Check status
curl http://localhost:3000/api/status

# Start MoleculeViewer
curl -X POST http://localhost:3000/api/start/moleculeviewer

# Stop all servers
curl -X POST http://localhost:3000/api/stop/all
```

## Files Overview

```
ChemParser/
├── START_LAUNCHER.bat         ← START HERE! Double-click this
├── launcher-server.js         ← Backend (Node.js) - controls servers
├── launcher.html              ← Frontend (Web UI) - the control panel
├── start_all.bat             ← Old method (still works)
├── stop_all.bat              ← Old method (still works)
└── status.bat                ← Old method (still works)
```

## Quick Reference

### To use the control panel:
1. **Double-click**: `START_LAUNCHER.bat`
2. **Wait**: 3 seconds
3. **Use**: The web control panel that opens

### To stop everything:
- **Option 1**: Click "Stop All" in control panel
- **Option 2**: Close the "ChemParser Launcher Control" window
- **Option 3**: Run `stop_all.bat`

### To check if launcher is running:
- Open: http://localhost:3000/health
- Should see: `{"status":"ok","service":"ChemParser Launcher Control",...}`

## Summary

**Before (Your Problem)**:
- Opened launcher.html directly
- Got useless popups
- Nothing actually started
- Very frustrating!

**After (The Solution)**:
- Run START_LAUNCHER.bat first
- Opens control panel automatically
- Click buttons → Servers actually start
- Real control panel that works!

---

**The key insight**: You need the launcher-server.js running **FIRST**, then the launcher.html can send it commands. The .bat file does both steps for you automatically.

Now you have a **legitimate control panel that actually works**! 🎉
