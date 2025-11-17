# 🎮 ChemParser Master Control Panel

## Quick Start

### Step 1: Start the Control Panel
Double-click: **`START_LAUNCHER.bat`**

This starts the control server on port 3000.

### Step 2: Open the Web Interface
The browser will open automatically, or go to:
**http://localhost:3000/launcher.html**

### Step 3: Control Your Servers
Click the buttons to:
- ▶ **Start** - Launch a server
- ⏹ **Stop** - Shut down a server
- 🔗 **Open** - Open server in browser
- ↻ **Refresh** - Update status

## Available Servers

### PubChem Server (Port 5002)
- 3D molecular viewer
- 2D structure images
- Chemical compound data

### MoleculeViewer (Port 5000)
- Main chemistry rendering engine
- Mol2ChemFig converter
- SMILES processing

## How It Works

1. **Launcher Control Server (Port 3000)**
   - Node.js backend that controls other servers
   - Provides API for start/stop/status
   - Must be running for buttons to work

2. **Web Control Panel**
   - Beautiful HTML interface
   - Real-time server status
   - One-click start/stop buttons

3. **Managed Servers**
   - PubChem (Node.js on 5002)
   - MoleculeViewer (Node.js on 5000)
   - More can be added

## Features

✅ **Visual Status Indicators**
- 🟢 Green = Running
- 🔴 Red = Stopped
- 🟡 Yellow = Checking

✅ **One-Click Controls**
- Start individual servers
- Stop individual servers
- Start ALL at once
- Stop ALL at once

✅ **Auto-Refresh**
- Status updates every 5 seconds
- See changes in real-time

✅ **Quick Links**
- Direct links to test pages
- Extension popup
- Health checks

## Troubleshooting

### "Control Server Disconnected"
- Run `START_LAUNCHER.bat`
- Wait for "ChemParser Launcher Control Server" message

### "Port already in use"
- Click Stop button first
- Or close the terminal window running that server

### Buttons don't work
- Check if control server is running (green status at top)
- Refresh the page
- Check browser console (F12) for errors

## Manual Control

If you prefer manual control, you can still use:
- `START_PUBCHEM_3D.bat` - Just PubChem server
- `START_ALL_SERVERS.bat` - All servers in separate windows

## Architecture

```
START_LAUNCHER.bat
    ↓
launcher-server.js (Port 3000)
    ↓
launcher.html (Web UI)
    ↓
    ├─→ Start/Stop PubChem (Port 5002)
    ├─→ Start/Stop MoleculeViewer (Port 5000)
    └─→ Start/Stop Other Servers
```

## Benefits

✨ **No more terminal juggling**
✨ **See all server status at once**
✨ **Start/stop with one click**
✨ **Beautiful modern interface**
✨ **Real-time monitoring**

Enjoy your master control panel! 🚀
