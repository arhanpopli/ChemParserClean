# ChemParser Master Launcher

A comprehensive control panel for managing all ChemParser servers and services.

## Overview

The ChemParser Master Launcher provides a unified interface to start, stop, and monitor all servers required for the ChemParser chemical structure visualization system.

## System Architecture

ChemParser consists of **4 main servers**:

### 1. MoleculeViewer Server (Port 5000)
- **Type:** Node.js
- **Purpose:** SMILES and nomenclature to SVG conversion
- **Health Check:** http://localhost:5000/health
- **Start Script:** `start_moleculeviewer.bat`

### 2. Mol2ChemFig Docker Backend (Port 8000)
- **Type:** Docker Container
- **Purpose:** Core mol2chemfig processing engine
- **Health Check:** http://localhost:8000
- **Start Script:** `start_docker.bat`
- **Requirements:** Docker Desktop must be running

### 3. Mol2ChemFig Server (Port 5001)
- **Type:** Flask (Python)
- **Purpose:** Flask wrapper for mol2chemfig with caching and OPSIN integration
- **Health Check:** http://localhost:5001/health
- **Start Script:** `start_mol2chemfig.bat`

### 4. PubChem Server (Port 5002)
- **Type:** Flask (Python)
- **Purpose:** PubChem image fetching and 3D model serving
- **Health Check:** http://localhost:5002/health
- **Start Script:** `start_pubchem.bat` (already exists)

## Quick Start

### Method 1: Use the Master Launcher (Recommended)

1. **Start all servers:**
   ```batch
   start_all.bat
   ```

2. **Open the launcher:**
   - The launcher will open automatically, or
   - Open `launcher.html` in your browser

3. **Monitor status:**
   - The launcher auto-refreshes every 5 seconds
   - Green badges = Running
   - Red badges = Stopped

### Method 2: Manual Server Control

#### Start Individual Servers:
```batch
start_moleculeviewer.bat   # MoleculeViewer (Node.js)
start_docker.bat           # Mol2ChemFig Docker Backend
start_mol2chemfig.bat      # Mol2ChemFig Server (Flask)
start_pubchem.bat          # PubChem Server (Flask)
```

#### Stop All Servers:
```batch
stop_all.bat
```

#### Check Server Status:
```batch
status.bat
```

## File Structure

```
ChemParser/
├── launcher.html              # Master control panel (HTML UI)
├── start_all.bat             # Start all servers
├── stop_all.bat              # Stop all servers
├── status.bat                # Check server status
├── start_moleculeviewer.bat  # Start MoleculeViewer only
├── start_docker.bat          # Start Docker backend only
├── start_mol2chemfig.bat     # Start Mol2ChemFig server only
├── start_pubchem.bat         # Start PubChem server only (existing)
├── LAUNCHER_README.md        # This file
└── docker-compose.yml        # Docker configuration
```

## Features

### Master Launcher UI (`launcher.html`)

- **Server Status Dashboard:** Real-time status indicators for all servers
- **One-Click Controls:** Start/Stop individual or all servers
- **Auto-Refresh:** Status updates every 5 seconds
- **Health Checks:** Automatic health endpoint monitoring
- **Quick Links:** Direct access to test pages and extension
- **Endpoint Reference:** Quick reference for all API endpoints
- **Troubleshooting Tips:** Built-in help for common issues

### Batch Scripts

#### `start_all.bat`
- Checks prerequisites (Python, Node.js, Docker)
- Creates `.env` file if missing
- Starts all servers in separate windows
- Opens the launcher automatically
- Provides status summary

#### `stop_all.bat`
- Gracefully stops all servers
- Stops Docker containers
- Kills Node.js and Python processes
- Cleans up resources

#### `status.bat`
- Checks all server health endpoints
- Lists Docker container status
- Shows running processes
- Provides quick action links

## Prerequisites

### Required Software

1. **Python 3.8+**
   - Download: https://www.python.org/downloads/
   - Required for Flask servers

2. **Node.js 14+**
   - Download: https://nodejs.org/
   - Required for MoleculeViewer server

3. **Docker Desktop** (Optional but recommended)
   - Download: https://www.docker.com/products/docker-desktop
   - Required for Mol2ChemFig Docker backend

### Python Dependencies

Install required packages:
```batch
pip install flask flask-cors requests
pip install rdkit  # For SMILES canonicalization
```

### Node.js Dependencies

```batch
cd MoleculeViewer
npm install
```

## Usage Guide

### Starting the System

1. **First Time Setup:**
   ```batch
   # Install dependencies (one-time)
   pip install -r requirements.txt
   cd MoleculeViewer && npm install && cd ..

   # Start Docker Desktop
   # (if using Mol2ChemFig Docker backend)
   ```

2. **Start All Servers:**
   ```batch
   start_all.bat
   ```

3. **Wait for Initialization:**
   - MoleculeViewer: ~2-3 seconds
   - Docker Backend: ~5-10 seconds (first time)
   - Mol2ChemFig Server: ~2-3 seconds
   - PubChem Server: ~2-3 seconds

4. **Verify Status:**
   - Check the launcher UI (opens automatically)
   - Or run `status.bat`

### Testing the System

After starting servers, test each component:

1. **MoleculeViewer:**
   - Visit: http://localhost:5000/img/smiles?smiles=CCO
   - Should show an SVG of ethanol

2. **Docker Backend:**
   - Visit: http://localhost:8000
   - Should show mol2chemfig backend info

3. **Mol2ChemFig Server:**
   - Open: http://localhost:5001
   - Should show test interface or server info

4. **PubChem Server:**
   - Visit: http://localhost:5002/pubchem/img/histamine
   - Should show PNG of histamine

### Using the Chrome Extension

1. Load the extension:
   - Open Chrome
   - Go to `chrome://extensions/`
   - Enable "Developer mode"
   - Click "Load unpacked"
   - Select the `chem-extension` folder

2. Test inline replacement:
   - Open any webpage or `TEST_EXTENSION.html`
   - Type: `chem:benzene:`
   - Should auto-replace with chemical structure

## API Endpoints Reference

### MoleculeViewer (Port 5000)

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/img/smiles?smiles=CCO` | Render SMILES as SVG |
| GET | `/img/nomenclature?nomenclature=benzene` | Render name as SVG |
| GET | `/health` | Health check |
| GET | `/cache-info` | Cache statistics |

### Mol2ChemFig Docker (Port 8000)

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/m2cf/submit` | Simple SMILES to ChemFig |
| POST | `/m2cf/apply` | SMILES with options |
| POST | `/m2cf/layers` | Layered SVG generation |
| GET | `/m2cf/search` | Name to SMILES lookup |

### Mol2ChemFig Server (Port 5001)

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/health` | Health check |
| POST | `/api/generate` | Generate molecule image |
| GET | `/api/search` | Search by name |
| POST | `/api/generate-3d` | Generate with 3D SMILES (OPSIN) |
| GET | `/api/opsin` | Convert name to 3D SMILES |
| GET | `/images/{hash}.svg` | Serve cached image |

### PubChem Server (Port 5002)

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/health` | Health check |
| GET | `/pubchem/img/{name}` | Direct PNG image |
| GET | `/pubchem/image?name={name}` | Image with metadata |
| GET | `/pubchem/3d-model?name={name}` | 3D SDF model |
| GET | `/pubchem/3d-viewer?name={name}` | Interactive 3D viewer |
| GET | `/pubchem/info?name={name}` | Compound information |

## Troubleshooting

### Port Already in Use

**Problem:** Server fails to start with "port already in use" error

**Solution:**
```batch
# Stop all servers
stop_all.bat

# Check what's using the port (e.g., 5000)
netstat -ano | findstr :5000

# Kill the process using Task Manager or:
taskkill /PID <process_id> /F

# Restart servers
start_all.bat
```

### Docker Not Starting

**Problem:** Mol2ChemFig Docker backend won't start

**Solutions:**
1. Make sure Docker Desktop is running
2. Check Docker daemon:
   ```batch
   docker ps
   ```
3. Restart Docker Desktop
4. Check port 8000 availability:
   ```batch
   netstat -ano | findstr :8000
   ```
5. Check `.env` file exists in project root

### Server Not Responding

**Problem:** Server shows as running but doesn't respond

**Solutions:**
1. Check server logs in terminal window
2. Restart the specific server:
   ```batch
   # Stop
   stop_all.bat

   # Start specific server
   start_moleculeviewer.bat  # or other server
   ```
3. Check firewall settings
4. Verify Python/Node.js installation

### CORS Errors

**Problem:** Browser shows CORS policy errors

**Solutions:**
1. Make sure all servers are running on localhost
2. Check browser console for specific error
3. Restart browser
4. Clear browser cache

### Cache Issues

**Problem:** Old results showing even after changes

**Solutions:**
1. Clear server cache:
   - MoleculeViewer: `DELETE http://localhost:5000/clear-cache`
   - Mol2ChemFig: `POST http://localhost:5001/api/cache/clear`
   - PubChem: `DELETE http://localhost:5002/clear-cache`

2. Clear browser cache:
   - Chrome: Ctrl+Shift+Delete
   - Select "Cached images and files"
   - Click "Clear data"

### Extension Not Working

**Problem:** Chrome extension doesn't replace text

**Solutions:**
1. Make sure all servers are running
2. Reload the extension:
   - Go to `chrome://extensions/`
   - Click reload button on ChemParser extension
3. Check extension popup for settings
4. Try on a test page: `TEST_EXTENSION.html`
5. Check browser console for errors

### Python/Node.js Not Found

**Problem:** Batch scripts can't find Python or Node.js

**Solutions:**
1. Install missing software:
   - Python: https://www.python.org/downloads/
   - Node.js: https://nodejs.org/

2. Add to PATH:
   - Windows: System Properties > Environment Variables
   - Add Python and Node.js installation paths

3. Restart terminal/command prompt

## Advanced Configuration

### Changing Ports

Edit the server files to change default ports:

1. **MoleculeViewer (5000):**
   - Edit `MoleculeViewer/server.js`
   - Change `const PORT = 5000;`

2. **Docker Backend (8000):**
   - Edit `.env` file
   - Change `BACKEND_HOST_PORT=8000`

3. **Mol2ChemFig Server (5001):**
   - Edit `mol2chemfig_server.py`
   - Change `app.run(... port=5001 ...)`

4. **PubChem Server (5002):**
   - Edit `pubchem_server.py`
   - Change `PORT = 5002`

**Note:** After changing ports, update launcher.html server configurations.

### Custom Docker Configuration

Edit `docker-compose.yml` to customize:
- Container ports
- Volume mappings
- Environment variables
- Resource limits

### Logging

Enable detailed logging:

1. **MoleculeViewer:** Already logs to console
2. **Flask Servers:** Set `debug=True` in server files
3. **Docker:** View logs with:
   ```batch
   docker-compose logs -f
   ```

## Development Tips

### Running Servers Individually

For development, you may want to run only specific servers:

```batch
# Only MoleculeViewer
start_moleculeviewer.bat

# Only Mol2ChemFig
start_docker.bat          # Start backend first
start_mol2chemfig.bat     # Then Flask wrapper

# Only PubChem
start_pubchem.bat
```

### Modifying the Launcher

The launcher is a static HTML file (`launcher.html`). To customize:

1. Edit server configurations in JavaScript
2. Modify UI styles in CSS
3. Add new features to the control panel
4. Extend with additional monitoring

**Note:** The launcher can't directly start/stop processes from the browser for security reasons. It provides buttons that guide users to use batch scripts.

## FAQ

**Q: Do I need all servers running?**
A: It depends on your use case:
- Extension only: Needs MoleculeViewer (5000) + Mol2ChemFig (5001)
- 3D features: Also needs PubChem (5002)
- Full features: All servers including Docker (8000)

**Q: Can I run this on Linux/Mac?**
A: The batch files are Windows-specific. For Linux/Mac, you can:
- Start servers manually using the commands in batch files
- Create shell script equivalents
- Use the launcher.html UI for monitoring

**Q: How do I update the system?**
A:
1. Stop all servers: `stop_all.bat`
2. Pull latest code
3. Update dependencies: `pip install -r requirements.txt`
4. Restart servers: `start_all.bat`

**Q: Where are cache files stored?**
A:
- MoleculeViewer: `MoleculeViewer/cache/moleculeviewer/`
- Mol2ChemFig: `cache/mol2chemfig/`
- PubChem: `pubchem-cache/`

**Q: How do I clear all caches?**
A:
```batch
# Delete cache directories
rmdir /s /q MoleculeViewer\cache\moleculeviewer
rmdir /s /q cache\mol2chemfig
rmdir /s /q pubchem-cache

# Or use API endpoints (servers must be running)
curl -X DELETE http://localhost:5000/clear-cache
curl -X POST http://localhost:5001/api/cache/clear
curl -X DELETE http://localhost:5002/clear-cache
```

## Support

For issues or questions:
1. Check this README
2. Review troubleshooting section
3. Check server logs in terminal windows
4. Verify all prerequisites are installed
5. Try `status.bat` to diagnose issues

## License

Part of the ChemParser project.
