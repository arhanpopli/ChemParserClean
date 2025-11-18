# ChemParser Quick Start Guide

## 🚀 Getting Started in 3 Steps

### Step 1: Install Prerequisites

Make sure you have:
- ✅ **Python 3.8+** - [Download](https://www.python.org/downloads/)
- ✅ **Node.js 14+** - [Download](https://nodejs.org/)
- ⚠️ **Docker Desktop** (Optional) - [Download](https://www.docker.com/products/docker-desktop)

### Step 2: Start All Servers

**Double-click this file:**
```
start_all.bat
```

Wait 10-15 seconds for all servers to initialize.

### Step 3: Open the Launcher

The launcher will open automatically, or open:
```
launcher.html
```

---

## 🎯 What Just Happened?

The launcher started **4 servers**:

| Server | Port | Purpose |
|--------|------|---------|
| 🧪 MoleculeViewer | 5000 | SMILES/Name → SVG |
| ⚗️ Docker Backend | 8000 | Mol2ChemFig Engine |
| 🔬 Mol2ChemFig | 5001 | Flask API Wrapper |
| 🌐 PubChem | 5002 | 3D Models & Images |

---

## ✅ Test It Works

### Test 1: MoleculeViewer
Open in browser:
```
http://localhost:5000/img/smiles?smiles=CCO
```
You should see an SVG image of ethanol.

### Test 2: PubChem
Open in browser:
```
http://localhost:5002/pubchem/img/histamine
```
You should see a PNG image of histamine.

### Test 3: Extension
1. Load extension: `chrome://extensions/` → Load unpacked → Select `chem-extension` folder
2. Open `TEST_EXTENSION.html`
3. Type: `chem:benzene:`
4. Should auto-replace with chemical structure

---

## 🛑 Stopping Servers

**Double-click this file:**
```
stop_all.bat
```

---

## 📊 Check Status

**Double-click this file:**
```
status.bat
```

Shows which servers are running.

---

## 🔧 Troubleshooting

### Server Won't Start?

1. **Port already in use:**
   ```batch
   stop_all.bat
   start_all.bat
   ```

2. **Docker not running:**
   - Start Docker Desktop
   - Wait for it to fully load
   - Run `start_all.bat` again

3. **Python/Node not found:**
   - Install Python and Node.js
   - Restart your computer
   - Try again

### Extension Not Working?

1. Check all servers are running (green badges in launcher)
2. Reload extension in Chrome
3. Check browser console for errors

---

## 📁 Files Overview

| File | Purpose |
|------|---------|
| `launcher.html` | Visual control panel |
| `start_all.bat` | Start all servers |
| `stop_all.bat` | Stop all servers |
| `status.bat` | Check server status |
| `start_moleculeviewer.bat` | Start only MoleculeViewer |
| `start_docker.bat` | Start only Docker backend |
| `start_mol2chemfig.bat` | Start only Mol2ChemFig |
| `start_pubchem.bat` | Start only PubChem |

---

## 🎓 Next Steps

### For Users:
1. ✅ Start servers with `start_all.bat`
2. ✅ Load Chrome extension
3. ✅ Test on any webpage by typing `chem:compound_name:`

### For Developers:
1. ✅ Read `LAUNCHER_README.md` for detailed docs
2. ✅ Check API endpoints in launcher
3. ✅ Explore test files (`test_*.html`)

---

## 💡 Pro Tips

- **Auto-start on login:** Create a shortcut to `start_all.bat` in Windows Startup folder
- **Multiple tests:** Keep launcher open to monitor all servers
- **Quick restart:** `stop_all.bat` → `start_all.bat`
- **View logs:** Check terminal windows opened by start scripts

---

## 📚 Full Documentation

For complete documentation, see:
- `LAUNCHER_README.md` - Complete launcher guide
- `.claude.md` - Project overview
- Server-specific READMEs in each directory

---

## ⚡ Common Commands

```batch
# Start everything
start_all.bat

# Stop everything
stop_all.bat

# Check status
status.bat

# Open launcher
start launcher.html
```

---

**Ready to use ChemParser!** 🎉

Need help? Check `LAUNCHER_README.md` for detailed troubleshooting.
