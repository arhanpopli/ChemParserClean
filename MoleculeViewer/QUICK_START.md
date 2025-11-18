# MoleculeViewer - Quick Start Guide

## How to Use

### Step 1: Start the Backend Server

Double-click this file:
```
START_BACKEND.bat
```

A black window will appear. **LEAVE IT OPEN!** Do not close this window.

You should see:
```
======================================================================
MOLECULEVIEWER BACKEND SERVER
======================================================================
Starting Flask server...
Listen on: http://localhost:5000
...
Running on all addresses (0.0.0.0)
Press CTRL+C to quit
```

### Step 2: Reload Chrome Extension

1. Open Chrome browser
2. Go to: `chrome://extensions/`
3. Find **"Chemistry Renderer"** extension
4. Click the **RELOAD** button (circular arrow icon)

### Step 3: Test in ChatGPT

1. Open ChatGPT
2. Type: `chem:benzene`
3. The molecule should appear inline!

### Step 4: Enjoy Dark Mode Support

Toggle dark mode in ChatGPT settings - the molecule lines automatically invert!

---

## Troubleshooting

### "Failed to fetch" error?
- Make sure the black backend window is still open
- If closed, double-click `START_BACKEND.bat` again

### Nothing showing up?
- Check that extension is reloaded at `chrome://extensions/`
- Make sure you're using a **NEW** ChatGPT conversation

### Want to stop the backend?
- Click the X button on the black window
- Or press `CTRL+C` in the window

---

## Files

- `START_BACKEND.bat` - Click this to start the backend
- `start_server_simple.py` - Backend server (called by the batch file)
- `app/api.py` - Flask backend logic
- `.env` - Configuration (backend URL)

That's it! Enjoy inline molecules! 🧬
