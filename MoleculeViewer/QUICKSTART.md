# MoleculeViewer - Quick Start Guide

## 🚀 Starting the Server

1. **Double-click `START.bat`** in this folder
2. A command window will open with the server running
3. **Keep that window open** while using the app

## 🌐 Using the Web Interface

1. Open your browser to: **http://localhost:5000**
2. Enter a SMILES string (e.g., `CCO`, `c1ccccc1`, `CC(=O)O`)
3. Click "Convert to SVG"
4. Your molecule will be rendered!

## ⚠️ Troubleshooting

**If you see "failed to fetch":**
- Make sure the black command window from `START.bat` is still open
- Refresh your browser page
- Check that nothing else is using port 5000

**To stop the server:**
- Close the black command window, or
- Press `Ctrl+C` in the command window

## 📝 What Does START.bat Do?

It runs this command:
```
python -m flask --app app.api run --host=0.0.0.0 --port=5000
```

This starts the Flask backend that converts SMILES to SVG images.

## 🎨 Features

- Convert SMILES to SVG molecules
- Lookup chemical names (benzene, aspirin, etc.)
- Customize rendering options (show carbons, methyls, etc.)
- 24-hour caching of generated molecules
- Dark mode support

---

**That's it!** Just double-click `START.bat` and you're ready to go.
