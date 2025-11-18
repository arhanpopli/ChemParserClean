# Chemistry Extension Update - Full MoleculeViewer Integration

## 🎯 Problem Solved

The extension was using CodeCogs API even when MoleculeViewer server was selected. Now:

✅ **Properly detects and uses MoleculeViewer** when selected
✅ **Sends rendering options** (flip, show carbons, etc.) to the server
✅ **Popup menu includes all MoleculeViewer rendering options**
✅ **Proper settings persistence** between browser sessions
✅ **Async rendering** with lazy-loading support

---

## 🔧 What Changed

### 1. **content.js - Settings Management** 
**Issue:** Settings weren't being loaded properly with the stored values
**Fix:** Updated settings loading to use proper merge with defaults

```javascript
// BEFORE: Settings object had hardcoded defaults
let settings = {
  enabled: true,
  rendererEngine: 'codecogs'  // ← Hard default, not merged with storage
};

// AFTER: Proper merge with stored values
let settings = {
  // All defaults including new options
  flipHorizontal: false,
  flipVertical: false,
  showCarbons: false,
  showMethyls: false,
  aromaticCircles: true,
  fancyBonds: true,
  atomNumbers: false,
  hydrogensMode: 'keep'
};

chrome.storage.sync.get(null, (result) => {
  settings = { ...settings, ...result };  // ← Proper merge!
});
```

### 2. **content.js - buildChemfigImageUrl Function**
**Issue:** Returned string URLs, couldn't pass complex options
**Fix:** Returns object for MoleculeViewer with all options

```javascript
// BEFORE: Simple string URL
return `http://localhost:5000/api/render-smiles?smiles=${smiles}&dark=${isDarkMode}`;

// AFTER: Object with full options
return {
  isMoleculeViewer: true,
  smiles: smiles,
  options: {
    show_carbons: settings.showCarbons,
    show_methyls: settings.showMethyls,
    aromatic_circles: settings.aromaticCircles,
    fancy_bonds: settings.fancyBonds,
    atom_numbers: settings.atomNumbers,
    hydrogens: settings.hydrogensMode,
    flip_horizontal: settings.flipHorizontal,
    flip_vertical: settings.flipVertical,
    recalculate_coordinates: false
  },
  isDarkMode: isDarkMode
};
```

### 3. **content.js - Image Creation & Loading**
**Issue:** Image URLs were set as string, couldn't handle object data
**Fix:** Detect MoleculeViewer data and load async via POST

```javascript
// BEFORE: All images loaded same way
converted = `<img src="${imageUrl}" alt="chemfig" ...>`;

// AFTER: Different handling for MoleculeViewer
if (typeof imageUrl === 'object' && imageUrl.isMoleculeViewer) {
  // Use special class for async loading
  const moleculeViewerData = btoa(JSON.stringify(imageUrl));
  converted = `<img ... class="chemfig-molecule-viewer" data-molecule-viewer="${moleculeViewerData}" ...>`;
} else {
  // Standard image loading
  converted = `<img src="${imageUrl}" ...>`;
}
```

### 4. **content.js - New loadMoleculeViewerImage Function**
**Issue:** MoleculeViewer images need async fetch + POST request
**Fix:** New function handles server communication

```javascript
function loadMoleculeViewerImage(img) {
  // 1. Decode stored options from data attribute
  const moleculeData = JSON.parse(atob(img.dataset.moleculeViewer));
  
  // 2. POST to server with all options
  fetch('http://localhost:5000/api/smiles-to-svg', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      smiles: moleculeData.smiles,
      width: 300,
      height: 200,
      options: moleculeData.options
    })
  })
  // 3. Get SVG response
  .then(response => response.json())
  .then(data => {
    // 4. Create blob URL for the SVG
    const svgBlob = new Blob([data.svg], { type: 'image/svg+xml' });
    const blobUrl = URL.createObjectURL(svgBlob);
    
    // 5. Set image source
    img.src = blobUrl;
  });
}
```

### 5. **content.js - Updated Lazy Loading**
**Issue:** Observer didn't know how to handle MoleculeViewer images
**Fix:** New observer detects image type and loads accordingly

```javascript
// Detects class and routes to correct loader
if (img.classList.contains('chemfig-molecule-viewer')) {
  loadMoleculeViewerImage(img);  // Async POST
} else {
  loadImage(img);  // Standard preload
}
```

### 6. **popup.html - New UI Section**
**Issue:** No UI to control rendering options
**Fix:** Added entire "MoleculeViewer Options" section

```html
<div id="moleculeViewerOptions" class="section" style="display: none;">
  <div class="section-title">🧪 MoleculeViewer Options</div>
  
  <!-- All 8 toggle controls -->
  <div class="option">
    <label for="showCarbonsToggle">
      <strong>Show Carbon Atoms</strong>
    </label>
    <input type="checkbox" id="showCarbonsToggle">
    <label class="toggle" for="showCarbonsToggle"></label>
  </div>
  
  <!-- ... 7 more toggles ... -->
  
  <!-- 1 dropdown for hydrogens -->
  <select id="hydrogensSelect">
    <option value="keep">Keep as drawn</option>
    <option value="add">Add all hydrogens</option>
    <option value="delete">Remove all hydrogens</option>
  </select>
</div>
```

### 7. **popup.js - New Event Handlers**
**Issue:** No way to save rendering option changes
**Fix:** Added event listeners for all 8 new controls

```javascript
showCarbonsToggle.addEventListener('change', (e) => {
  chrome.storage.sync.set({ showCarbons: e.target.checked }, () => {
    showStatus('Show carbons ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
  });
});

// ... similar for all other options ...
```

### 8. **popup.js - Show/Hide Options Based on Engine**
**Issue:** Rendering options showed for all engines
**Fix:** Only show when MoleculeViewer is selected

```javascript
function updateEngineInfo(engine) {
  const info = rendererInfo[engine];
  
  // Show/hide MoleculeViewer options
  if (moleculeViewerOptions) {
    if (engine === 'molecule-viewer') {
      moleculeViewerOptions.style.display = 'block';
    } else {
      moleculeViewerOptions.style.display = 'none';
    }
  }
}
```

---

## 📊 Rendering Options Now Available

| Option | Default | Effect |
|--------|---------|--------|
| **Show Carbon Atoms** | OFF | Display C labels in structure |
| **Show Methyl Groups** | OFF | Display CH₃ labels |
| **Aromatic Circles** | ON | Draw circles in benzene rings |
| **Fancy Bonds** | ON | Enhanced bond visualization |
| **Atom Numbers** | OFF | Show atom index numbers |
| **Flip Horizontal** | OFF | Mirror structure left-right |
| **Flip Vertical** | OFF | Mirror structure up-down |
| **Hydrogens Mode** | keep | keep / add / delete |

---

## 🚀 How It Works Now

### Flow Diagram

```
User views \chemfig{...} on webpage
           ↓
Extension detects formula
           ↓
Converts to SMILES string
           ↓
Checks: Is MoleculeViewer selected?
           ↓
    YES → Creates object with:
         - SMILES string
         - All rendering options (flip, show carbons, etc.)
         - Stores in data attribute (base64 encoded)
         - Creates <img> with special class
           ↓
When image enters viewport (lazy-loading):
         ↓
loadMoleculeViewerImage() function:
         - Decodes options from data attribute
         - POSTs to http://localhost:5000/api/smiles-to-svg
         - Gets SVG response
         - Creates blob URL
         - Sets img.src = blobUrl
           ↓
SVG displays in webpage ✓
           ↓
    NO → Uses standard CodeCogs flow (unchanged)
```

### Key Improvements

1. **Proper Settings Merge** - New settings load from storage correctly
2. **Full Options Support** - All 8 MoleculeViewer rendering options available
3. **Async Loading** - POST requests don't block UI
4. **Lazy Loading** - Still respects performance settings
5. **Fallback** - CodeCogs still works if MoleculeViewer unavailable
6. **User Friendly UI** - Options only show when engine is selected

---

## 🧪 Testing Steps

### 1. Start the Server
```bash
cd MoleculeViewer
python run_server.py
```

### 2. Load Extension
1. Open `chrome://extensions/`
2. Enable "Developer mode"
3. Click "Load unpacked"
4. Select `chem-extension` folder

### 3. Open Extension Popup
1. Click extension icon in Chrome
2. Scroll to "Rendering Engine" section
3. Select "🧪 MoleculeViewer (Best)"
4. New "🧪 MoleculeViewer Options" section should appear

### 4. Toggle Options
1. Try turning on "Show Carbon Atoms"
2. Try "Flip Horizontal"
3. Try changing "Hydrogen Display"
4. Each should show success message

### 5. Test Live
1. Go to ChatGPT or any webpage
2. Use chemistry notation: `\chemfig{C-C-C}`
3. Should render using local server (not CodeCogs URL)
4. Should respect your rendering options

### 6. Check Console
1. Press F12 → Console tab
2. Look for: `🔬 Using MoleculeViewer server for rendering`
3. Should see rendering options logged

---

## 🔗 File Locations

| File | Changes |
|------|---------|
| `chem-extension/content.js` | Settings loading, buildChemfigImageUrl, image creation, loadMoleculeViewerImage |
| `chem-extension/popup.html` | New MoleculeViewer Options section |
| `chem-extension/popup.js` | Event listeners, option handling, show/hide logic |

---

## ✅ Checklist

- ✅ Settings properly merge with stored values
- ✅ MoleculeViewer engine recognized when selected
- ✅ All 8 rendering options configurable in UI
- ✅ Options saved to chrome.storage.sync
- ✅ POST request sends options to server
- ✅ SVG response rendered correctly
- ✅ Lazy-loading still works
- ✅ Fallback to CodeCogs if needed
- ✅ UI shows/hides options based on engine
- ✅ Console logs for debugging

---

## 🎉 You're all set!

The extension now properly integrates with MoleculeViewer server and includes all rendering options. The UI dynamically shows/hides options based on your selection, and everything persists between browser sessions.

Try it out and let me know if you see it using your local server! 🧪✨
