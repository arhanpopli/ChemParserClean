# Size Controls Integration Instructions

## Step 1: Add Size Controls Code to content.js

Insert the entire contents of `size-controls.js` into `content.js` after the performance monitoring system (around line 142, after `window.chemRendererDebug`).

## Step 2: Modify Image Replacement Calls

Replace all instances of:
```javascript
img.parentNode.replaceChild(svgImg, img);
```

With:
```javascript
// Get settings for size controls
chrome.storage.sync.get({
  saveSizePerImage: false,
  saveSizeBySMILES: false
}, async (settings) => {
  await wrapImageWithSizeControls(svgImg, img, moleculeData, settings);
});
```

This affects the following lines in content.js:
- Line 909 (loadMoleculeViewerImage)
- Line 1118 (loadMol2chemfigImage - chemfig fallback)
- Line 1145 (loadMol2chemfigImage - PDF placeholder)
- Line 1215 (loadMol2chemfigImage - data URI)
- Line 1262 (loadMol2chemfigImage - raw SVG)
- Line 1332 (loadMol2chemfigImage - svglink)
- Line 1340 (loadMol2chemfigImage - fallback)
- Line 1344 (loadMol2chemfigImage - error)

## Step 3: Add Settings to popup.html

Add this section to popup.html after the "Developer Mode" section (after line 466):

```html
<!-- Image Size Controls -->
<div class="section">
  <div class="section-title">📏 Image Size Controls</div>

  <div class="option">
    <label for="saveSizePerImageToggle">
      <strong>Save Size Per Page</strong>
      <small>Remember image size for each page separately</small>
    </label>
    <input type="checkbox" id="saveSizePerImageToggle">
    <label class="toggle" for="saveSizePerImageToggle"></label>
  </div>

  <div class="option option-border">
    <label for="saveSizeBySMILESToggle">
      <strong>Save Size By SMILES</strong>
      <small>Use same size for all molecules with same SMILES (overrides per-page setting)</small>
    </label>
    <input type="checkbox" id="saveSizeBySMILESToggle">
    <label class="toggle" for="saveSizeBySMILESToggle"></label>
  </div>

  <div class="info-box">
    <strong>How it works:</strong> Use the up/down arrows in the bottom-left corner of each molecule image to adjust its size. Your preferences will be saved based on the options above.
  </div>
</div>
```

## Step 4: Add Settings Handlers to popup.js

Add this code to popup.js after the devMode toggle event listener (after line 147):

```javascript
// Image size control settings
const saveSizePerImageToggle = document.getElementById('saveSizePerImageToggle');
const saveSizeBySMILESToggle = document.getElementById('saveSizeBySMILESToggle');

if (saveSizePerImageToggle) {
  saveSizePerImageToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ saveSizePerImage: e.target.checked }, () => {
      showStatus('Save size per page ' + (e.target.checked ? 'enabled' : 'disabled') + '.', 'success');
      // Disable SMILES option if per-page is enabled
      if (e.target.checked && saveSizeBySMILESToggle) {
        saveSizeBySMILESToggle.disabled = false;
      }
    });
  });
}

if (saveSizeBySMILESToggle) {
  saveSizeBySMILESToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ saveSizeBySMILES: e.target.checked }, () => {
      showStatus('Save size by SMILES ' + (e.target.checked ? 'enabled' : 'disabled') + '. This will apply the same size to all molecules with the same SMILES.', 'success');
    });
  });
}
```

And add these to the initial settings load (around line 46):

```javascript
saveSizePerImage: false,
saveSizeBySMILES: false
```

And load them (around line 97):

```javascript
if (saveSizePerImageToggle) saveSizePerImageToggle.checked = settings.saveSizePerImage;
if (saveSizeBySMILESToggle) saveSizeBySMILESToggle.checked = settings.saveSizeBySMILES;
```

## Notes

- The size controls will appear when hovering over molecule images
- Size changes are saved to chrome.storage.local
- The "Save Size Per Page" option saves different sizes for the same molecule on different pages
- The "Save Size By SMILES" option saves the same size for all instances of a molecule with the same SMILES across all pages
