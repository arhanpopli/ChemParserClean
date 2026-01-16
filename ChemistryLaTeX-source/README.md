# ChemistryLaTeX - Source Code for Mozilla Review

## Overview
ChemistryLaTeX is a browser extension that renders chemical structures (molecules, proteins, minerals) directly in AI chat interfaces like ChatGPT and Claude.

---

## Build Environment Requirements

- **Operating System:** Windows, macOS, or Linux
- **Node.js:** v18 or later
- **npm:** v9 or later (comes with Node.js)

---

## Build Instructions

### Step 1: Install Dependencies
```bash
npm install
```

This installs `esbuild` (the only build dependency).

### Step 2: Build the Firefox Extension
```bash
node build-firefox.js
```

### Step 3: Output Location
The built extension will be in the `dist/` folder:
```
dist/
├── content.min.js
├── popup.min.js
├── background.min.js
├── size-controls.min.js
├── popup.html
├── manifest.json
├── assets/
│   ├── icon16.png
│   ├── icon48.png
│   ├── icon128.png
│   └── MetaMask-mUSD-Icon.svg
└── physics/
    ├── popup-aninmation.min.js
    ├── physics-colliders.min.js
    ├── vendor-matter.min.js
    └── ...
```

---

## Comparing Built Output to Submitted Extension

After running the build, compare the `dist/` folder contents with the extension ZIP.

### Important Notes:

1. **JavaScript files** - Should be functionally identical. The minified code logic should match.

2. **manifest.json** - May have minor formatting differences. Key settings (permissions, scripts) should match.

3. **popup.html** - Should be nearly identical (inline handlers removed, script refs updated).

### Key verification points:
- The minified JS contains the same code logic as the source
- No additional hidden code in the extension
- CSP bypass wrapper is clearly identifiable at the start of content.min.js

---

## Build Process Details

### What the build script does:

1. **Minifies JS files with esbuild**
   - Standard minification only
   - NO custom obfuscation or name remapping
   - Original variable names are preserved (just shortened by esbuild)

2. **Prepends Firefox CSP bypass to content.min.js**
   - Required because Firefox content scripts are subject to page CSP
   - Routes fetch() calls through the background script

3. **Processes popup.html**
   - Updates script references to use .min.js versions
   - Removes inline event handlers (blocked by Firefox CSP)
   - Adds CSS hover styles as replacement

4. **Creates Firefox-compatible manifest.json**
   - Adds `browser_specific_settings.gecko` section
   - Converts `service_worker` to `scripts` array

---

## File Mapping

| Source File | → | Built Output |
|-------------|---|--------------|
| content.js | → | dist/content.min.js |
| popup.js | → | dist/popup.min.js |
| background.js | → | dist/background.min.js |
| size-controls.js | → | dist/size-controls.min.js |
| popup.html | → | dist/popup.html (modified) |
| manifest.json | → | dist/manifest.json (modified) |
| physics/*.js | → | dist/physics/*.min.js |
| assets/* | → | dist/assets/* (copied) |

---

## Third-Party Libraries

| Library | File | Source | License |
|---------|------|--------|---------|
| Matter.js | physics/vendor-matter.min.js | https://brm.io/matter-js/ | MIT |
| svg-path-properties | physics/vendor-svg-path-properties.cjs | npm package | MIT |

These are included as-is (already minified by their authors).

---

## External Server

The extension communicates with `https://server-chemistryrenderer.vercel.app` to:
- Render chemical structures as SVG images
- Look up molecule data from PubChem, RCSB PDB, COD

This server code is NOT part of the extension and is not included in this source package.

---

## Contact

Email: quintessenlabs@gmail.com
