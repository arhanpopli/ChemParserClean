/**
 * Create complete source code package for Mozilla review
 * Includes everything needed to reproduce the Firefox build
 */

const fs = require('fs');
const path = require('path');

const SOURCE_DIR = __dirname;
const OUTPUT_DIR = path.join(path.dirname(__dirname), 'ChemistryLaTeX-source');

// Files to include in source package
const FILES_TO_COPY = [
    // Main extension JS files (unminified source)
    'content.js',
    'popup.js',
    'background.js',
    'size-controls.js',

    // HTML and config
    'popup.html',
    'manifest.json',

    // Physics files for popup animation
    'physics/popup-aninmation.js',
    'physics/physics-colliders.js',
    'physics/vendor-matter.min.js',
    'physics/vendor-svg-path-properties-wrapper.js',
    'physics/vendor-svg-path-properties.cjs',

    // Assets
    'assets/icon16.png',
    'assets/icon48.png',
    'assets/icon128.png',
    'assets/MetaMask-mUSD-Icon.svg',

    // Documentation
    'USAGE.md',

    // Build files (required by Mozilla)
    'build-firefox.js',
    'package.json',
    'package-lock.json',
];

// Clean output directory
if (fs.existsSync(OUTPUT_DIR)) {
    fs.rmSync(OUTPUT_DIR, { recursive: true });
}
fs.mkdirSync(OUTPUT_DIR, { recursive: true });

console.log('📦 Creating source code package for Mozilla review...\n');

let copied = 0;
for (const file of FILES_TO_COPY) {
    const srcPath = path.join(SOURCE_DIR, file);
    const destPath = path.join(OUTPUT_DIR, file);

    if (fs.existsSync(srcPath)) {
        const destDir = path.dirname(destPath);
        if (!fs.existsSync(destDir)) {
            fs.mkdirSync(destDir, { recursive: true });
        }
        fs.copyFileSync(srcPath, destPath);
        console.log(`   ✅ ${file}`);
        copied++;
    } else {
        console.log(`   ⚠️ ${file} (not found)`);
    }
}

// Create detailed README for Mozilla reviewers
const README = `# ChemistryLaTeX - Source Code for Mozilla Review

## Overview
ChemistryLaTeX is a browser extension that renders chemical structures (molecules, proteins, minerals) directly in AI chat interfaces like ChatGPT and Claude.

## Build Environment Requirements

- **Operating System:** Windows, macOS, or Linux
- **Node.js:** v18 or later
- **npm:** v9 or later (comes with Node.js)

## Build Instructions

### Step 1: Install Dependencies
\`\`\`bash
npm install
\`\`\`

This will install esbuild (the only build dependency).

### Step 2: Build Firefox Extension
\`\`\`bash
node build-firefox.js
\`\`\`

This creates the Firefox extension in a \`ChemistryLaTeX-firefox\` folder.

### Step 3: Compare Output
The built extension will be in \`../ChemistryLaTeX-firefox/\`

## Build Process Details

The build script (build-firefox.js) does the following:

1. **Copies source files** to output directory
2. **Minifies JS files with esbuild** - Standard minification, NO custom obfuscation
3. **Modifies manifest.json** for Firefox compatibility
4. **Injects CSP bypass wrapper** - Required because Firefox content scripts are subject to page CSP

### Minification
We use \`esbuild\` for minification. The minified code:
- Preserves original variable/function names (no custom renaming)
- Removes console.log statements
- Reduces file size

### No Obfuscation
The Firefox build does NOT apply any custom obfuscation. Variable names in the minified output match the source.

## File Structure

| Source File | Built Output | Description |
|-------------|--------------|-------------|
| content.js | content.min.js | Main content script - renders molecules |
| popup.js | popup.min.js | Extension popup UI |
| background.js | background.min.js | Background script for CSP bypass |
| size-controls.js | size-controls.min.js | Image resize controls |
| popup.html | popup.html | Popup HTML (modified for CSP) |
| manifest.json | manifest.json | Firefox-modified manifest |
| physics/*.js | physics/*.min.js | Physics animation (Matter.js based) |

## Third-Party Libraries

| Library | File | Source | License |
|---------|------|--------|---------|
| Matter.js | vendor-matter.min.js | https://brm.io/matter-js/ | MIT |
| svg-path-properties | vendor-svg-path-properties.cjs | npm package | MIT |

## External Dependencies

The extension communicates with \`https://server-chemistryrenderer.vercel.app\` to:
- Render chemical structures as SVG
- Look up molecule data from PubChem, RCSB, COD

This server code is NOT included as it runs externally and is not part of the extension.

## Questions?
Contact: quintessenlabs@gmail.com
`;

fs.writeFileSync(path.join(OUTPUT_DIR, 'README.md'), README);
console.log('   ✅ README.md (build instructions for reviewers)');

console.log(`\n✅ Source package created with ${copied} files`);
console.log(`📁 Output: ${OUTPUT_DIR}`);
