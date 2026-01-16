const esbuild = require('esbuild');
const fs = require('fs');
const path = require('path');

const SOURCE_DIR = __dirname;
const OUTPUT_DIR = path.join(path.dirname(__dirname), 'exp_ChemistryLaTeX_Ind');

// Load name mappings for obfuscation-lite
const nameMappings = require('./name-mappings.js');

const COPYRIGHT = `/**
 * ChemistryLaTeX - Independent Client-Side Edition
 * Copyright (c) ${new Date().getFullYear()} Arhan Popli. All rights reserved.
 * Unauthorized copying, modification, or distribution of this software is prohibited.
 */
`;

if (fs.existsSync(OUTPUT_DIR)) fs.rmSync(OUTPUT_DIR, { recursive: true });
fs.mkdirSync(OUTPUT_DIR, { recursive: true });

console.log('🧪 Building exp_ChemistryLaTeX_Ind (Independent Client-Side)...\n');

function copyFile(src, dest) {
    const destDir = path.dirname(dest);
    if (!fs.existsSync(destDir)) fs.mkdirSync(destDir, { recursive: true });
    fs.copyFileSync(src, dest);
}

function copyDir(src, dest) {
    if (!fs.existsSync(dest)) fs.mkdirSync(dest, { recursive: true });
    for (const entry of fs.readdirSync(src, { withFileTypes: true })) {
        const srcPath = path.join(src, entry.name);
        const destPath = path.join(dest, entry.name);
        entry.isDirectory() ? copyDir(srcPath, destPath) : copyFile(srcPath, destPath);
    }
}

function addCopyright(filePath) {
    const content = fs.readFileSync(filePath, 'utf8');
    fs.writeFileSync(filePath, COPYRIGHT + content);
}

/**
 * Apply name mappings to compiled code
 * Replaces variable/function names with fun alternatives
 */
function applyNameMappings(filePath, scope) {
    let content = fs.readFileSync(filePath, 'utf8');
    let changes = 0;

    // Apply variable mappings
    for (const [original, replacement, mappingScope] of nameMappings.variables) {
        if (mappingScope === 'all' || mappingScope === scope) {
            // Use word boundary matching to avoid partial replacements
            const regex = new RegExp(`\\b${original}\\b`, 'g');
            const matches = content.match(regex);
            if (matches) {
                content = content.replace(regex, replacement);
                changes += matches.length;
            }
        }
    }

    // Apply string mappings (for traces like "SmilesDrawer")  
    for (const [original, replacement, mappingScope] of nameMappings.strings) {
        if (mappingScope === 'all' || mappingScope === scope) {
            const regex = new RegExp(original, 'gi');
            const matches = content.match(regex);
            if (matches) {
                content = content.replace(regex, replacement);
                changes += matches.length;
            }
        }
    }

    fs.writeFileSync(filePath, content);
    return changes;
}


async function build() {
    try {
        // === CONTENT.MIN.JS (keep console logs for developers) ===
        console.log('📦 Building content.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'content.js')],
            bundle: false,
            minify: true,
            // Note: NO drop:['console'] - independent variant keeps logs for debugging
            outfile: path.join(OUTPUT_DIR, 'content.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        // Note: NO name remapping for independent variant - developers need readable code
        addCopyright(path.join(OUTPUT_DIR, 'content.min.js'));
        console.log('   ✅ content.min.js (no name remapping - developer build)');


        // === SIZE-CONTROLS.MIN.JS ===
        console.log('📦 Building size-controls.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'size-controls.js')],
            bundle: false,
            minify: true,
            outfile: path.join(OUTPUT_DIR, 'size-controls.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'size-controls.min.js'));
        console.log('   ✅ size-controls.min.js');

        // === POPUP.MIN.JS ===
        console.log('📦 Building popup.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'popup.js')],
            bundle: false,
            minify: true,
            outfile: path.join(OUTPUT_DIR, 'popup.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        // Note: NO name remapping for independent variant
        addCopyright(path.join(OUTPUT_DIR, 'popup.min.js'));
        console.log('   ✅ popup.min.js (no name remapping - developer build)');


        // === PHYSICS (NO copyright - not our code) ===
        console.log('📦 Building physics...');
        fs.mkdirSync(path.join(OUTPUT_DIR, 'physics'), { recursive: true });

        // Copy vendor files as-is (already minified, not our code)
        copyFile(path.join(SOURCE_DIR, 'physics', 'vendor-matter.min.js'), path.join(OUTPUT_DIR, 'physics', 'vendor-matter.min.js'));
        copyFile(path.join(SOURCE_DIR, 'physics', 'vendor-svg-path-properties-wrapper.js'), path.join(OUTPUT_DIR, 'physics', 'vendor-svg-path-properties-wrapper.js'));
        copyFile(path.join(SOURCE_DIR, 'physics', 'vendor-svg-path-properties.cjs'), path.join(OUTPUT_DIR, 'physics', 'vendor-svg-path-properties.cjs'));

        // Minify custom physics files (NO copyright - uses Matter.js which isn't ours)
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'physics', 'physics-colliders.js')],
            bundle: false, minify: true,
            outfile: path.join(OUTPUT_DIR, 'physics', 'physics-colliders.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });

        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'physics', 'popup-aninmation.js')],
            bundle: false, minify: true,
            outfile: path.join(OUTPUT_DIR, 'physics', 'popup-aninmation.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });
        console.log('   ✅ physics files (no copyright - vendor code)');

        // === BACKGROUND.MIN.JS ===
        console.log('📦 Building background.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'background.js')],
            bundle: false, minify: true,
            outfile: path.join(OUTPUT_DIR, 'background.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });
        // Note: NO name remapping for independent variant
        addCopyright(path.join(OUTPUT_DIR, 'background.min.js'));
        console.log('   ✅ background.min.js (no name remapping - developer build)');


        // === NO VIEWERS NEEDED (using external iframes: MolView, Mol* websites) ===
        console.log('📦 Skipping viewers (using external URLs)...');

        // === ASSETS (download professional icons from Iconify CDN) ===
        console.log('📦 Creating assets with beaker icon...');
        fs.mkdirSync(path.join(OUTPUT_DIR, 'assets'), { recursive: true });

        // Download Material Design beaker icon from Iconify CDN and convert to PNG
        // Using white icons which work well on both light and dark toolbars
        const https = require('https');
        const sharp = require('sharp');

        const downloadAndConvertIcon = async (size) => {
            return new Promise((resolve, reject) => {
                // Use white color for visibility on both light and dark toolbars
                const url = `https://api.iconify.design/mdi/beaker.svg?width=${size}&height=${size}&color=white`;
                https.get(url, (res) => {
                    let data = '';
                    res.on('data', chunk => data += chunk);
                    res.on('end', async () => {
                        try {
                            // Convert SVG to PNG using sharp
                            const pngPath = path.join(OUTPUT_DIR, 'assets', `icon${size}.png`);
                            await sharp(Buffer.from(data))
                                .resize(size, size)
                                .png()
                                .toFile(pngPath);
                            resolve(pngPath);
                        } catch (err) {
                            reject(err);
                        }
                    });
                }).on('error', reject);
            });
        };

        // Download and convert all icon sizes to PNG
        try {
            await Promise.all([
                downloadAndConvertIcon(16),
                downloadAndConvertIcon(48),
                downloadAndConvertIcon(128)
            ]);
            console.log('   ✅ Downloaded and converted beaker icons');
        } catch (e) {
            console.log('   ⚠️ Could not download/convert icons:', e.message);
            // Fallback: create PNG from embedded SVG
            const fallbackSvg = (size) => `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}" viewBox="0 0 24 24"><path fill="white" d="M3 3h2v11c0 2.21 1.79 4 4 4s4-1.79 4-4V3h2c.55 0 1-.45 1-1s-.45-1-1-1H3c-.55 0-1 .45-1 1s.45 1 1 1m4 11V3h4v11c0 1.1-.9 2-2 2s-2-.9-2-2m11.5 0c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5m-1 3c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5m3-3c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5m-1.38 5.38c-1.34 1.34-3.5 1.34-4.85 0l4.85-4.85c1.35 1.34 1.35 3.51 0 4.85z"/></svg>`;
            for (const size of [16, 48, 128]) {
                try {
                    await sharp(Buffer.from(fallbackSvg(size)))
                        .resize(size, size)
                        .png()
                        .toFile(path.join(OUTPUT_DIR, 'assets', `icon${size}.png`));
                } catch (err2) {
                    console.log(`   ⚠️ Fallback failed for size ${size}`);
                }
            }
        }

        // Copy MetaMask donation icon (used in popup.js donation modal)
        copyFile(path.join(SOURCE_DIR, 'assets', 'MetaMask-mUSD-Icon.svg'), path.join(OUTPUT_DIR, 'assets', 'MetaMask-mUSD-Icon.svg'));

        // Copy SmilesDrawer (required for client-side rendering)
        copyFile(path.join(SOURCE_DIR, 'assets', 'smiles-drawer.min.js'), path.join(OUTPUT_DIR, 'assets', 'smiles-drawer.min.js'));
        console.log('   ✅ SmilesDrawer copied');

        // Copy MineralNames.js (fallback when COD is down)
        copyFile(path.join(SOURCE_DIR, 'MineralNames.js'), path.join(OUTPUT_DIR, 'MineralNames.js'));
        console.log('   ✅ MineralNames.js copied');

        console.log('   ✅ assets complete');

        // === POPUP.HTML ===
        console.log('📦 Updating popup.html...');
        let html = fs.readFileSync(path.join(SOURCE_DIR, 'popup.html'), 'utf8');
        html = html.replace(/popup\.js/g, 'popup.min.js');
        html = html.replace(/physics-colliders\.js/g, 'physics-colliders.min.js');
        html = html.replace(/popup-aninmation\.js/g, 'popup-aninmation.min.js');
        fs.writeFileSync(path.join(OUTPUT_DIR, 'popup.html'), html);
        console.log('   ✅ popup.html');

        // === USAGE.md only (styles.css not used - styles are inline in popup.html) ===
        if (fs.existsSync(path.join(SOURCE_DIR, 'USAGE.md'))) {
            copyFile(path.join(SOURCE_DIR, 'USAGE.md'), path.join(OUTPUT_DIR, 'USAGE.md'));
        }

        // === MANIFEST ===
        console.log('📦 Creating manifest.json...');
        const manifest = JSON.parse(fs.readFileSync(path.join(SOURCE_DIR, 'manifest.json'), 'utf8'));
        manifest.name = 'exp_ChemistryLaTeX_Ind';
        manifest.description = 'Renders chemical structures, proteins, and minerals in ChatGPT and AI assistants. Independent version - no external server required.';
        manifest.background.service_worker = 'background.min.js';
        manifest.action.default_title = 'exp_ChemistryLaTeX_Ind Settings';
        // Add SmilesDrawer before content.js for client-side rendering
        manifest.content_scripts[0].js = ['assets/smiles-drawer.min.js', 'content.min.js', 'size-controls.min.js'];

        // Chrome Manifest V3 requires PNG icons
        // Using white icons for visibility on both light and dark toolbars
        manifest.icons = {
            "16": "assets/icon16.png",
            "48": "assets/icon48.png",
            "128": "assets/icon128.png"
        };
        manifest.action.default_icon = manifest.icons;
        manifest.web_accessible_resources[0].resources = [
            "assets/smiles-drawer.min.js",
            "MineralNames.js",
            "USAGE.md"
        ];
        fs.writeFileSync(path.join(OUTPUT_DIR, 'manifest.json'), JSON.stringify(manifest, null, 2));
        console.log('   ✅ manifest.json');

        console.log('\n🎉 Build complete!');
        console.log('📜 Copyright: Arhan Popli (on your code only)');
        console.log('🔓 No name remapping - developer readable');
        console.log('📝 console.log KEPT for debugging');


        const getSize = (dir) => {
            let s = 0;
            for (const f of fs.readdirSync(dir, { withFileTypes: true })) {
                const p = path.join(dir, f.name);
                s += f.isDirectory() ? getSize(p) : fs.statSync(p).size;
            }
            return s;
        };
        console.log(`📊 Total: ${(getSize(OUTPUT_DIR) / 1024 / 1024).toFixed(2)} MB`);

    } catch (e) {
        console.error('❌ Build failed:', e);
        process.exit(1);
    }
}

build();
