const esbuild = require('esbuild');
const fs = require('fs');
const path = require('path');

const SOURCE_DIR = path.dirname(__dirname);  // Parent folder (chem-extension-server)
const OUTPUT_DIR = path.join(path.dirname(SOURCE_DIR), 'dev_ChemistryLaTeX');

// DEV BUILD: No name remapping - just minified with console logs for debugging

const COPYRIGHT = `/**
 * ChemistryLaTeX (Developer Build)
 * Copyright (c) ${new Date().getFullYear()} Arhan Popli. All rights reserved.
 * DEBUG MODE: Console logs ENABLED
 */
`;

if (fs.existsSync(OUTPUT_DIR)) fs.rmSync(OUTPUT_DIR, { recursive: true });
fs.mkdirSync(OUTPUT_DIR, { recursive: true });

console.log('🧪 Building ChemistryLaTeX (DEV MODE - console logs ENABLED)...\n');

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

// DEV BUILD: Name remapping disabled - use production build for obfuscation


async function build() {
    try {
        // === CONTENT.MIN.JS (WITHOUT console.log removal for dev) ===
        console.log('📦 Building content.min.js (with console logs)...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'content.js')],
            bundle: false,
            minify: true,
            // NO drop: ['console', 'debugger'] - keep logs for debugging
            outfile: path.join(OUTPUT_DIR, 'content.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'content.min.js'));
        console.log('   ✅ content.min.js (console logs ENABLED, no name remapping)');


        // === SIZE-CONTROLS.MIN.JS ===
        console.log('📦 Building size-controls.min.js (with console logs)...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'size-controls.js')],
            bundle: false,
            minify: true,
            // NO drop: ['console', 'debugger']
            outfile: path.join(OUTPUT_DIR, 'size-controls.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'size-controls.min.js'));
        console.log('   ✅ size-controls.min.js');

        // === POPUP.MIN.JS ===
        console.log('📦 Building popup.min.js (with console logs)...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'popup.js')],
            bundle: false,
            minify: true,
            // NO drop: ['console', 'debugger']
            outfile: path.join(OUTPUT_DIR, 'popup.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'popup.min.js'));
        console.log('   ✅ popup.min.js (console logs ENABLED, no name remapping)');


        // === PHYSICS ===
        console.log('📦 Building physics...');
        fs.mkdirSync(path.join(OUTPUT_DIR, 'physics'), { recursive: true });

        copyFile(path.join(SOURCE_DIR, 'physics', 'vendor-matter.min.js'), path.join(OUTPUT_DIR, 'physics', 'vendor-matter.min.js'));
        copyFile(path.join(SOURCE_DIR, 'physics', 'vendor-svg-path-properties-wrapper.js'), path.join(OUTPUT_DIR, 'physics', 'vendor-svg-path-properties-wrapper.js'));
        copyFile(path.join(SOURCE_DIR, 'physics', 'vendor-svg-path-properties.cjs'), path.join(OUTPUT_DIR, 'physics', 'vendor-svg-path-properties.cjs'));

        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'physics', 'physics-colliders.js')],
            bundle: false, minify: true,
            // NO drop: ['console', 'debugger']
            outfile: path.join(OUTPUT_DIR, 'physics', 'physics-colliders.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });

        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'physics', 'popup-aninmation.js')],
            bundle: false, minify: true,
            // NO drop: ['console', 'debugger']
            outfile: path.join(OUTPUT_DIR, 'physics', 'popup-aninmation.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });
        console.log('   ✅ physics files');

        // === BACKGROUND.MIN.JS ===
        console.log('📦 Building background.min.js (with console logs)...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'background.js')],
            bundle: false, minify: true,
            // NO drop: ['console', 'debugger']
            outfile: path.join(OUTPUT_DIR, 'background.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'background.min.js'));
        console.log('   ✅ background.min.js (console logs ENABLED, no name remapping)');


        // === ASSETS ===
        console.log('📦 Copying assets from source...');
        fs.mkdirSync(path.join(OUTPUT_DIR, 'assets'), { recursive: true });

        // Copy all icons from source assets directory
        const assetsDir = path.join(SOURCE_DIR, 'assets');
        if (fs.existsSync(assetsDir)) {
            for (const file of fs.readdirSync(assetsDir)) {
                copyFile(path.join(assetsDir, file), path.join(OUTPUT_DIR, 'assets', file));
            }
            console.log('   ✅ assets copied from source');
        } else {
            console.log('   ⚠️ No assets directory found in source');
        }

        // === POPUP.HTML ===
        console.log('📦 Updating popup.html...');
        let html = fs.readFileSync(path.join(SOURCE_DIR, 'popup.html'), 'utf8');
        html = html.replace(/popup\.js/g, 'popup.min.js');
        html = html.replace(/physics-colliders\.js/g, 'physics-colliders.min.js');
        html = html.replace(/popup-aninmation\.js/g, 'popup-aninmation.min.js');
        fs.writeFileSync(path.join(OUTPUT_DIR, 'popup.html'), html);
        console.log('   ✅ popup.html');

        // === USAGE.md ===
        if (fs.existsSync(path.join(SOURCE_DIR, 'USAGE.md'))) {
            copyFile(path.join(SOURCE_DIR, 'USAGE.md'), path.join(OUTPUT_DIR, 'USAGE.md'));
        }

        // === MANIFEST (with DEV name) ===
        console.log('📦 Creating manifest.json (DEV)...');
        const manifest = JSON.parse(fs.readFileSync(path.join(SOURCE_DIR, 'manifest.json'), 'utf8'));
        manifest.name = 'dev_ChemistryLaTeX';  // DEV NAME
        manifest.description = '[DEV] Renders chemical structures with console logs enabled for debugging.';
        manifest.background.service_worker = 'background.min.js';
        manifest.content_scripts[0].js = ['content.min.js', 'size-controls.min.js'];

        manifest.icons = {
            "16": "assets/icon16.png",
            "48": "assets/icon48.png",
            "128": "assets/icon128.png"
        };
        manifest.action.default_icon = manifest.icons;
        manifest.web_accessible_resources[0].resources = [
            "USAGE.md"
        ];
        fs.writeFileSync(path.join(OUTPUT_DIR, 'manifest.json'), JSON.stringify(manifest, null, 2));
        console.log('   ✅ manifest.json (DEV)');

        console.log('\n🎉 DEV Build complete!');
        console.log('📜 Copyright: Arhan Popli');
        console.log('🔍 Console logs ENABLED for debugging');
        console.log('❌ NO name remapping (pure debug build)');


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
