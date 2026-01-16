const esbuild = require('esbuild');
const fs = require('fs');
const path = require('path');

const SOURCE_DIR = path.dirname(__dirname);  // Parent folder (chem-extension-server)
const OUTPUT_DIR = path.join(path.dirname(SOURCE_DIR), 'ChemistryLaTeX-compiled');

// NOTE: Name mapping feature removed for compatibility with existing users' settings

const COPYRIGHT = `/**
 * ChemistryLaTeX
 * Copyright (c) ${new Date().getFullYear()} Arhan Popli. All rights reserved.
 * Unauthorized copying, modification, or distribution of this software is prohibited.
 */
`;

if (fs.existsSync(OUTPUT_DIR)) fs.rmSync(OUTPUT_DIR, { recursive: true });
fs.mkdirSync(OUTPUT_DIR, { recursive: true });

console.log('🧪 Building ChemistryLaTeX...\n');

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


async function build() {
    try {
        // === CONTENT.MIN.JS (with console.log removed) ===
        console.log('📦 Building content.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'content.js')],
            bundle: false,
            minify: true,
            drop: ['console', 'debugger'],
            outfile: path.join(OUTPUT_DIR, 'content.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'content.min.js'));
        console.log('   ✅ content.min.js');


        // === SIZE-CONTROLS.MIN.JS ===
        console.log('📦 Building size-controls.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'size-controls.js')],
            bundle: false,
            minify: true,
            drop: ['console', 'debugger'],
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
            drop: ['console', 'debugger'],
            outfile: path.join(OUTPUT_DIR, 'popup.min.js'),
            target: ['chrome90'],
            legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'popup.min.js'));
        console.log('   ✅ popup.min.js');


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
            drop: ['console', 'debugger'],
            outfile: path.join(OUTPUT_DIR, 'physics', 'physics-colliders.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });

        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'physics', 'popup-aninmation.js')],
            bundle: false, minify: true,
            drop: ['console', 'debugger'],
            outfile: path.join(OUTPUT_DIR, 'physics', 'popup-aninmation.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });
        console.log('   ✅ physics files (no copyright - vendor code)');

        // === BACKGROUND.MIN.JS ===
        console.log('📦 Building background.min.js...');
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'background.js')],
            bundle: false, minify: true,
            drop: ['console', 'debugger'],
            outfile: path.join(OUTPUT_DIR, 'background.min.js'),
            target: ['chrome90'], legalComments: 'none',
        });
        addCopyright(path.join(OUTPUT_DIR, 'background.min.js'));
        console.log('   ✅ background.min.js');


        // === NO VIEWERS NEEDED (using external iframes: MolView, Mol* websites) ===
        console.log('📦 Skipping viewers (using external URLs)...');

        // === ASSETS (copy pre-generated icons from source) ===
        console.log('📦 Copying assets...');
        fs.mkdirSync(path.join(OUTPUT_DIR, 'assets'), { recursive: true });

        // Copy pre-generated icons (created by generate-icons.js)
        // These are blue beaker icons visible on both light and dark toolbars
        for (const size of [16, 48, 128]) {
            const srcIcon = path.join(SOURCE_DIR, 'assets', `icon${size}.png`);
            const destIcon = path.join(OUTPUT_DIR, 'assets', `icon${size}.png`);
            if (fs.existsSync(srcIcon)) {
                copyFile(srcIcon, destIcon);
            } else {
                console.log(`   ⚠️ icon${size}.png not found - run 'node generate-icons.js' first`);
            }
        }
        console.log('   ✅ Copied beaker icons');

        // Copy MetaMask donation icon (used in popup.js donation modal)
        copyFile(path.join(SOURCE_DIR, 'assets', 'MetaMask-mUSD-Icon.svg'), path.join(OUTPUT_DIR, 'assets', 'MetaMask-mUSD-Icon.svg'));
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
        manifest.name = 'ChemistryLaTeX';
        manifest.description = 'Renders chemical structures, proteins, and minerals directly in ChatGPT, Claude, and other AI assistants.';
        manifest.background.service_worker = 'background.min.js';
        manifest.content_scripts[0].js = ['content.min.js', 'size-controls.min.js'];

        // Chrome Manifest V3 requires PNG icons
        // Using white icons for visibility on both light and dark toolbars
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
        console.log('   ✅ manifest.json');

        console.log('\n🎉 Build complete!');
        console.log('📜 Copyright: Arhan Popli (on your code only)');
        console.log('🚫 console.log removed from all files');


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

// Run build, then automatically build dev and Firefox versions
build().then(() => {
    const { execSync } = require('child_process');

    // Build dev version
    console.log('\n' + '='.repeat(50));
    console.log('📦 Also building dev version...\n');
    try {
        execSync('node build-dev.js', {
            cwd: __dirname,
            stdio: 'inherit'
        });
    } catch (e) {
        console.error('❌ Dev build failed:', e.message);
    }

    // Build Firefox version
    console.log('\n' + '='.repeat(50));
    console.log('🦊 Also building Firefox version...\n');
    try {
        execSync('node build-firefox.js', {
            cwd: __dirname,
            stdio: 'inherit'
        });
    } catch (e) {
        console.error('❌ Firefox build failed:', e.message);
    }

    console.log('\n' + '='.repeat(50));
    console.log('✅ All builds complete! (Chrome, Dev, Firefox)');
    console.log('='.repeat(50));
});
