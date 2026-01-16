/**
 * Firefox Build Script for ChemistryLaTeX
 * FOR MOZILLA REVIEWERS - Standalone build from source
 * 
 * This script builds the Firefox extension directly from source files.
 * No Chrome build required.
 * 
 * Run: npm install && node build-firefox.js
 */

const fs = require('fs');
const path = require('path');
const esbuild = require('esbuild');

const SOURCE_DIR = __dirname;
const OUTPUT_DIR = path.join(__dirname, 'dist');

// Firefox extension configuration
const FIREFOX_CONFIG = {
  extensionId: '{fddce64f-6901-4661-a45a-c1c43d0eeada}',
  minFirefoxVersion: '140.0'
};

console.log('🦊 Building ChemistryLaTeX for Firefox...\n');

/**
 * Copy file with directory creation
 */
function copyFile(src, dest) {
  const destDir = path.dirname(dest);
  if (!fs.existsSync(destDir)) fs.mkdirSync(destDir, { recursive: true });
  fs.copyFileSync(src, dest);
}

/**
 * Copy directory recursively
 */
function copyDir(src, dest) {
  if (!fs.existsSync(dest)) fs.mkdirSync(dest, { recursive: true });
  for (const entry of fs.readdirSync(src, { withFileTypes: true })) {
    const srcPath = path.join(src, entry.name);
    const destPath = path.join(dest, entry.name);
    entry.isDirectory() ? copyDir(srcPath, destPath) : copyFile(srcPath, destPath);
  }
}

/**
 * CSP Bypass Wrapper for Firefox
 * Prepended to content.min.js to route fetch() through background script
 */
const FIREFOX_CSP_BYPASS_CODE = `
/* Firefox CSP Bypass - Routes fetch through background script */
(function(){
  const _originalFetch = window.fetch;
  const CHEMSERVER = 'server-chemistryrenderer.vercel.app';
  
  window.fetch = async function(url, opts) {
    const urlStr = typeof url === 'string' ? url : url.toString();
    
    if (urlStr.includes(CHEMSERVER)) {
      return new Promise((resolve, reject) => {
        if (typeof chrome !== 'undefined' && chrome.runtime && chrome.runtime.sendMessage) {
          chrome.runtime.sendMessage(
            { type: 'FETCH_TEXT', url: urlStr, options: opts || {} },
            (response) => {
              if (chrome.runtime.lastError) {
                _originalFetch(url, opts).then(resolve).catch(reject);
                return;
              }
              if (response && response.success) {
                const headers = new Headers();
                headers.set('content-type', 'image/svg+xml');
                resolve(new Response(response.data, { status: 200, statusText: 'OK', headers }));
              } else {
                reject(new Error(response?.error || 'Background fetch failed'));
              }
            }
          );
        } else {
          _originalFetch(url, opts).then(resolve).catch(reject);
        }
      });
    }
    return _originalFetch(url, opts);
  };
})();
`;

/**
 * Fix popup.html for Firefox CSP
 */
function fixPopupHtml(content) {
  // Remove inline event handlers
  content = content.replace(/\s+onmouseenter="[^"]*"/gi, '');
  content = content.replace(/\s+onmouseleave="[^"]*"/gi, '');
  content = content.replace(/\s+onclick="[^"]*"/gi, '');

  // Add CSS hover styles
  const hoverStyles = `
    /* Firefox CSP Fix */
    #clearCacheBtn:hover { background: #e74c3c !important; color: white !important; }
    #reloadAllBtn:hover { background: #2c3e50 !important; color: white !important; }
  `;
  const lastStyleIndex = content.lastIndexOf('</style>');
  if (lastStyleIndex > -1) {
    content = content.slice(0, lastStyleIndex) + hoverStyles + content.slice(lastStyleIndex);
  }

  return content;
}

/**
 * Create Firefox manifest from source manifest
 */
function createFirefoxManifest(sourceManifest) {
  const manifest = JSON.parse(sourceManifest);

  // Add Firefox-specific settings
  manifest.browser_specific_settings = {
    gecko: {
      id: FIREFOX_CONFIG.extensionId,
      strict_min_version: FIREFOX_CONFIG.minFirefoxVersion,
      data_collection_permissions: { required: ["none"] }
    }
  };

  // Convert service_worker to scripts array
  if (manifest.background?.service_worker) {
    manifest.background = { scripts: [manifest.background.service_worker] };
  }

  return manifest;
}

/**
 * Main build function
 */
async function build() {
  try {
    // Clean output
    if (fs.existsSync(OUTPUT_DIR)) {
      fs.rmSync(OUTPUT_DIR, { recursive: true });
    }
    fs.mkdirSync(OUTPUT_DIR, { recursive: true });

    // === Build JS files with esbuild (minified, no custom obfuscation) ===
    console.log('📦 Building JavaScript files...');

    // content.min.js
    await esbuild.build({
      entryPoints: [path.join(SOURCE_DIR, 'content.js')],
      bundle: false,
      minify: true,
      drop: ['console', 'debugger'],
      outfile: path.join(OUTPUT_DIR, 'content.min.js'),
      target: ['firefox140'],
      legalComments: 'none',
    });
    // Prepend CSP bypass wrapper
    let contentJs = fs.readFileSync(path.join(OUTPUT_DIR, 'content.min.js'), 'utf8');
    fs.writeFileSync(path.join(OUTPUT_DIR, 'content.min.js'), FIREFOX_CSP_BYPASS_CODE + '\n' + contentJs);
    console.log('   ✅ content.min.js');

    // popup.min.js
    await esbuild.build({
      entryPoints: [path.join(SOURCE_DIR, 'popup.js')],
      bundle: false,
      minify: true,
      drop: ['console', 'debugger'],
      outfile: path.join(OUTPUT_DIR, 'popup.min.js'),
      target: ['firefox140'],
      legalComments: 'none',
    });
    console.log('   ✅ popup.min.js');

    // background.min.js
    await esbuild.build({
      entryPoints: [path.join(SOURCE_DIR, 'background.js')],
      bundle: false,
      minify: true,
      drop: ['console', 'debugger'],
      outfile: path.join(OUTPUT_DIR, 'background.min.js'),
      target: ['firefox140'],
      legalComments: 'none',
    });
    console.log('   ✅ background.min.js');

    // size-controls.min.js
    await esbuild.build({
      entryPoints: [path.join(SOURCE_DIR, 'size-controls.js')],
      bundle: false,
      minify: true,
      drop: ['console', 'debugger'],
      outfile: path.join(OUTPUT_DIR, 'size-controls.min.js'),
      target: ['firefox140'],
      legalComments: 'none',
    });
    console.log('   ✅ size-controls.min.js');

    // === Physics files ===
    console.log('📦 Building physics files...');
    const physicsDir = path.join(OUTPUT_DIR, 'physics');
    fs.mkdirSync(physicsDir, { recursive: true });

    // Minify physics JS files
    for (const file of ['popup-aninmation.js', 'physics-colliders.js']) {
      const srcPath = path.join(SOURCE_DIR, 'physics', file);
      if (fs.existsSync(srcPath)) {
        await esbuild.build({
          entryPoints: [srcPath],
          bundle: false,
          minify: true,
          outfile: path.join(physicsDir, file.replace('.js', '.min.js')),
          target: ['firefox140'],
          legalComments: 'none',
        });
      }
    }
    // Copy vendor files as-is
    for (const file of ['vendor-matter.min.js', 'vendor-svg-path-properties-wrapper.js', 'vendor-svg-path-properties.cjs']) {
      const srcPath = path.join(SOURCE_DIR, 'physics', file);
      if (fs.existsSync(srcPath)) {
        copyFile(srcPath, path.join(physicsDir, file));
      }
    }
    console.log('   ✅ physics files');

    // === Assets ===
    console.log('📦 Copying assets...');
    copyDir(path.join(SOURCE_DIR, 'assets'), path.join(OUTPUT_DIR, 'assets'));
    console.log('   ✅ assets');

    // === popup.html ===
    console.log('📦 Processing popup.html...');
    let popupHtml = fs.readFileSync(path.join(SOURCE_DIR, 'popup.html'), 'utf8');
    popupHtml = popupHtml.replace(/popup\.js/g, 'popup.min.js');
    popupHtml = popupHtml.replace(/physics-colliders\.js/g, 'physics-colliders.min.js');
    popupHtml = popupHtml.replace(/popup-aninmation\.js/g, 'popup-aninmation.min.js');
    popupHtml = fixPopupHtml(popupHtml);
    fs.writeFileSync(path.join(OUTPUT_DIR, 'popup.html'), popupHtml);
    console.log('   ✅ popup.html');

    // === manifest.json ===
    console.log('📦 Creating Firefox manifest...');
    const sourceManifest = fs.readFileSync(path.join(SOURCE_DIR, 'manifest.json'), 'utf8');
    const firefoxManifest = createFirefoxManifest(sourceManifest);
    fs.writeFileSync(path.join(OUTPUT_DIR, 'manifest.json'), JSON.stringify(firefoxManifest, null, 2));
    console.log('   ✅ manifest.json');

    // === Summary ===
    console.log('\n🦊 Build complete!');
    console.log(`📁 Output: ${OUTPUT_DIR}`);
    console.log('\n📝 The built extension is in the "dist" folder.');
    console.log('   Compare dist/ contents with the submitted extension ZIP.');

  } catch (e) {
    console.error('❌ Build failed:', e);
    process.exit(1);
  }
}

build();
