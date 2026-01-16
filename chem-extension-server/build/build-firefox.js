/**
 * Firefox Build Script for ChemistryLaTeX
 * Creates a Firefox-compatible version from the Chrome build
 * 
 * IMPORTANT: For Firefox, we rebuild JS files WITHOUT custom name remapping
 * so Mozilla reviewers can read the code. We still minify with esbuild.
 * 
 * Run: node build-firefox.js
 * (Automatically called by build.js after Chrome/dev builds)
 */

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');
const esbuild = require('esbuild');

const SOURCE_DIR = path.dirname(__dirname);  // Parent folder (chem-extension-server)
const CHROME_BUILD_DIR = path.join(path.dirname(SOURCE_DIR), 'ChemistryLaTeX-compiled');
const FIREFOX_OUTPUT_DIR = path.join(path.dirname(SOURCE_DIR), 'ChemistryLaTeX-firefox');
const FIREFOX_ZIP_PATH = path.join(FIREFOX_OUTPUT_DIR, 'ChemistryLaTeX-firefox.zip');

// Firefox extension configuration
const FIREFOX_CONFIG = {
    // Unique extension ID - fresh for new Mozilla submission
    extensionId: '{fddce64f-6901-4661-a45a-c1c43d0eeada}',
    minFirefoxVersion: '140.0',  // Required for data_collection_permissions support
    updateUrl: null
};

console.log('🦊 Building ChemistryLaTeX for Firefox...\n');

/**
 * Copy a file with directory creation
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
 * Get directory size recursively
 */
function getSize(dir) {
    let s = 0;
    for (const f of fs.readdirSync(dir, { withFileTypes: true })) {
        const p = path.join(dir, f.name);
        s += f.isDirectory() ? getSize(p) : fs.statSync(p).size;
    }
    return s;
}

/**
 * Create a ZIP file from the Firefox extension directory
 */
function createZipFile(sourceDir, zipPath) {
    if (fs.existsSync(zipPath)) {
        fs.unlinkSync(zipPath);
    }

    try {
        const tempZipPath = path.join(path.dirname(sourceDir), 'ChemistryLaTeX-firefox-temp.zip');
        if (fs.existsSync(tempZipPath)) {
            fs.unlinkSync(tempZipPath);
        }

        const srcEscaped = sourceDir.replace(/\\/g, '\\\\');
        const destEscaped = tempZipPath.replace(/\\/g, '\\\\');

        execSync(`powershell -Command "Compress-Archive -Path '${srcEscaped}\\*' -DestinationPath '${destEscaped}' -Force"`, {
            stdio: 'pipe',
            encoding: 'utf8',
            timeout: 30000
        });

        if (fs.existsSync(tempZipPath)) {
            fs.copyFileSync(tempZipPath, zipPath);
            fs.unlinkSync(tempZipPath);
            return true;
        }

        return false;
    } catch (e) {
        console.error('   ⚠️ PowerShell zip error:', e.message);
        return false;
    }
}

// Icons are now pre-generated in source assets/ and copied during Chrome build
// Firefox build just copies from Chrome build, so no icon generation needed here

/**
 * Create Firefox-compatible manifest from Chrome manifest
 */
function createFirefoxManifest(chromeManifestPath) {
    const manifest = JSON.parse(fs.readFileSync(chromeManifestPath, 'utf8'));

    // Firefox-specific settings
    manifest.browser_specific_settings = {
        gecko: {
            id: FIREFOX_CONFIG.extensionId,
            strict_min_version: FIREFOX_CONFIG.minFirefoxVersion,
            // Required privacy declaration - object format required by Firefox
            // "none" means extension doesn't collect any user data
            data_collection_permissions: {
                required: ["none"]
            }
        }
    };

    if (FIREFOX_CONFIG.updateUrl) {
        manifest.browser_specific_settings.gecko.update_url = FIREFOX_CONFIG.updateUrl;
    }

    // Convert service_worker to scripts array (Firefox requirement)
    if (manifest.background?.service_worker) {
        const serviceWorkerFile = manifest.background.service_worker;
        manifest.background = {
            scripts: [serviceWorkerFile]
        };
        console.log(`   → Converted service_worker to scripts: [${serviceWorkerFile}]`);
    }

    return manifest;
}

/**
 * CSP Bypass Wrapper for Firefox
 * This code is prepended to content.min.js to intercept fetch() calls
 * and route them through the background script (which bypasses page CSP)
 */
const FIREFOX_CSP_BYPASS_CODE = `
/* Firefox CSP Bypass - Routes fetch through background script */
(function(){
  const _originalFetch = window.fetch;
  const CHEMSERVER = 'server-chemistryrenderer.vercel.app';
  
  window.fetch = async function(url, opts) {
    const urlStr = typeof url === 'string' ? url : url.toString();
    
    // Only intercept requests to our server
    if (urlStr.includes(CHEMSERVER)) {
      console.log('[Firefox CSP Bypass] Routing through background:', urlStr);
      
      return new Promise((resolve, reject) => {
        if (typeof chrome !== 'undefined' && chrome.runtime && chrome.runtime.sendMessage) {
          chrome.runtime.sendMessage(
            { type: 'FETCH_TEXT', url: urlStr, options: opts || {} },
            (response) => {
              if (chrome.runtime.lastError) {
                console.error('[Firefox CSP Bypass] Background error:', chrome.runtime.lastError.message);
                // Try original fetch as fallback
                _originalFetch(url, opts).then(resolve).catch(reject);
                return;
              }
              
              if (response && response.success) {
                // Create a fake Response object
                const headers = new Headers();
                headers.set('content-type', 'image/svg+xml');
                
                resolve(new Response(response.data, {
                  status: 200,
                  statusText: 'OK',
                  headers: headers
                }));
              } else {
                reject(new Error(response?.error || 'Background fetch failed'));
              }
            }
          );
        } else {
          // No chrome.runtime, try original
          _originalFetch(url, opts).then(resolve).catch(reject);
        }
      });
    }
    
    // Pass through for all other URLs
    return _originalFetch(url, opts);
  };
  
  console.log('[ChemistryLaTeX] Firefox CSP bypass installed');
})();
`;

/**
 * Apply Firefox-specific code fixes to JavaScript files
 */
function applyFirefoxCodeFixes(filePath, fileName) {
    if (!fs.existsSync(filePath)) return 0;
    if (!filePath.endsWith('.js')) return 0;

    let content = fs.readFileSync(filePath, 'utf8');
    let changes = 0;

    // ========================================
    // Special handling for content.min.js - Prepend CSP bypass
    // ========================================
    if (fileName === 'content.min.js') {
        // Prepend the Firefox CSP bypass code
        content = FIREFOX_CSP_BYPASS_CODE + '\n' + content;
        changes++;
        console.log('   → Injected Firefox CSP bypass wrapper');
    }

    // ========================================
    // Fix URL filtering for Firefox internal pages
    // ========================================
    const urlPatterns = [
        {
            find: /!tab\.url\.startsWith\('chrome:\/\/'\)\s*&&\s*!tab\.url\.startsWith\('edge:\/\/'\)/g,
            replace: "!tab.url.startsWith('chrome://') && !tab.url.startsWith('edge://') && !tab.url.startsWith('about:') && !tab.url.startsWith('moz-extension://')"
        }
    ];

    for (const pattern of urlPatterns) {
        if (pattern.find.test(content)) {
            content = content.replace(pattern.find, pattern.replace);
            changes++;
        }
    }

    // ========================================
    // Fix extension context detection
    // ========================================
    const contextPatterns = [
        {
            find: /event\.filename\.includes\('chrome-extension:\/\/'\)/g,
            replace: "(event.filename.includes('chrome-extension://') || event.filename.includes('moz-extension://'))"
        }
    ];

    for (const pattern of contextPatterns) {
        if (pattern.find.test(content)) {
            content = content.replace(pattern.find, pattern.replace);
            changes++;
        }
    }

    if (changes > 0) {
        fs.writeFileSync(filePath, content);
    }

    return changes;
}

/**
 * Fix popup.html for Firefox CSP compliance
 * Removes inline event handlers that Firefox blocks
 */
function fixPopupHtmlForFirefox(popupPath) {
    if (!fs.existsSync(popupPath)) return 0;

    let content = fs.readFileSync(popupPath, 'utf8');
    let changes = 0;

    // Remove inline onmouseenter/onmouseleave handlers
    // These are blocked by Firefox's strict CSP
    const inlineHandlerPatterns = [
        // Pattern: onmouseenter="this.style.background='#e74c3c'; this.style.color='white';"
        /\s+onmouseenter="[^"]*"/gi,
        // Pattern: onmouseleave="this.style.background='transparent'; this.style.color='#e74c3c';"
        /\s+onmouseleave="[^"]*"/gi,
        // Pattern: onclick="..." (if any)
        /\s+onclick="[^"]*"/gi,
    ];

    for (const pattern of inlineHandlerPatterns) {
        const matches = content.match(pattern);
        if (matches) {
            content = content.replace(pattern, '');
            changes += matches.length;
        }
    }

    // Add CSS hover styles instead (more reliable and CSP-compliant)
    // Find the closing </style> tag and add our hover styles before it
    const hoverStyles = `
    /* Firefox CSP Fix: CSS-based hover effects */
    #clearCacheBtn:hover {
      background: #e74c3c !important;
      color: white !important;
    }
    #reloadAllBtn:hover {
      background: #2c3e50 !important;
      color: white !important;
    }
  `;

    // Insert before the last </style>
    const lastStyleIndex = content.lastIndexOf('</style>');
    if (lastStyleIndex > -1 && changes > 0) {
        content = content.slice(0, lastStyleIndex) + hoverStyles + content.slice(lastStyleIndex);
    }

    if (changes > 0) {
        fs.writeFileSync(popupPath, content);
    }

    return changes;
}

/**
 * Main Firefox build function
 */
async function buildFirefox() {
    try {
        // Check if Chrome build exists
        if (!fs.existsSync(CHROME_BUILD_DIR)) {
            console.error('❌ Chrome build not found! Run build.js first.');
            console.log('   Expected:', CHROME_BUILD_DIR);
            process.exit(1);
        }

        // Clean previous Firefox build
        if (fs.existsSync(FIREFOX_OUTPUT_DIR)) {
            fs.rmSync(FIREFOX_OUTPUT_DIR, { recursive: true });
        }

        // Copy Chrome build to Firefox directory
        console.log('📦 Copying Chrome build to Firefox directory...');
        copyDir(CHROME_BUILD_DIR, FIREFOX_OUTPUT_DIR);
        console.log('   ✅ Copied all files');

        // Create Firefox-compatible manifest
        console.log('📦 Creating Firefox manifest...');
        const chromeManifestPath = path.join(FIREFOX_OUTPUT_DIR, 'manifest.json');
        const firefoxManifest = createFirefoxManifest(chromeManifestPath);
        fs.writeFileSync(chromeManifestPath, JSON.stringify(firefoxManifest, null, 2));
        console.log('   ✅ manifest.json (Firefox-modified)');
        console.log(`   → Extension ID: ${FIREFOX_CONFIG.extensionId}`);
        console.log(`   → Min Firefox: ${FIREFOX_CONFIG.minFirefoxVersion}`);

        // Icons are pre-generated in source assets/ and copied from Chrome build
        console.log('📦 Verifying icons...');
        const assetsDir = path.join(FIREFOX_OUTPUT_DIR, 'assets');
        const icon128 = path.join(assetsDir, 'icon128.png');
        if (fs.existsSync(icon128)) {
            console.log('   ✅ Icons copied from Chrome build (blue beaker)');
        } else {
            console.log('   ⚠️ Icons missing - run full build.js first');
        }

        // Rebuild JS files from source WITHOUT name remapping (for Mozilla review)
        // This makes the code readable while still being minified
        console.log('📦 Rebuilding JS files (no name remapping for Mozilla review)...');

        // content.min.js - minified but readable variable names
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'content.js')],
            bundle: false,
            minify: true,
            drop: ['console', 'debugger'],
            outfile: path.join(FIREFOX_OUTPUT_DIR, 'content.min.js'),
            target: ['firefox140'],
            legalComments: 'none',
        });
        console.log('   ✅ content.min.js (minified, no name remapping)');

        // popup.min.js
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'popup.js')],
            bundle: false,
            minify: true,
            drop: ['console', 'debugger'],
            outfile: path.join(FIREFOX_OUTPUT_DIR, 'popup.min.js'),
            target: ['firefox140'],
            legalComments: 'none',
        });
        console.log('   ✅ popup.min.js (minified, no name remapping)');

        // background.min.js
        await esbuild.build({
            entryPoints: [path.join(SOURCE_DIR, 'background.js')],
            bundle: false,
            minify: true,
            drop: ['console', 'debugger'],
            outfile: path.join(FIREFOX_OUTPUT_DIR, 'background.min.js'),
            target: ['firefox140'],
            legalComments: 'none',
        });
        console.log('   ✅ background.min.js (minified, no name remapping)');

        // Apply Firefox-specific code fixes
        console.log('📦 Applying Firefox code compatibility fixes...');
        const jsFiles = [
            'content.min.js',
            'background.min.js',
            'popup.min.js'
        ];

        let totalFixes = 0;
        for (const jsFile of jsFiles) {
            const filePath = path.join(FIREFOX_OUTPUT_DIR, jsFile);
            const fixes = applyFirefoxCodeFixes(filePath, jsFile);
            if (fixes > 0) {
                console.log(`   ✅ ${jsFile} (${fixes} fixes applied)`);
                totalFixes += fixes;
            } else {
                console.log(`   ✅ ${jsFile} (no changes needed)`);
            }
        }

        // Fix popup.html inline handlers for Firefox CSP
        console.log('📦 Fixing popup.html for Firefox CSP...');
        const popupHtmlPath = path.join(FIREFOX_OUTPUT_DIR, 'popup.html');
        const popupFixes = fixPopupHtmlForFirefox(popupHtmlPath);
        if (popupFixes > 0) {
            console.log(`   ✅ popup.html (${popupFixes} inline handlers removed)`);
            console.log('   → Added CSS hover styles as replacement');
        } else {
            console.log('   ✅ popup.html (no changes needed)');
        }

        // Create ZIP file for easy loading
        console.log('📦 Creating ZIP file for Firefox...');
        const zipSuccess = createZipFile(FIREFOX_OUTPUT_DIR, FIREFOX_ZIP_PATH);
        if (zipSuccess && fs.existsSync(FIREFOX_ZIP_PATH)) {
            const zipSize = fs.statSync(FIREFOX_ZIP_PATH).size;
            console.log(`   ✅ ChemistryLaTeX-firefox.zip (${(zipSize / 1024).toFixed(1)} KB)`);
        } else {
            console.log('   ⚠️ ZIP creation skipped or failed, but extension files are ready');
            console.log('   → You can still load the extension by selecting manifest.json');
        }

        // Summary
        console.log('\n🦊 Firefox build complete!');
        console.log(`📊 Total size: ${(getSize(FIREFOX_OUTPUT_DIR) / 1024 / 1024).toFixed(2)} MB`);
        console.log(`📁 Output: ${FIREFOX_OUTPUT_DIR}`);
        console.log(`📦 ZIP: ${FIREFOX_ZIP_PATH}`);
        console.log('\n📝 To test in Firefox:');
        console.log('   1. Open Firefox → about:debugging#/runtime/this-firefox');
        console.log('   2. Click "Load Temporary Add-on..."');
        console.log(`   3. Select: ${FIREFOX_ZIP_PATH}`);
        console.log('      OR select manifest.json from the folder');

    } catch (e) {
        console.error('❌ Firefox build failed:', e);
        process.exit(1);
    }
}

// Run the build
buildFirefox();
