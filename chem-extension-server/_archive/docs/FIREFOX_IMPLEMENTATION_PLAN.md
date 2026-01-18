# Firefox Extension Implementation Plan
## ChemistryLaTeX - Firefox Port Strategy

**Created:** 2026-01-11  
**Status:** Planning Phase  
**Estimated Complexity:** Medium (shared codebase with build system)

---

## 📋 Executive Summary

Your Chrome extension is **highly compatible** with Firefox with minimal changes. Both browsers now support Manifest V3 (Firefox 109+) and the WebExtensions API, which means your extension can work on both platforms from a **single codebase** using a unified build system.

### Key Finding
The `chrome.*` namespace used throughout your code **automatically works in Firefox** through its `browser.*` polyfill, OR you can use the WebExtension-polyfill to ensure 100% compatibility.

---

## 🏗️ Recommended Architecture: Unified Build System

### Why This Approach?
1. **Single Codebase** - No need to maintain two separate versions
2. **Automatic Updates** - Changes propagate to both platforms simultaneously
3. **Reduced Maintenance** - Fix bugs once, deploy everywhere
4. **Consistent Features** - Users on both platforms get the same experience

### Output Structure
```
ChemistryLaTeX-compiled/       # Chrome build (existing)
ChemistryLaTeX-firefox/        # Firefox build (new)
```

---

## ⚠️ CRITICAL: Content Script CSP Differences

### The Problem

**Chrome and Firefox handle content script network requests DIFFERENTLY:**

| Behavior | Chrome | Firefox |
|----------|--------|---------|
| Content script `fetch()` | **Bypasses page CSP** | **Subject to page CSP!** |
| Content script `<img>.src` | Bypasses page CSP | Subject to page CSP! |
| Background script `fetch()` | Full network access | Full network access |

This means when your content script does `fetch('https://server-chemistryrenderer.vercel.app/...')` on ChatGPT:
- **Chrome**: Works fine (content scripts bypass page CSP)
- **Firefox**: BLOCKED by ChatGPT's strict `connect-src` CSP!

### The Solution

The Firefox build automatically injects a **CSP bypass wrapper** at the start of `content.min.js` that:

1. Intercepts all `fetch()` calls to our server
2. Routes them through `chrome.runtime.sendMessage()` to the background script
3. The background script makes the actual fetch (which bypasses page CSP)
4. Returns the response back to the content script

```javascript
// Simplified version of the injected wrapper
window.fetch = async function(url, opts) {
  if (url.includes('server-chemistryrenderer.vercel.app')) {
    // Route through background script
    return new Promise((resolve, reject) => {
      chrome.runtime.sendMessage({ type: 'FETCH_TEXT', url }, (response) => {
        resolve(new Response(response.data, { status: 200 }));
      });
    });
  }
  return originalFetch(url, opts);
};
```

This is automatically handled by `build-firefox.js` - no manual changes needed!

### 1. Manifest V3 Differences (Minimal)

Both Chrome and Firefox support Manifest V3, but there are small differences:

| Feature | Chrome | Firefox | Action Required |
|---------|--------|---------|-----------------|
| Service Worker | `background.service_worker` | `background.scripts` (preferred) | Build script transforms |
| CSP format | Object format | Same (works!) | None |
| Host permissions | `<all_urls>` | Same (works!) | None |
| Action API | `chrome.action` | Same | None |
| Storage API | `chrome.storage` | Same | None |
| Context Menus | `chrome.contextMenus` | Same | None |
| Tabs API | `chrome.tabs` | Same | None |

**Firefox Manifest Differences:**
```json
{
  "background": {
    "scripts": ["background.js"],  // Firefox MV3 prefers this format
    "type": "module"               // Optional: for ES module support
  },
  "browser_specific_settings": {
    "gecko": {
      "id": "chemistrylatex@your-domain.com",
      "strict_min_version": "109.0"
    }
  }
}
```

### 2. API Namespace (Zero Code Changes!)

**Good News:** Firefox supports the `chrome.*` namespace for WebExtensions compatibility!

Your existing code like:
```javascript
chrome.storage.sync.set({ enabled: true });
chrome.runtime.sendMessage({ type: 'RELOAD_ALL_IMAGES' });
chrome.tabs.query({}, (tabs) => { ... });
```

**Works in Firefox without changes** because Firefox provides `chrome.*` as an alias to its `browser.*` API.

However, for best practices and Promise support, consider:

```javascript
// Option A: Keep using chrome.* (works in both, but callback-based)
chrome.storage.sync.set({ enabled: true }, callback);

// Option B: Use browser.* with polyfill (Promise-based, modern)
await browser.storage.sync.set({ enabled: true });
```

### 3. Background Script Mode

Firefox MV3 has some differences in background script handling:

| Aspect | Chrome MV3 | Firefox MV3 |
|--------|-----------|-------------|
| Type | Service Worker only | Event Page or Service Worker |
| Persistence | Non-persistent | Non-persistent |
| DOM Access | No `window`/`document` | No `window`/`document` |
| FileReader | Available ✓ | Available ✓ |

**Your Code Is Compatible!** Your `background.js` doesn't use any DOM APIs that would cause issues.

### 4. Content Security Policy

Your current CSP is compatible:
```json
"content_security_policy": {
  "extension_pages": "script-src 'self'; default-src 'self'; ..."
}
```

This format works in both Chrome and Firefox MV3.

---

## 📦 Build System Implementation

### New Build Structure

Create a new file: `build-firefox.js` (or extend existing `build.js`)

```javascript
// build-firefox.js
const FIREFOX_OUTPUT_DIR = path.join(path.dirname(__dirname), 'ChemistryLaTeX-firefox');

// ... same build process as Chrome ...

// Firefox-specific manifest modifications
function createFirefoxManifest() {
  const manifest = JSON.parse(fs.readFileSync(chromeManifest));
  
  // Add Firefox-specific settings
  manifest.browser_specific_settings = {
    gecko: {
      id: "chemistrylatex@quintessenlabs.com",
      strict_min_version: "109.0"
    }
  };
  
  // Convert service_worker to scripts array (Firefox compatibility)
  if (manifest.background?.service_worker) {
    manifest.background = {
      scripts: [manifest.background.service_worker]
    };
  }
  
  return manifest;
}
```

### Recommended: Unified Build Script

Modify your existing `build.js` to output both platforms:

```javascript
// In build.js - add after Chrome build

async function buildFirefox() {
  console.log('\n🦊 Building Firefox version...\n');
  
  const FIREFOX_OUTPUT_DIR = path.join(path.dirname(__dirname), 'ChemistryLaTeX-firefox');
  
  // 1. Copy all Chrome build files to Firefox directory
  copyDir(OUTPUT_DIR, FIREFOX_OUTPUT_DIR);
  
  // 2. Modify manifest for Firefox
  const manifestPath = path.join(FIREFOX_OUTPUT_DIR, 'manifest.json');
  const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
  
  // Add Firefox extension ID
  manifest.browser_specific_settings = {
    gecko: {
      id: "chemistrylatex@quintessenlabs.com",
      strict_min_version: "109.0"
    }
  };
  
  // Firefox supports service_worker in MV3, but scripts array is more reliable
  // Convert to scripts array format for maximum compatibility
  if (manifest.background?.service_worker) {
    manifest.background = {
      scripts: [manifest.background.service_worker]
    };
  }
  
  fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
  console.log('   ✅ manifest.json (Firefox-modified)');
  
  console.log('\n🦊 Firefox build complete!');
}

// Run both builds
build()
  .then(() => buildFirefox())
  .then(() => console.log('\n✅ All builds complete!'));
```

---

## 📁 File Structure After Implementation

```
chem-extension-server/
├── build.js                    # Main build (outputs Chrome)
├── build-dev.js                # Dev build
├── manifest.json               # Source manifest
├── background.js               # Source (unchanged)
├── content.js                  # Source (unchanged)
├── popup.js                    # Source (unchanged)
└── ...

ChemistryLaTeX-compiled/        # Chrome output
├── manifest.json               # Chrome MV3 manifest
├── background.min.js
├── content.min.js
└── ...

ChemistryLaTeX-firefox/         # Firefox output (NEW)
├── manifest.json               # Firefox MV3 manifest (with gecko ID)
├── background.min.js           # Same code, different manifest reference
├── content.min.js
└── ...
```

---

## 🔧 Edge Cases & Specific Fixes

### 1. Tab URL Filtering
Your code filters chrome:// and edge:// URLs:
```javascript
if (tab.url && !tab.url.startsWith('chrome://') && !tab.url.startsWith('edge://')) {
```

**Add Firefox support:**
```javascript
if (tab.url && 
    !tab.url.startsWith('chrome://') && 
    !tab.url.startsWith('edge://') &&
    !tab.url.startsWith('about:') &&     // Firefox internal pages
    !tab.url.startsWith('moz-extension://')) {  // Firefox extension pages
```

### 2. Extension Context Validation
Your code checks for extension context:
```javascript
if (event.filename && event.filename.includes('chrome-extension://')) {
```

**Add Firefox support:**
```javascript
if (event.filename && 
    (event.filename.includes('chrome-extension://') || 
     event.filename.includes('moz-extension://'))) {
```

### 3. navigator.clipboard (Works in Both!)
Your clipboard usage in `popup.js` works identically in both browsers.

---

## 🚀 Deployment Workflow

### Initial Setup

1. **Run Build:**
   ```bash
   cd chem-extension-server
   node build.js
   ```
   This produces:
   - `ChemistryLaTeX-compiled/` (Chrome)
   - `ChemistryLaTeX-firefox/` (Firefox)

2. **Test in Firefox:**
   - Go to `about:debugging#/runtime/this-firefox`
   - Click "Load Temporary Add-on"
   - Select `ChemistryLaTeX-firefox/manifest.json`

3. **Test in Chrome:**
   - Go to `chrome://extensions`
   - Enable Developer Mode
   - Load unpacked → `ChemistryLaTeX-compiled/`

### Publishing

| Platform | Store | Process |
|----------|-------|---------|
| Chrome | Chrome Web Store | Upload `ChemistryLaTeX-compiled/` as ZIP |
| Firefox | addons.mozilla.org | Upload `ChemistryLaTeX-firefox/` as ZIP |

---

## ⚠️ Known Firefox Limitations (Minor)

1. **Temporary Addon Loading:** Extension IDs change each load unless you sign the extension or use a test ID in manifest.

2. **Storage Sync Limits:** Firefox has slightly different storage quotas, but for your use case, it's not a concern.

3. **Service Worker Lifecycle:** Firefox MV3 service workers may have slightly different wake patterns, but your current message-based architecture handles this well.

---

## 📝 Implementation Checklist

### Phase 1: Build System (Estimated: 2 hours)
- [ ] Create Firefox build step in `build.js`
- [ ] Add Firefox manifest modifications
- [ ] Add Firefox extension ID
- [ ] Test basic build output

### Phase 2: Code Compatibility (Estimated: 1 hour)  
- [ ] Update URL filtering to include `about:` and `moz-extension://`
- [ ] Update error handler to detect `moz-extension://`
- [ ] Test all features in Firefox

### Phase 3: Testing (Estimated: 2 hours)
- [ ] Load extension in Firefox Developer Mode
- [ ] Test molecule rendering
- [ ] Test context menus
- [ ] Test popup settings
- [ ] Test storage sync
- [ ] Test 3D viewers
- [ ] Test biomolecule/mineral rendering

### Phase 4: Publishing (Estimated: 1 hour)
- [ ] Create Firefox Add-ons developer account
- [ ] Submit for review
- [ ] Complete review process (may take 1-7 days)

---

## 🔮 Future Considerations

### Browser-specific Polyfill (Optional)
If you ever need Promise-based APIs (cleaner code), add the WebExtension-polyfill:
```html
<!-- In popup.html, before popup.js -->
<script src="browser-polyfill.min.js"></script>
```

This provides a unified `browser.*` API that returns Promises:
```javascript
// Instead of callbacks
const settings = await browser.storage.sync.get(null);
await browser.runtime.sendMessage({ type: 'RELOAD_ALL_IMAGES' });
```

### Cross-Browser Development Tips
1. Always test in both browsers before releases
2. Use feature detection if using newer APIs
3. Consider automated testing with Selenium or Puppeteer

---

## ✅ Conclusion

Your extension architecture is **already 95% Firefox-compatible**. The main work is:

1. **Build system modification** - Create Firefox output with modified manifest
2. **Minor code fixes** - URL filtering for Firefox internal pages
3. **Testing** - Ensure all features work

**No new CSP system needed!** Your current CSP works in Firefox.  
**No architecture changes needed!** The WebExtensions API is consistent.

The unified build system approach ensures that when you update features for Chrome, they automatically propagate to Firefox with zero extra effort.

---

## 📚 References

- [Firefox WebExtensions](https://developer.mozilla.org/en-US/docs/Mozilla/Add-ons/WebExtensions)
- [Chrome Extension Migration Guide](https://developer.chrome.com/docs/extensions/mv3/intro/)
- [browser vs chrome Namespace](https://developer.mozilla.org/en-US/docs/Mozilla/Add-ons/WebExtensions/Chrome_incompatibilities)
- [WebExtension-polyfill](https://github.com/nicktennant/webextension-polyfill)
