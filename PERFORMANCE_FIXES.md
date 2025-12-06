# Performance & Lazy Loading Fixes

## Changes Made

### 1. Disabled Lazy Loading by Default
**Problem:** Lazy loading was causing extremely slow image loads with requestIdleCallback delays of up to 1 second
**Solution:** Set `performanceMode: false` by default

```javascript
// OLD:
performanceMode: true,
maxVisibleSVGs: 5,

// NEW:
performanceMode: false,  // Disabled by default - lazy loading causes slow image loads
maxVisibleSVGs: 50,  // Increased from 5 to allow more images to load
```

### 2. Removed requestIdleCallback Delays
**Problem:** Images were waiting for browser idle time before loading, causing 1+ second delays
**Solution:** Load immediately when images enter viewport

```javascript
// OLD:
requestIdleCallback(() => {
  loadImage(img);
}, { timeout: 1000 });  // 1 second delay!

// NEW:
loadImage(img);  // Immediate load when in viewport
```

### 3. Increased Concurrent Load Limit
**Problem:** Only 3 images could load simultaneously, causing queue bottleneck
**Solution:** Increased to 10 concurrent loads

```javascript
// OLD:
const maxConcurrentLoads = 3;

// NEW:
const maxConcurrentLoads = 10;
```

### 4. Expanded Viewport Margin
**Problem:** Images only started loading when very close to viewport
**Solution:** Increased rootMargin from 300px to 500px

```javascript
// OLD:
rootMargin: '300px'

// NEW:
rootMargin: '500px'  // Start loading earlier for smoother experience
```

### 5. Added Immediate Loading Mode
**Problem:** Even with performanceMode OFF, lazy loading was still active
**Solution:** Created a pass-through observer that loads immediately

```javascript
if (!settings.performanceMode) {
  // Set up immediate loading - no lazy loading delays
  window.chemfigObserver = {
    observe: (img) => {
      if (img.dataset.src) {
        img.src = img.dataset.src;
        img.dataset.loaded = 'true';
        img.classList.remove('chemfig-loading');
      }
    },
    unobserve: () => {}
  };
  return;  // Skip lazy loading setup entirely
}
```

## Performance Comparison

### Before (with lazy loading):
- ❌ 1+ second delay per image (requestIdleCallback timeout)
- ❌ Maximum 3 concurrent loads
- ❌ Images loaded 300px before viewport
- ❌ Even with performanceMode OFF, lazy loading was active
- ❌ Result: Very slow, choppy loading

### After (immediate loading):
- ✅ No delays - images load immediately
- ✅ Up to 10 concurrent loads
- ✅ Images preload 500px before viewport
- ✅ performanceMode OFF = true immediate loading
- ✅ Result: Fast, smooth loading

## Testing

### Test Immediate Loading (Recommended)
1. Open extension settings
2. Disable "Performance Mode" (if it's enabled)
3. Reload the page
4. Images should load instantly with no delays

### Test Lazy Loading (For slow connections)
1. Open extension settings
2. Enable "Performance Mode"
3. Reload the page
4. Images load as you scroll (but faster than before)

## Debug Commands

Check current performance mode:
```javascript
window.chemRendererDebug.getSettings().performanceMode
```

Toggle performance mode:
```javascript
window.chemRendererDebug.togglePerformanceMode()
```

Force reload all images:
```javascript
window.chemRendererDebug.scanPage()
```

## RCSB Image Loading

The slow RCSB image loading was caused by:
1. Lazy loading delays (now fixed)
2. Search API fetch delays (unavoidable on ChatGPT due to CSP)
3. Concurrent load limits (now increased)

### Expected Behavior Now:
- ✅ Local HTML files: Images load immediately
- ✅ Most websites: Images load quickly
- ⚠️ ChatGPT: Search API blocked by CSP (expected - security restriction)

## Recommendations

### For Fast Loading (Default):
```javascript
performanceMode: false
maxVisibleSVGs: 50
```

### For Memory-Constrained Devices:
```javascript
performanceMode: true
maxVisibleSVGs: 20
```

### For Very Slow Connections:
```javascript
performanceMode: true
maxVisibleSVGs: 10
```

## Files Modified
- `chem-extension/content.js`:
  - Line 1462: `performanceMode: false` (was `true`)
  - Line 1463: `maxVisibleSVGs: 50` (was `5`)
  - Line 2254-2273: Added immediate loading mode when performanceMode is OFF
  - Line 2278: `maxConcurrentLoads: 10` (was `3`)
  - Line 2280-2294: Removed requestIdleCallback delays
  - Line 2306: `rootMargin: '500px'` (was `'300px'`)

## Summary

**The extension now loads images immediately by default, with no artificial delays or restrictions.** Lazy loading can still be enabled via settings for devices with limited memory or slow connections.
