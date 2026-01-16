# Content.js Function Audit Report

**Total Lines:** 6,643  
**Total Functions/Sections:** ~85 functions + code blocks  
**Date:** 2026-01-02

---

## SECTION 1: CORE SERVER COMMUNICATION (Lines 1-87)
**Status: ✅ KEEP - Essential**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `buildServerSvgUrl()` | 28-75 | Builds URL to ChemTex server for SVG rendering | ✅ KEEP - Core |
| `fetchFromChemTexServer()` | 77-87 | Fetches JSON data from server (pdbid/codid) | ✅ KEEP - Core |

---

## SECTION 2: THEMING SYSTEM (Lines 88-181)
**Status: ✅ KEEP - Used for SVG recoloring**

| Item | Lines | Purpose | Keep/Remove |
|------|-------|---------|-------------|
| `MARKER_COLORS` | 93-112 | Gray placeholder colors in server SVGs | ✅ KEEP |
| `THEME_COLORS` | 113-158 | Actual colors for 5 themes (light/dark/matrix/github/carbon) | ✅ KEEP |
| `applyThemeColors()` | 160-181 | Replace marker colors with theme colors | ✅ KEEP - Core |

---

## SECTION 3: ERROR HANDLING & LOGGING (Lines 182-280)
**Status: ⚠️ REVIEW - Potentially verbose**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| Global error handlers | 182-218 | Catch crashes silently | ✅ KEEP |
| `log.info/success/warn/error/debug/inject` | 220-267 | Console logging with prefixes | ⚠️ SIMPLIFY - Verbose |
| `logToStorage()` | 275-280 | Stores logs for debug | ⚠️ CONSIDER REMOVE - Dev only |

---

## SECTION 4: CACHING SYSTEM (Lines 281-438)
**Status: ✅ KEEP - Essential for performance**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `generateImageCacheKey()` | 299-314 | Creates cache key from SMILES+options | ✅ KEEP |
| `getCachedImage()` | 316-322 | Get cached rendered SVG | ✅ KEEP |
| `setCachedImage()` | 324-337 | Store rendered SVG in cache | ✅ KEEP |
| `loadSmilesCache()` | 350-366 | Load SMILES cache from storage | ✅ KEEP |
| `saveSmilesCache()` | 368-384 | Save SMILES cache to storage | ✅ KEEP |
| `getCachedSmiles()` | 385-416 | Get cached SMILES for name | ✅ KEEP |
| `setCachedSmiles()` | 418-438 | Store SMILES in cache | ✅ KEEP |

---

## SECTION 5: DEDUPLICATION (Lines 439-481)
**Status: ✅ KEEP - Prevents duplicate API calls**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `pendingSearches` Map | 439-447 | Track in-flight searches | ✅ KEEP |
| `getOrCreatePendingSearch()` | 448-481 | Deduplicate concurrent searches | ✅ KEEP |

---

## SECTION 6: CACHE CLEARING (Lines 483-516)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `clearAllCaches()` | 483-503 | Clear all caches + bust timestamp | ✅ KEEP |
| `stripStereochemistry()` | 505-516 | Remove @/@@/\/ from SMILES | ✅ KEEP |

---

## SECTION 7: INSTANT SETTINGS APPLICATION (Lines 517-1018)
**Status: ✅ KEEP - Core re-rendering**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `getAffectedImageTypes()` | 529-597 | Determine which images need re-render | ✅ KEEP |
| `reRenderAllMolecules()` | 599-608 | Trigger re-render (delegates) | ✅ KEEP |
| `addLoadingIndicator()` | 621-676 | Show spinner on image | ✅ KEEP |
| `removeLoadingIndicator()` | 678-694 | Remove spinner | ✅ KEEP |
| `showPendingIndicator()` | 696-700 | Show update indicator | ⚠️ REVIEW - Unused? |
| `hidePendingIndicator()` | 702-705 | Hide update indicator | ⚠️ REVIEW - Unused? |
| `applyAverageSizeScaling()` | 707-756 | CSS scale transform for size slider | ✅ KEEP |
| `lazyReRenderMolecules()` | 758-1018 | Main lazy re-render logic | ✅ KEEP - Core |

---

## SECTION 8: LAZY RE-RENDER OBSERVER (Lines 1019-1154)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `setupLazyReRenderObserver()` | 1032-1062 | IntersectionObserver for re-render | ✅ KEEP |
| `processImageIfNeeded()` | 1064-1123 | Process single image re-render | ✅ KEEP |
| `triggerVisibleImagesReRender()` | 1125-1154 | Force re-render visible images | ✅ KEEP |

---

## SECTION 9: SINGLE MOLECULE RE-RENDER (Lines 1156-1261)
**Status: ✅ KEEP - Recently refactored**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `reRenderSingleMolecule()` | 1156-1229 | Re-render one molecule via server | ✅ KEEP - Core |
| `reloadAllImages()` | 1232-1261 | Force reload all images | ✅ KEEP |

---

## SECTION 10: CSS-ONLY SETTINGS (Lines 1263-1305)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `applyCssOnlySettings()` | 1263-1305 | Apply CSS changes without re-fetch | ✅ KEEP |

---

## SECTION 11: MESSAGE HANDLERS (Lines 1306-1379)
**Status: ✅ KEEP**

| Handler | Lines | Purpose | Keep/Remove |
|---------|-------|---------|-------------|
| `APPLY_SETTINGS` | 1314-1343 | Handle popup settings | ✅ KEEP |
| `RELOAD_ALL_IMAGES` | 1346-1356 | Handle reload command | ✅ KEEP |
| `CLEAR_CACHE` | 1358-1374 | Handle cache clear | ✅ KEEP |

---

## SECTION 12: DARK MODE DETECTION (Lines 1380-1522)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `isDarkModeEnabled()` | 1387-1431 | Detect page dark mode | ✅ KEEP |
| `getCSSFilterForTheme()` | 1433-1489 | Get CSS filter for theme | ✅ KEEP |
| `getEffectiveTheme()` | 1491-1522 | Resolve auto-adapt theme | ✅ KEEP |

---

## SECTION 13: BACKGROUND FETCH (Lines 1523-1606)
**Status: ⚠️ PARTIAL - backgroundFetchJSON removed**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `isExtensionContextValid()` | 1532-1543 | Check if extension valid | ✅ KEEP |
| `backgroundFetchBlob()` | 1545-1585 | Fetch blob via background | ✅ KEEP - Used for RCSB images |
| `directFetchBlob()` | 1587-1606 | Direct blob fetch fallback | ✅ KEEP |

---

## SECTION 14: SMILES BRIDGE (Lines 1607-1693)
**Status: ✅ KEEP - Recently refactored**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `smilesBridge()` | 1615-1663 | Convert name→SMILES via server | ✅ KEEP - Core |
| `getPubChemCID()` | 1668-1693 | Get CID via smilesBridge | ✅ KEEP |

---

## SECTION 15: PERFORMANCE MONITORING (Lines 1694-1755)
**Status: ⚠️ REVIEW - Dev/debug only**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `chemRendererPerformance` object | 1701-1755 | Track load times, counts | ⚠️ OPTIONAL - Debug only |
| `recordLoad()` | 1714-1723 | Record load duration | ⚠️ OPTIONAL |
| `recordStructure()` | 1725-1727 | Count structures | ⚠️ OPTIONAL |
| `recordFormula()` | 1729-1731 | Count formulas | ⚠️ OPTIONAL |
| `getStats()` | 1733-1745 | Get performance stats | ⚠️ OPTIONAL |
| `logStats()` | 1747-1755 | Log stats to console | ⚠️ OPTIONAL |

---

## SECTION 16: SIZE CONTROLS (Lines 1756-2143)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `getImageKey()` | 1773-1778 | Get storage key for image | ✅ KEEP |
| `getPageImageKey()` | 1780-1784 | Get page-specific key | ✅ KEEP |
| `loadImageSize()` | 1786-1820 | Load saved size | ✅ KEEP |
| `saveImageSize()` | 1822-1842 | Save image size | ✅ KEEP |
| `createSizeControls()` | 1844-1911 | Create +/- buttons | ✅ KEEP |
| `adjustImageSize()` | 1913-1940 | Handle size change | ✅ KEEP |
| `wrapImageWithSizeControls()` | 1942-2083 | Wrap image in container | ✅ KEEP |
| `applyScaleToImage()` | 2085-2143 | Apply CSS scale | ✅ KEEP |

---

## SECTION 17: DEBUG INTERFACE (Lines 2144-2431)
**Status: ⚠️ REVIEW - Many verbose functions**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `getLogs()` | 2147-2150 | Get log history | ⚠️ DEBUG |
| `getSettings()` | 2151-2154 | Show settings | ⚠️ DEBUG |
| `scanPage()` | 2155-2158 | Manual scan | ⚠️ DEBUG |
| `testFormulas()` | 2159-2172 | Test formula parsing | ❌ REMOVE - Debug only |
| `getCurrentFormulas()` | 2173-2192 | List detected formulas | ❌ REMOVE - Debug only |
| `rotateFormulas()` | 2193-2205 | Rotate all images | ⚠️ OPTIONAL |
| `getRotationHelp()` | 2206-2227 | Help text | ❌ REMOVE - Verbose |
| `togglePerformanceMode()` | 2228-2238 | Toggle perf mode | ⚠️ DEBUG |
| `setMaxVisibleSVGs()` | 2239-2243 | Set max SVGs | ⚠️ DEBUG |
| `getPerformanceStats()` | 2244-2254 | Get perf stats | ⚠️ DEBUG |
| `setLayoutMode()` | 2255-2265 | Set layout | ⚠️ DEBUG |
| `toggleLayoutMode()` | 2266-2272 | Toggle layout | ⚠️ DEBUG |
| `getLayoutSettings()` | 2273-2282 | Get layout | ⚠️ DEBUG |
| `clearMemory()` | 2283-2313 | Clear caches | ⚠️ DEBUG |
| `checkMemory()` | 2314-2329 | Memory usage | ⚠️ DEBUG |
| `toggleCarbonLabels()` | 2330-2341 | Toggle carbons | ⚠️ DEBUG |
| `setSizePreset()` | 2342-2353 | Size preset | ⚠️ DEBUG |
| `getRenderingHelp()` | 2355-2382 | Help text | ❌ REMOVE - Verbose |
| `clearCache()` | 2388-2392 | Clear cache | ⚠️ DEBUG |
| `viewSmilesCache()` | 2394-2403 | View cache | ⚠️ DEBUG |
| `viewImageCache()` | 2405-2412 | View cache | ⚠️ DEBUG |
| `getCacheStats()` | 2414-2429 | Cache stats | ⚠️ DEBUG |

---

## SECTION 18: SETTINGS INITIALIZATION (Lines 2432-2660)
**Status: ✅ KEEP**

| Item | Lines | Purpose | Keep/Remove |
|------|-------|---------|-------------|
| `settings` object | 2445-2491 | Default settings | ✅ KEEP |
| Storage listener | 2492-2582 | Load/merge settings | ✅ KEEP |
| Message handlers | 2614-2660 | Handle popup messages | ✅ KEEP |

---

## SECTION 19: INITIALIZATION (Lines 2661-2863)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `initializeRenderer()` | 2662-2695 | Initialize extension | ✅ KEEP |
| `injectStyles()` | 2697-2863 | Inject CSS styles | ✅ KEEP |

---

## SECTION 20: FLAG PARSING (Lines 2864-3332)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `parseChemFlags()` | ~2864+ | Parse chem:: flags | ✅ KEEP - Core |
| `stripFlagsFromName()` | 3255-3267 | Remove flags from name | ✅ KEEP |
| `applyFlagOverrides()` | 3269-3300 | Apply flag overrides | ✅ KEEP |
| `wrapImageWithRotationContainer()` | 3302-3332 | Rotation wrapper | ✅ KEEP |

---

## SECTION 21: HOVER CONTROLS (Lines 3334-3803)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `addHoverControls()` | 3334-3783 | Add molecule name + 3D button | ✅ KEEP - Core |
| `getInvertFilter()` | 3785-3793 | Dark mode filter | ✅ KEEP |
| `applyLayoutMode()` | 3795-3803 | Apply layout mode | ✅ KEEP |

---

## SECTION 22: URL BUILDERS (Lines 3805-3913)
**Status: ✅ KEEP - Recently refactored**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `buildRCSBImageUrl()` | 3805-3817 | Build server-proxied RCSB URL | ✅ KEEP |
| `buildMolViewEmbedUrl()` | 3819-3905 | Build MolView iframe URL | ✅ KEEP |

---

## SECTION 23: LAZY LOADING CORE (Lines 3914-5144)
**Status: ✅ KEEP - Main rendering engine**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `setupLazyLoading()` | 3914-5144 | Main lazy loading setup | ✅ KEEP - Core |
| `.decrementActiveLoads()` | 3935-3943 | Track active loads | ✅ KEEP |
| `.loadImage()` | 3978-4001 | Load single image | ✅ KEEP |
| `.renderBiomolecule2D()` | 4003-4106 | Render biomolecule | ✅ KEEP |
| `.renderClientSide()` | 4108-4332 | Main render function | ✅ KEEP - Core |
| `.show3DViewerInline()` | 4335-4807 | Show 3D MolView | ✅ KEEP |
| `.applySharpenFilter()` | 4812-4897 | Sharpen images | ⚠️ REVIEW - Used? |
| `.autoCropCanvas()` | 4899-4952 | Auto-crop | ⚠️ REVIEW - Used? |
| `.loadMoleculeImage()` | 4955-4976 | Load molecule | ✅ KEEP |
| `.processLoadQueue()` | 4982-5003 | Process queue | ✅ KEEP |
| `._queueMoleculeLoad()` | 5063-5075 | Queue image | ✅ KEEP |
| `.preloadNearbyMolecules()` | 5091-5131 | Preload nearby | ⚠️ OPTIONAL |

---

## SECTION 24: SCAN AND RENDER (Lines 5145-5375)
**Status: ✅ KEEP - Core**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `scanAndRender()` | 5152-5174 | Debounced page scan | ✅ KEEP |
| `scanAndRenderImmediate()` | 5176-5375 | Main page scanning | ✅ KEEP - Core |

---

## SECTION 25: FORMULA WRAPPING (Lines 5378-5781)
**Status: ✅ KEEP - Core**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `wrapChemicalFormulas()` | 5378-5781 | Parse chem:: tags | ✅ KEEP - Core |

---

## SECTION 26: MUTATION OBSERVER (Lines 5784-5946)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `observePageChanges()` | 5784-5920 | Watch DOM for changes | ✅ KEEP |
| `cleanupObservers()` | 5922-5946 | Cleanup on unload | ✅ KEEP |

---

## SECTION 27: DARK MODE SUPPORT (Lines 5947-6018)
**Status: ✅ KEEP**

| Function | Lines | Purpose | Keep/Remove |
|----------|-------|---------|-------------|
| `updateMoleculeColors()` | 5959-6018 | Update colors on theme change | ✅ KEEP |

---

## SECTION 28: CONTEXT MENU HANDLERS (Lines 6046-6643)
**Status: ✅ KEEP**

| Handler | Lines | Purpose | Keep/Remove |
|---------|-------|---------|-------------|
| `INSPECT_MOLECULE` | 6049-6120 | Right-click render as molecule | ✅ KEEP |
| `RENDER_SMILES` | 6125-6202 | Right-click render as SMILES | ✅ KEEP |
| `RENDER_BIOMOLECULE` | 6203-6273 | Right-click render as biomolecule | ✅ KEEP |
| `RENDER_MINERAL` | 6278-6348 | Right-click render as mineral | ✅ KEEP |
| `RERENDER_IMAGE` | 6353-6482 | Re-render existing image | ✅ KEEP |
| `EDIT_FLAGS` | 6486-6639 | Flag editor dialog | ✅ KEEP |

---

# SUMMARY

## Functions to REMOVE (verbose/unused):
1. `testFormulas()` - Debug only
2. `getCurrentFormulas()` - Debug only  
3. `getRotationHelp()` - Verbose help text
4. `getRenderingHelp()` - Verbose help text
5. `showPendingIndicator()` - Appears unused
6. `hidePendingIndicator()` - Appears unused

## Functions to REVIEW (debug/optional):
- Performance monitoring functions (~40 lines)
- Debug interface functions (~150 lines)
- `applySharpenFilter()` - May be unused
- `autoCropCanvas()` - May be unused

## Estimated removable lines: ~200-250 lines

## Core functions to KEEP:
- All caching functions
- All rendering functions  
- All lazy loading functions
- All message handlers
- All theme/dark mode functions
- All size control functions
- All context menu handlers
