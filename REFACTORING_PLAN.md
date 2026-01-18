# Content.js Refactoring Plan - Code Reduction Analysis

**Current Stats:**
- Total Lines: 4,403
- Total Functions: 133
- Target: Reduce by ~20-30% (~880-1,320 lines)

---

## 🔴 CRITICAL: Dead/Broken Code to Remove

### 1. Layout Functions Reference Non-Existent `applyLayoutMode()`
**Lines:** 1547-1574
- [`setLayoutMode`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1547-L1557) - Calls `applyLayoutMode()` which doesn't exist
- [`toggleLayoutMode`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1558-L1564) - Also calls non-existent function
- [`getLayoutSettings`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1565-L1574) - References `settings.layoutMode` which was removed

**Action:** Remove all 3 functions (28 lines)
**Reason:** `layoutMode` setting was removed, `applyLayoutMode()` doesn't exist

---

## 🟡 MAJOR: Consolidation Opportunities

### 2. Duplicate Performance Stats Functions
**Issue:** Two `getPerformanceStats` in debug interface
- First: [Line 1444-1446](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1444-L1446) - Simple wrapper
- Second: [Line 1536-1546](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1536-L1546) - Detailed implementation

**Action:** Keep second one (detailed), remove first (3 lines)

### 3. Monster Function: `setupLazyLoading` (1,336 lines!)
**Lines:** 2453-3788
**Issue:** Contains 11 nested helper functions. This is ~30% of entire file!

**Nested Functions Inside:**
- `loadImage` (24 lines) - [2493-2516](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L2493-L2516)
- `renderBiomolecule2D` (104 lines) - [2518-2621](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L2518-L2621)
- `renderClientSide` (605 lines!) - [2623-3227](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L2623-L3227)
- `show3DViewerInline` (318 lines) - [3229-3546](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L3229-L3546)
- `applySharpenFilter` (86 lines) - [3551-3636](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L3551-L3636)
- `autoCropCanvas` (54 lines) - [3638-3691](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L3638-L3691)
- `loadMoleculeImage` (22 lines) - [3694-3715](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L3694-L3715)

**Action:** Extract to module-level functions
**Estimated Savings:** ~50 lines (from reduced nesting/scope management)

### 4. Redundant Cache Functions
**Issue:** Debug interface has duplicate cache functions

In `chemRendererDebug`:
- [`clearCache`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1680-L1684) - Calls `clearAllCaches()`
- [`viewSmilesCache`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1686-L1695) - Wrapper around direct access
- [`viewImageCache`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1697-L1704) - Another wrapper
- [`getCacheStats`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1706-L1721) - Stats wrapper

**Solution:** These can reference existing functions directly instead of wrapping
**Savings:** ~30 lines

---

## 🟢 MEDIUM: Simplification Targets

### 5. Excessive Debug Interface Functions
**Lines:** 1433-1724 (292 lines)

Many functions are trivial wrappers:
- [`getLogs`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1436-L1439) - Returns `logHistory`
- [`getSettings`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1440-L1443) - Returns `settings`
- [`scanPage`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1447-L1450) - Calls `scanAndRender()`
- [`rotateFormulas`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1485-L1497) - Simple DOM manipulation
- [`togglePerformanceMode`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1520-L1530) - 11 lines
- [`setMaxVisibleSVGs`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1531-L1535) - 5 lines

**Action:** Simplify to property accessors where possible
**Savings:** ~40 lines

### 6. Redundant Help Functions
- [`getRotationHelp`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1498-L1519) - 22 lines of console output
- [`getRenderingHelp`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1647-L1674) - 28 lines of console output

**Action:** Consolidate into single help function or remove (user can read docs)
**Savings:** ~30 lines

### 7. Size Control Constants
**Lines:** 1053-1069 (17 lines)

Multiple size-related constants that overlap:
```javascript
const SIZE_STEP = 20;
const MIN_SIZE = 100;
const MAX_SIZE = 800;
const DEFAULT_WIDTH = 400;
const DEFAULT_HEIGHT = 350;
const BASE_SIZE = 150;
const SIZE_SCALE_FACTOR = 0.5;
const MAX_DEFAULT_SIZE = 400;
```

**Action:** Consolidate into single `SIZE_CONFIG` object
**Savings:** ~5 lines

### 8. Repetitive Size Control Functions
- [`getImageKey`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1070-L1075) - 6 lines
- [`getPageImageKey`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1077-L1081) - 5 lines

These can be combined into single function with optional pageUrl parameter
**Savings:** ~3 lines

---

## 🔵 MINOR: Code Quality Improvements

### 9. Performance Monitoring Redundancy
**Lines:** 991-1052

The `chemRendererPerformance` object tracks metrics that are rarely used:
- `recordLoad` - Tracks every image load time
- `loadTimes` array - Limited to 50 but still memory overhead

**Action:** Simplify to basic counters only
**Savings:** ~15 lines

### 10. Legacy Test Cases
**Lines:** 1451-1464

[`testFormulas`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1451-L1464) function with hardcoded test strings including chemfig examples

**Action:** Remove or update with SmilesDrawer examples
**Savings:** ~14 lines

### 11. Unused Memory Functions
- [`clearMemory`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1575-L1605) - 31 lines, aggressive cache clearing
- [`checkMemory`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1606-L1621) - 16 lines, uses deprecated `performance.memory`

**Action:** Simplify or remove (modern browsers handle memory well)
**Savings:** ~30 lines

### 12. Redundant Carbon Toggle
[`toggleCarbonLabels`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1622-L1633) - 12 lines

This duplicates functionality in settings system
**Action:** Remove (use popup settings instead)
**Savings:** ~12 lines

### 13. Size Preset Function
[`setSizePreset`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1634-L1645) - 12 lines

Minimal functionality, setting is auto-applied
**Action:** Remove from debug interface
**Savings:** ~12 lines

---

## 📊 Estimated Total Savings

| Category | Lines Saved |
|----------|-------------|
| Dead layout functions | 28 |
| Duplicate performance stats | 3 |
| setupLazyLoading refactor | 50 |
| Cache function wrappers | 30 |
| Debug interface simplification | 40 |
| Help functions consolidation | 30 |
| Size constants consolidation | 5 |
| Performance monitoring simplify | 15 |
| Test functions cleanup | 14 |
| Memory functions removal | 30 |
| Carbon toggle removal | 12 |
| Size preset removal | 12 |
| **TOTAL** | **269 lines** |

---

## 🎯 Implementation Priority

### Phase 1: Remove Dead Code (31 lines)
1. ✅ Remove `setLayoutMode`, `toggleLayoutMode`, `getLayoutSettings`
2. ✅ Remove duplicate `getPerformanceStats` (first instance)

### Phase 2: Extract setupLazyLoading (~50 lines)
1. Move nested functions to module level
2. Keep only IntersectionObserver setup in `setupLazyLoading`

### Phase 3: Simplify Debug Interface (~134 lines)
1. Remove wrapper functions
2. Consolidate help functions
3. Remove redundant toggles/setters

### Phase 4: Code Quality (~54 lines)
1. Consolidate constants
2. Simplify performance monitoring
3. Remove memory functions
4. Update test cases

**Total Estimated Reduction: 269+ lines (~6%)**

---

## ⚠️ Notes for Implementation Agent

1. **DO NOT REMOVE:**
   - Log statements (`log.info`, `log.debug`, etc.)
   - Comments (they help readability)
   - Core functionality (only remove wrappers/duplicates)

2. **VERIFY BEFORE REMOVING:**
   - Check if function is called anywhere else with grep
   - Ensure no external references (popup.js, background.js)

3. **TEST AFTER EACH PHASE:**
   - Load extension and test molecule rendering
   - Check console for errors
   - Verify settings still work

4. **EXTRACTED FUNCTIONS:**
   When extracting from `setupLazyLoading`, maintain closure access to `settings` and `log` objects

---

**Report Generated:** 2025-12-11
**For:** ChemTex v6.0 Code Reduction Project
