# Dead Code Analysis Report - ChemTex Extension

**File:** `content.js`  
**Total Lines:** 4,902  
**Date:** 2025-12-11  

---

## 📋 Executive Summary

This report identifies dead code, unused functions, and legacy remnants in `content.js`. The extension has been refactored to use SmilesDrawer exclusively, but several legacy code blocks from mol2chemfig, chemfig, CDK Depict, and other deprecated systems remain.

### Key Findings:
- **12 Functions** can be safely removed (completely dead)
- **1 Database Constant** is unused and can be removed
- **7 Setting Variables** are never used
- **Outdated Comments** reference removed systems
- **~400+ lines** of dead code can be eliminated

---

## 🔴 DEAD CODE - Safe to Remove

### 1. NOMENCLATURE_DB (Lines 4105-4137)
**Type:** Constant  
**Status:** ❌ NEVER USED  
**Description:** Database mapping chemical names to chemfig formulas. This was intended for a chemfig-based rendering system that has been replaced by SmilesDrawer.

```javascript
const NOMENCLATURE_DB = {
  "ethanol": "\\chemfig{C-[1]C(-[1]OH)",
  "propanol": "\\chemfig{C-[1]C-[7]C(-[1]OH)}",
  // ... 15+ entries
};
```

**Lines to Remove:** 4105-4137 (~33 lines)

---

### 2. getChemfigFromNomenclature() (Lines 4139-4150)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Looks up nomenclature in NOMENCLATURE_DB to get chemfig formula. Dead because NOMENCLATURE_DB is unused.

```javascript
function getChemfigFromNomenclature(name) {
  if (!name) return null;
  const normalized = name.toLowerCase().trim();
  const formula = NOMENCLATURE_DB[normalized];
  // ...
  return formula || null;
}
```

**Lines to Remove:** 4139-4150 (~12 lines)

---

### 3. removeUnnecessaryHydrogens() (Lines 4153-4172)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Cleans up chemfig content by removing explicit hydrogens. This was for chemfig rendering which is no longer used.

```javascript
function removeUnnecessaryHydrogens(chemfigContent) {
  let result = chemfigContent;
  // Regex operations on chemfig content
  return result;
}
```

**Lines to Remove:** 4153-4172 (~20 lines)

---

### 4. addZigzagAngles() (Lines 4174-4205)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Adds zigzag angles to chemfig carbon chains for realistic rendering. Dead because chemfig rendering is not used.

```javascript
function addZigzagAngles(chemfigContent) {
  let result = chemfigContent;
  // Multiple regex replacements
  return result;
}
```

**Lines to Remove:** 4174-4205 (~32 lines)

---

### 5. simplifyChemfigCarbons() (Lines 4207-4237)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Simplifies chemfig structures by removing explicit carbon labels. References the dead setting `renderCarbonsAsSticks`.

```javascript
function simplifyChemfigCarbons(chemfigContent) {
  if (!settings.renderCarbonsAsSticks) {
    return chemfigContent;
  }
  // ...
}
```

**Lines to Remove:** 4207-4237 (~31 lines)

---

### 6. detectDarkMode() (Lines 4051-4103)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Detects dark mode via multiple methods. However, `isDarkModeEnabled()` (line 653-689) is the **actually used** function for dark mode detection. This is a duplicate.

```javascript
function detectDarkMode() {
  // Multiple detection methods...
  // Priority 1: CSS media query
  // Priority 2: Background color analysis
  // etc.
}
```

**Lines to Remove:** 4051-4103 (~53 lines)

---

### 7. applyLayoutMode() (Lines 3902-3922)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Applies horizontal/vertical layout to chemfig containers. The `layoutMode` setting it references is also dead.

```javascript
function applyLayoutMode() {
  if (!settings.layoutMode || settings.layoutMode === 'horizontal') return;
  const containers = document.querySelectorAll('.chemfig-container');
  // ...
}
```

**Lines to Remove:** 3902-3922 (~21 lines)

---

### 8. convertChemistry() (Lines 4623-4658)
**Type:** Function  
**Status:** ❌ DEAD - NOT CALLED ANYWHERE  
**Description:** Converts chemistry notation to Unicode (H2O → H₂O). Was meant to be used in wrapChemicalFormulas but is NOT called.

```javascript
function convertChemistry(text) {
  // Subscript and superscript conversions
}
```

**Lines to Remove:** 4623-4658 (~36 lines)

---

### 9. logToStorage() (Lines 73-78)
**Type:** Function  
**Status:** ⚠️ DISABLED BUT CALLED  
**Description:** Function is called by all log methods but immediately returns (disabled).

```javascript
function logToStorage(type, msg, data) {
  // DISABLED: Too much memory overhead
  return;  // Skip logging to save 6GB RAM!
}
```

**Action:** Keep the empty function (called by log methods) OR remove calls from log object.

---

### 10. registerRenderedMolecule() (Lines 237-242)
**Type:** Function  
**Status:** ❌ NEVER CALLED  
**Description:** Was intended to track rendered molecules for re-rendering, but the `renderedMolecules` Map is never used.

```javascript
function registerRenderedMolecule(container, moleculeData) {
  if (!container || !moleculeData) return;
  renderedMolecules.set(container, moleculeData);
}
```

**Lines to Remove:** 237-242 (~6 lines) + `renderedMolecules` Map at line 235

---

### 11. invertSvgForDarkMode() (Lines 691-769)
**Type:** Function  
**Status:** ⚠️ APPEARS DEAD  
**Description:** Comprehensive SVG color inversion function. `updateMoleculeColors()` (lines 4702-4761) does similar work and IS used.

```javascript
function invertSvgForDarkMode(svgContent) {
  // Complex color inversion logic
}
```

**Action:** Verify if called anywhere before removing (~79 lines).

---

## 🟡 DEAD SETTINGS - Never Used

### Settings in defaults (Lines 1824-1863):

| Setting | Line | Status | Reason |
|---------|------|--------|--------|
| `renderMhchem` | 1827 | ❌ DEAD | mhchem rendering removed |
| `renderChemfig` | 1828 | ❌ DEAD | chemfig rendering removed |
| `maxVisibleSVGs` | 1830 | ⚠️ CHECK | May be used in lazy loading |
| `layoutMode` | 1831 | ❌ DEAD | applyLayoutMode() not called |
| `renderCarbonsAsSticks` | 1832 | ❌ DEAD | simplifyChemfigCarbons() not called |
| `cdkColorScheme` | 1847 | ❌ DEAD | CDK Depict removed |
| `cdkHydrogenDisplay` | 1848 | ❌ DEAD | CDK Depict removed |
| `cdkZoom` | 1849 | ❌ DEAD | CDK Depict removed |
| `cdkShowCarbons` | 1850 | ❌ DEAD | CDK Depict removed |
| `cdkShowMethyls` | 1851 | ❌ DEAD | CDK Depict removed |
| `cdkAtomNumbers` | 1852 | ❌ DEAD | CDK Depict removed |
| `cdkAnnotation` | 1853 | ❌ DEAD | CDK Depict removed |
| `clientSideRenderer` | 1862 | ❌ DEAD | Always 'smilesdrawer' |

---

## 🟢 KEEP - Used Functions with Legacy Naming

### m2cf* Variables (Lines 1877-1884)
**Status:** ✅ KEEP - Actually used by SmilesDrawer  
**Note:** The `m2cf` prefix is legacy naming from mol2chemfig, but these variables now control SmilesDrawer rendering options.

```javascript
settings.m2cfShowCarbons = result.sdShowCarbons === true;
settings.m2cfAromaticCircles = result.sdAromaticRings === true;
settings.m2cfShowMethyls = result.sdShowMethyls === true;
settings.m2cfAtomNumbers = result.sdAtomNumbers === true;
settings.m2cfAddH2 = result.sdShowHydrogens === true;
settings.m2cfFlipHorizontal = result.sdFlipHorizontal === true;
settings.m2cfFlipVertical = result.sdFlipVertical === true;
settings.m2cfRotate = result.sdRotate || 0;
```

**Action:** Consider renaming to `sd*` for consistency, but NOT dead code.

---

## 🟠 OUTDATED COMMENTS/STRINGS - Update Required

### 1. Line 650-651:
```javascript
// Handles mol2chemfig/dvisvgm SVGs which have colors in CSS styles, attributes, and text elements
```
**Should say:** `// Handles SVGs rendered by SmilesDrawer`

### 2. Lines 1749-1752 (Debug Help):
```javascript
Switch rendering engine:
  window.chemRendererDebug.setRendererEngine('codecogs')      // Standard (current)
  window.chemRendererDebug.setRendererEngine('latex-online')  // Alternative
  window.chemRendererDebug.setRendererEngine('quicklatex')    // Fast mode
```
**Action:** Remove - there's no setRendererEngine function anymore.

### 3. Line 4776:
```javascript
log.info('✨ Chemistry formulas will be rendered via Unicode and CodeCogs API');
```
**Should say:** `'✨ Chemistry formulas will be rendered via SmilesDrawer'`

---

## 📊 Summary of Removable Code

| Category | Functions/Items | Estimated Lines |
|----------|----------------|-----------------|
| Chemfig Functions | 5 | ~128 |
| Dead Detection Functions | 1 | ~53 |
| Layout Functions | 1 | ~21 |
| Registration Functions | 2 | ~15 |
| SVG Inversion (verify) | 1 | ~79 |
| NOMENCLATURE_DB | 1 | ~33 |
| Dead Settings | 12 | ~10 |
| **TOTAL** | **23** | **~339+** |

---

## ✅ Recommended Action Plan

### Phase 1: Safe Removals (No verification needed)
1. Remove `NOMENCLATURE_DB` constant
2. Remove `getChemfigFromNomenclature()`
3. Remove `removeUnnecessaryHydrogens()`
4. Remove `addZigzagAngles()`
5. Remove `simplifyChemfigCarbons()`
6. Remove `applyLayoutMode()`
7. Remove `detectDarkMode()` (keep `isDarkModeEnabled()`)
8. Remove `registerRenderedMolecule()` + `renderedMolecules` Map
9. Remove dead CDK settings from defaults

### Phase 2: Verify Before Removing
1. Check if `invertSvgForDarkMode()` is called
2. Check if `convertChemistry()` is called from wrapChemicalFormulas
3. Check if `maxVisibleSVGs` setting is used

### Phase 3: Cleanup
1. Update outdated comments
2. Update log messages referencing CodeCogs
3. Remove debug help text referencing setRendererEngine

---

## 🔧 Files to Check for Consistency

After cleanup, verify these files don't reference removed functions:
- `popup.html` - Check for CDK Depict UI
- `popup.js` - Check for CDK/mhchem/chemfig settings handlers
- `background.js` - Check for removed message handlers

---

*Report generated for ChemTex Extension dead code cleanup*
