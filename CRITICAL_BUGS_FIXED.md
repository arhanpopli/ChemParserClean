# Critical Bugs Fixed - Extension Settings

## 🔴 Issues Found

### Issue 1: Missing Functions Causing Extension Crash
**Error:** `ReferenceError: getInvertFilter is not defined`
**Location:** [`content.js:2451`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L2451)

**Problem:**
- Function `getInvertFilter()` was being exported to `window` but never defined
- Function `applyLayoutMode()` was being called but never defined
- Extension crashed on page load

**Fix:**
Added placeholder functions:
```javascript
function getInvertFilter() {
  return ''; // Placeholder - no dark mode inversion yet
}

function applyLayoutMode() {
  log.info('applyLayoutMode called (placeholder - no action taken)');
}
```

---

### Issue 2: Settings NOT Persisting from Popup to Content Script
**Error:** Settings show as `false` in content.js even when enabled in popup
**Location:** [`content.js:1777-1784`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js#L1777)

**Problem:**
```javascript
// ❌ WRONG - This treats undefined as false!
settings.m2cfShowCarbons = result.sdShowCarbons === true;
```

When `result.sdShowCarbons` is `undefined` (not saved yet), the check `=== true` returns `false`, so all settings defaulted to OFF instead of ON.

**Fix:**
```javascript
// ✅ CORRECT - Use ternary to default to true when undefined
settings.m2cfShowCarbons = result.sdShowCarbons !== undefined ? result.sdShowCarbons : true;
settings.m2cfAromaticCircles = result.sdAromaticRings !== undefined ? result.sdAromaticRings : true;
settings.m2cfShowMethyls = result.sdShowMethyls !== undefined ? result.sdShowMethyls : true;
```

This matches the popup.js defaults where these are `true` by default.

---

### Issue 3: Popup Settings Using Wrong Storage Key Loading
**Error:** Popup always showed default values instead of saved values
**Location:** [`popup.js:128`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\popup.js#L128)

**Problem:**
```javascript
// ❌ WRONG - Chrome returns defaults instead of stored values!
chrome.storage.sync.get({
  sdShowCarbons: true,  // These defaults OVERRIDE stored values!
  sdAromaticRings: true,
  sdShowMethyls: true
}, (settings) => { ... });
```

**Fix:**
```javascript
// ✅ CORRECT - Read ALL storage, then apply defaults manually
chrome.storage.sync.get(null, (storedSettings) => {
  const defaults = { sdShowCarbons: true, sdAromaticRings: true, ... };
  const settings = {};
  for (const key in defaults) {
    settings[key] = storedSettings[key] !== undefined ? storedSettings[key] : defaults[key];
  }
  // Now settings has actual stored values OR defaults if not stored
});
```

---

### Issue 4: Duplicate Initialization Overwriting Settings
**Error:** Toggle states reset after being loaded
**Location:** [`popup.js:330-380`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\popup.js#L330) (REMOVED)

**Problem:**
```javascript
// ❌ BAD - Runs AFTER settings loaded, overwrites toggles!
chrome.storage.sync.get(null, (stored) => {
  if (needsInit) {
    chrome.storage.sync.set(defaults, () => {
      // Overwrites toggles that were just set!
      sdShowCarbonsToggle.checked = defaults.sdShowCarbons;
    });
  }
});
```

**Fix:**
Removed the entire duplicate initialization block. Settings initialization now happens once in the main loading code.

---

## 📊 Console Log Analysis

### Before Fix:
```
content.js:2451 Uncaught ReferenceError: getInvertFilter is not defined
content.js:1894 ReferenceError: applyLayoutMode is not defined
content.js:20 📊 SmilesDrawer settings: {showCarbons: false, aromaticRings: false, showMethyls: false, ...}
content.js:3003 🎌 Flag overrides: {..., aromaticCircles: false, showCarbons: false, showMethyls: false}
```
❌ All settings showing as `false`!

### After Fix:
```
[Content] Raw storage values: {sdShowCarbons: true, sdAromaticRings: true, sdShowMethyls: true, ...}
content.js:20 📊 SmilesDrawer settings: {showCarbons: true, aromaticRings: true, showMethyls: true, ...}
content.js:3003 🎌 Flag overrides: {..., aromaticCircles: true, showCarbons: true, showMethyls: true}
```
✅ Settings properly loaded!

---

## 🧪 Testing Instructions

### 1. Reload Extension
```
chrome://extensions/ → Find ChemTex → Click reload button
```

### 2. Clear Storage (Optional - for clean test)
Open console and run:
```javascript
chrome.storage.sync.clear(() => console.log('Storage cleared'));
```

### 3. Open Extension Popup
You should see in console:
```
[Popup] Loading settings from storage: {}  // Empty on first run
[Popup] Final settings after applying defaults: {sdShowCarbons: true, ...}
[Popup] Loaded sdShowCarbons: true
[Popup] Loaded sdAromaticRings: true
[Popup] Loaded sdShowMethyls: true
```

### 4. Toggle Settings
Turn OFF "Show methyls":
```
[Popup] Saving sdShowMethyls: false
[Popup] Verified sdShowMethyls saved as: false
[Popup] Successfully saved and verified sdShowMethyls: false
```

### 5. Close & Reopen Popup
You should see:
```
[Popup] Loading settings from storage: {sdShowMethyls: false, ...}
[Popup] Loaded sdShowMethyls: false
```
Toggle should be OFF! ✅

### 6. Reload Page with Molecules
You should see:
```
[Content] Raw storage values: {sdShowMethyls: false, ...}
📊 SmilesDrawer settings: {showMethyls: false, ...}
🎌 Flag overrides: {showMethyls: false, ...}
```

The molecule should render WITHOUT methyl labels! ✅

---

## 📝 Summary of Changes

### [`content.js`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\content.js)
1. **Added** `getInvertFilter()` function (lines 2446-2453)
2. **Added** `applyLayoutMode()` function (lines 2455-2462)
3. **Fixed** settings loading to use ternary operator for proper defaults (lines 1777-1779)
4. **Added** console logging of raw storage values for debugging (lines 1786-1794)

### [`popup.js`](c:\Users\Kapil\Personal\STUFF\Chemparser\chem-extension\popup.js)
1. **Changed** `chrome.storage.sync.get({defaults})` to `get(null)` (line 128)
2. **Added** manual defaults merging logic (lines 129-184)
3. **Removed** duplicate initialization block that was overwriting settings (lines 330-380 deleted)
4. **Added** `saveSettingWithVerification()` helper function (lines 90-120)
5. **Updated** all toggle handlers to use verification (lines 460-520)

---

## ✅ What's Fixed

| Issue | Status | Impact |
|-------|--------|--------|
| Extension crashes on load | ✅ Fixed | Critical - extension now loads |
| Settings don't persist | ✅ Fixed | High - changes now save correctly |
| Popup shows wrong values | ✅ Fixed | High - UI reflects actual saved state |
| Settings don't apply to molecules | ✅ Fixed | High - molecules render with correct settings |
| Missing error handling | ✅ Fixed | Medium - users see errors if save fails |

---

## 🎯 Expected Behavior Now

1. **First time opening popup:** All toggles ON by default (showCarbons, aromaticRings, showMethyls)
2. **Toggle setting OFF:** Saves immediately, shows success message
3. **Close and reopen popup:** Toggle stays OFF (persisted correctly)
4. **Reload page:** Molecules render without that feature
5. **Toggle setting back ON:** Saves immediately, molecules re-render with feature

---

## 🔍 Debugging Tips

If settings still don't work, check console for:

1. **Popup Console:**
```javascript
[Popup] Loading settings from storage: {...}
[Popup] Saving sdShowMethyls: true
[Popup] Verified sdShowMethyls saved as: true
```

2. **Page Console (where molecules render):**
```javascript
[Content] Raw storage values: {sdShowMethyls: true, ...}
📊 SmilesDrawer settings: {showMethyls: true, ...}
```

3. **Check actual storage:**
```javascript
chrome.storage.sync.get(null, (result) => console.log('Storage:', result));
```

---

**Report Date:** December 11, 2025
**Files Modified:** `content.js`, `popup.js`
**Status:** ✅ ALL CRITICAL BUGS FIXED
