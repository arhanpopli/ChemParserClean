# Settings Persistence Fix

## Problem
Settings weren't being saved when toggling options in the popup. When reopening the popup, toggles would revert to their default state.

## Root Cause
When `chrome.storage.sync.get()` is called with a defaults object, it returns the defaults for any keys that don't exist in storage. However, these defaults weren't being saved to storage, so every time the popup opened, it would use the defaults again instead of the user's saved preferences.

## Solution
Added one-time initialization code that:
1. Checks if key settings exist in storage (`sdShowCarbons`, `sdAromaticRings`, `sdShowMethyls`)
2. If they don't exist, saves the defaults to storage
3. This ensures future changes are properly persisted

## Changes Made

### popup.js
1. **Added error handling** to storage operations with `chrome.runtime.lastError` checks
2. **Added console logging** to track when settings are saved and loaded
3. **Added initialization code** to save defaults on first run

### How to Test
1. Open extension popup
2. Toggle "Show Carbon Atoms" ON
3. Close popup
4. Reopen popup
5. **Expected**: "Show Carbon Atoms" should still be ON
6. Check browser console for logs:
   - `[Popup] Setting sdShowCarbons to: true`
   - `[Popup] Successfully saved sdShowCarbons: true`
   - `[Popup] Loaded sdShowCarbons: true`

## If Issues Persist

Run this in the popup console to manually initialize settings:

```javascript
chrome.storage.sync.set({
  sdShowCarbons: true,
  sdAromaticRings: true,
  sdShowMethyls: true,
  sdAtomNumbers: false,
  sdShowHydrogens: false,
  sdFlipHorizontal: false,
  sdFlipVertical: false,
  sdGradientColors: false,
  showTags: false
}, () => {
  console.log('Settings manually initialized');
  location.reload();
});
```

## Debug Commands

### View all stored settings:
```javascript
chrome.storage.sync.get(null, (result) => {
  console.table(result);
});
```

### Clear all settings (reset to defaults):
```javascript
chrome.storage.sync.clear(() => {
  console.log('Storage cleared');
  location.reload();
});
```
