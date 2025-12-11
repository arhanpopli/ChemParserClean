# Storage Debug Instructions

## Test if settings are saving:

1. Open the extension popup
2. Open browser console (F12)
3. Toggle "Show Carbon Atoms"
4. Check console for:
   - `[Popup] Setting sdShowCarbons to: true`
   - `[Popup] Successfully saved sdShowCarbons: true`

5. Close and reopen the popup
6. Check console for:
   - `[Popup] Loading settings: {...}`
   - `[Popup] Loaded sdShowCarbons: true`

## If settings show as `undefined`:

The defaults in popup.js don't match what's in storage. Run this in the popup console:

```javascript
chrome.storage.sync.get(null, (result) => {
  console.log('All storage:', result);
});
```

## To manually fix storage:

```javascript
chrome.storage.sync.set({
  sdShowCarbons: true,
  sdAromaticRings: true,
  sdShowMethyls: true
}, () => {
  console.log('Settings manually set');
  location.reload();
});
```

## To clear all storage and start fresh:

```javascript
chrome.storage.sync.clear(() => {
  console.log('Storage cleared');
  location.reload();
});
```
