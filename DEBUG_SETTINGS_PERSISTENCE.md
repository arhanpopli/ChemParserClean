# Emergency Debugging Guide - Settings Not Persisting

## 🔴 Current Issue
You reported:
- Only "Aromatic Rings" toggle shows as enabled in popup
- When you enable other settings (methyls, atom numbers, etc.) and close the popup, they don't stay enabled
- When you try to disable "Aromatic Rings", it doesn't turn off

## 🧪 Step-by-Step Debugging

### Step 1: Check What's Actually in Storage

1. Open the extension popup
2. Press F12 to open DevTools
3. Go to Console tab
4. Run this command:
```javascript
chrome.storage.sync.get(null, (result) => {
  console.log('=== STORAGE CONTENTS ===');
  console.log(JSON.stringify(result, null, 2));
});
```

**What to look for:**
- If you see only `{"sdAromaticRings": true}` → Storage is corrupted/incomplete
- If you see many settings → Storage is OK, popup loading is broken

---

### Step 2: Test Saving a Setting

1. In the popup console, run:
```javascript
chrome.storage.sync.set({ sdShowMethyls: true }, () => {
  console.log('Saved sdShowMethyls');
  chrome.storage.sync.get(['sdShowMethyls'], (r) => {
    console.log('Verified:', r.sdShowMethyls);
  });
});
```

2. Close and reopen popup
3. Run the Step 1 command again
4. Check if `sdShowMethyls: true` is in storage

**What to look for:**
- If it's saved → Saving works, popup loading is broken
- If it's not saved → Storage API is broken

---

### Step 3: Clear Storage and Start Fresh

**⚠️ WARNING: This will reset ALL extension settings!**

1. In console, run:
```javascript
chrome.storage.sync.clear(() => {
  console.log('Storage cleared');
  chrome.storage.sync.get(null, (r) => console.log('Storage now:', r));
});
```

2. Reload the extension:
   - Go to `chrome://extensions/`
   - Find ChemTex
   - Click reload icon

3. Open popup again
4. Check console logs - should see:
```
[Popup] ========== SETTINGS LOADING DEBUG ==========
[Popup] Raw storage (storedSettings): {}
[Popup] Final merged settings: {sdShowCarbons: true, sdAromaticRings: true, ...}
```

---

### Step 4: Test Toggle Behavior

1. Open popup
2. Watch console closely
3. Click "Show Methyls" toggle **ON**
4. Should see:
```
[Popup] Saving sdShowMethyls: true
[Popup] Verified sdShowMethyls saved as: true
```

5. **DON'T CLOSE POPUP YET**
6. In console, run:
```javascript
chrome.storage.sync.get(['sdShowMethyls'], (r) => console.log('Check:', r));
```
7. Should show `{sdShowMethyls: true}`

8. **NOW close popup**
9. **Reopen popup**
10. Watch console logs:
```
[Popup] Raw storage (storedSettings): {sdShowMethyls: true, sdAromaticRings: true, ...}
[Popup] Set sdShowMethylsToggle.checked to: true | Actual value now: true
```

11. **VERIFY**: Is the "Show Methyls" toggle visually ON?

---

### Step 5: Check for Event Listener Issues

If toggles save but don't show as enabled when popup reopens, the issue might be CSS or HTML:

1. Open popup
2. Inspect the toggle element:
   - Right-click on "Show Methyls" toggle
   - Click "Inspect"
3. In Elements tab, look for the `<input type="checkbox">` element
4. Check its attributes:
   - Should have `checked` attribute if enabled
   - Should NOT have `checked` attribute if disabled

5. In console, run:
```javascript
const toggle = document.getElementById('sdShowMethylsToggle');
console.log('Element:', toggle);
console.log('Checked property:', toggle.checked);
console.log('Checked attribute:', toggle.getAttribute('checked'));
console.log('Disabled?:', toggle.disabled);
```

---

## 🔍 Expected Console Output (Normal Behavior)

### When Opening Popup (First Time):
```
[Popup] ========== SETTINGS LOADING DEBUG ==========
[Popup] Raw storage (storedSettings): {}
[Popup] Defaults object: {sdShowCarbons: true, sdAromaticRings: true, sdShowMethyls: true, ...}
[Popup] Final merged settings: {sdShowCarbons: true, sdAromaticRings: true, sdShowMethyls: true, ...}
[Popup] SmilesDrawer settings specifically: {
  sdShowCarbons: true,
  sdAromaticRings: true,
  sdShowMethyls: true,
  sdAtomNumbers: false,
  sdShowHydrogens: false,
  sdFlipHorizontal: false,
  sdFlipVertical: false
}
[Popup] ================================================
[Popup] ========== TOGGLE ELEMENT CHECK ==========
[Popup] sdShowCarbonsToggle exists? true
[Popup] sdAromaticRingsToggle exists? true
[Popup] sdShowMethylsToggle exists? true
...
[Popup] ================================================
[Popup] Set sdShowCarbonsToggle.checked to: true | Actual value now: true
[Popup] Set sdAromaticRingsToggle.checked to: true | Actual value now: true
[Popup] Set sdShowMethylsToggle.checked to: true | Actual value now: true
```

### When Clicking Toggle OFF:
```
[Popup] Saving sdShowMethyls: false
[Popup] Verified sdShowMethyls saved as: false
[Popup] Successfully saved and verified sdShowMethyls: false
```

### When Reopening Popup After Toggle OFF:
```
[Popup] Raw storage (storedSettings): {sdShowMethyls: false, sdShowCarbons: true, sdAromaticRings: true, ...}
[Popup] Final merged settings: {sdShowMethyls: false, ...}
[Popup] Set sdShowMethylsToggle.checked to: false | Actual value now: false
```

---

## 🐛 Common Issues & Fixes

### Issue 1: "Raw storage shows correct values but toggles don't reflect them"
**Cause:** HTML elements not found or CSS hiding them
**Fix:**
1. Check if toggle elements exist (Step 5)
2. Check if toggles are disabled: `toggle.disabled = false;`
3. Check CSS for `display: none` or `visibility: hidden`

### Issue 2: "Settings save but disappear when popup reopens"
**Cause:** Chrome sync quota exceeded (100KB limit)
**Fix:**
1. Check storage size:
```javascript
chrome.storage.sync.getBytesInUse(null, (bytes) => {
  console.log('Storage used:', bytes, 'bytes');
});
```
2. If > 100KB, clear old data

### Issue 3: "Only sdAromaticRings persists, others don't"
**Cause:** Old code saved only that one setting
**Fix:**
1. Clear storage (Step 3)
2. Reload extension
3. Test each toggle individually

### Issue 4: "Toggle shows as ON but console says it's OFF"
**Cause:** Event listener attached multiple times
**Fix:**
1. Check for duplicate event listeners
2. Reload extension
3. Open popup fresh

---

## 🎯 Quick Fix (Nuclear Option)

If nothing works:

1. **Backup important data** (if any)

2. **Uninstall extension:**
   - `chrome://extensions/`
   - Click "Remove" on ChemTex

3. **Clear ALL extension data:**
```javascript
chrome.storage.sync.clear();
chrome.storage.local.clear();
```

4. **Reload VS Code** (or wherever you're developing)

5. **Reinstall extension:**
   - Load unpacked from folder

6. **Test fresh:**
   - Open popup
   - All toggles should show defaults (ON for carbons/rings/methyls, OFF for others)
   - Toggle one setting
   - Close and reopen
   - Verify it persisted

---

## 📊 Share These Console Logs

Please copy and share:

1. **When you first open popup:**
```
Look for: [Popup] ========== SETTINGS LOADING DEBUG ==========
```

2. **When you click a toggle:**
```
Look for: [Popup] Saving <settingName>: <value>
Look for: [Popup] Verified <settingName> saved as: <value>
```

3. **When you reopen popup:**
```
Look for: [Popup] Raw storage (storedSettings): {...}
Look for: [Popup] Set <toggle>Toggle.checked to: <value> | Actual value now: <value>
```

4. **Storage contents:**
```javascript
chrome.storage.sync.get(null, (r) => console.log(JSON.stringify(r, null, 2)));
```

---

**Debug Date:** December 11, 2025
**Files:** `popup.js` (with enhanced logging), `content.js`
