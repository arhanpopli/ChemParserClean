## Extension Debugging Guide

The extension stopped working after the recent changes. Here's how to debug:

### Step 1: Check Browser Console
1. Open ChatGPT (or any page where you use the extension)
2. Open Developer Tools (F12)
3. Go to the Console tab
4. Look for any JavaScript errors (red text)

### Step 2: Check Extension Console
1. Go to `chrome://extensions/`
2. Enable "Developer mode" (top right)
3. Find "Chemistry Formula Renderer"
4. Click "Inspect views: service worker" or "background page"
5. Check for errors in the console

### Step 3: Reload the Extension
1. Go to `chrome://extensions/`
2. Find "Chemistry Formula Renderer"
3. Click the reload icon (circular arrow)
4. Try using the extension again

### Common Issues to Check

#### Issue 1: Syntax Errors
- Check if all functions are properly closed
- Check if all brackets/parentheses match
- Look for missing semicolons or commas

#### Issue 2: Undefined Functions
- The `calculateMoleculeDefaultSize` function was added to `content.js`
- Make sure it's defined before it's used

#### Issue 3: Autocorrection Notice Disabled
- Lines 2724 and 2904 in `content.js` have the autocorrection disabled
- This should NOT cause the extension to break

### What Changed

1. **content.js** (Lines 2724, 2904):
   - Commented out `showAutocorrectNotice()` calls
   
2. **content.js** (Lines 799-828):
   - Added `calculateMoleculeDefaultSize()` function
   - This function calculates parabolic sizing

3. **size-controls.js**:
   - Modified but this file is NOT loaded by the extension
   - Changes here should NOT affect anything

### Quick Fix

If the extension is completely broken, try reverting just the autocorrection changes:

1. Open `content.js`
2. Find lines 2724-2726 and uncomment:
   ```javascript
   if (searchData.wasCorrected) {
     showAutocorrectNotice(img, searchData.originalQuery, searchData.correctedName);
   }
   ```
3. Find lines 2904-2906 and uncomment the same
4. Reload the extension

This will restore the autocorrection notice but should fix any breaking issues.
