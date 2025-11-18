# 🧪 Extension Testing Checklist

## Pre-Test Setup

- [ ] MoleculeViewer server running: `python run_server.py`
- [ ] Extension loaded unpacked in Chrome
- [ ] All files updated: `content.js`, `popup.html`, `popup.js`

---

## Part 1: UI & Settings

### Popup Display
- [ ] Click extension icon → popup opens
- [ ] Rendering Engine dropdown visible
- [ ] "🧪 MoleculeViewer (Best)" option available
- [ ] MoleculeViewer Options section **hidden** by default

### Selecting MoleculeViewer
- [ ] Click Rendering Engine dropdown
- [ ] Select "🧪 MoleculeViewer (Best)"
- [ ] MoleculeViewer Options section **now visible**
- [ ] Success message shows: "Switched to 🧪 MoleculeViewer Server..."
- [ ] Section shows all 8 controls:
  - [ ] Show Carbon Atoms (toggle)
  - [ ] Show Methyl Groups (toggle)
  - [ ] Aromatic Circles (toggle)
  - [ ] Fancy Bonds (toggle)
  - [ ] Atom Numbers (toggle)
  - [ ] Flip Horizontal (toggle)
  - [ ] Flip Vertical (toggle)
  - [ ] Hydrogen Display (dropdown)

### Settings Persistence
- [ ] Toggle "Show Carbon Atoms" ON
- [ ] Success message: "Show carbons enabled..."
- [ ] Close popup completely
- [ ] Reopen popup
- [ ] "Show Carbon Atoms" is still ON
- [ ] Settings persisted! ✓

---

## Part 2: Rendering

### Default Rendering
- [ ] Go to ChatGPT or webpage with chemistry
- [ ] Type or find: `\chemfig{C-C-C}` (propane)
- [ ] Wait 2-3 seconds
- [ ] SVG structure appears
- [ ] F12 → Console: search for `🔬 Using MoleculeViewer`
- [ ] Found! ✓

### Check SMILES Conversion
- [ ] F12 → Console
- [ ] Look for: `Converted chemfig → SMILES: CCC`
- [ ] Found! ✓

### Check Rendering Options Sent
- [ ] F12 → Console
- [ ] Look for: `Rendering options:`
- [ ] Should show all options being sent:
  ```
  {
    show_carbons: false,
    show_methyls: false,
    aromatic_circles: true,
    fancy_bonds: true,
    atom_numbers: false,
    hydrogens: 'keep',
    flip_horizontal: false,
    flip_vertical: false
  }
  ```
- [ ] Found! ✓

---

## Part 3: Option Testing

### Test Show Carbon Atoms
- [ ] Popup: Toggle "Show Carbon Atoms" ON
- [ ] Success message: "Show carbons enabled. Reload page to apply."
- [ ] Reload the webpage (F5)
- [ ] Structure now shows C labels
- [ ] Expected: "C-C-C" visible for propane
- [ ] Verified! ✓

### Test Show Methyl Groups
- [ ] Popup: Toggle "Show Methyl Groups" ON
- [ ] Toggle "Show Carbon Atoms" OFF
- [ ] Reload page
- [ ] Structure now shows CH₃ labels
- [ ] Verified! ✓

### Test Aromatic Circles
- [ ] Use structure: `\chemfig{*6(=(-)-=(-)-=(-)-)}` (benzene)
- [ ] Toggle "Aromatic Circles" ON
- [ ] Reload page
- [ ] Benzene ring shows circle inside
- [ ] Toggle OFF → reload → circle gone
- [ ] Verified! ✓

### Test Flip Options
- [ ] Popup: Toggle "Flip Horizontal" ON
- [ ] Reload page
- [ ] Structure is flipped left-right
- [ ] Toggle OFF → reload → back to normal
- [ ] Toggle "Flip Vertical" ON
- [ ] Reload page
- [ ] Structure is flipped up-down
- [ ] Verified! ✓

### Test Hydrogen Display
- [ ] Popup: "Hydrogen Display" dropdown
- [ ] Try each option:
  - [ ] "keep" - hydrogen as drawn
  - [ ] "add" - all hydrogens shown
  - [ ] "delete" - no hydrogens shown
- [ ] Reload page each time
- [ ] Structure changes appropriately
- [ ] Verified! ✓

### Test Atom Numbers
- [ ] Popup: Toggle "Atom Numbers" ON
- [ ] Reload page
- [ ] Structure shows numbered atoms (1, 2, 3...)
- [ ] Toggle OFF → reload → numbers gone
- [ ] Verified! ✓

---

## Part 4: Complex Structures

### Test Multiple Structures
- [ ] Use multiple chemfig on same page
- [ ] All should render from server
- [ ] Check console: multiple `✅ MoleculeViewer SVG loaded` messages
- [ ] Verified! ✓

### Test Lazy Loading
- [ ] Enable performance mode in popup (should be on by default)
- [ ] Load page with 10+ structures
- [ ] Scroll down slowly
- [ ] Structures load only when visible
- [ ] Console shows: `Loading SVG (#1)`, `Loading SVG (#2)`, etc.
- [ ] Max 3 concurrent loads (check console)
- [ ] Verified! ✓

### Test with Different Molecules
- [ ] `\chemfig{C-C-OH}` (ethanol)
- [ ] `\chemfig{C(=O)C}` (acetone)
- [ ] `\chemfig{*6(=(-NH2)-=(-)-=(-)-)}` (aniline)
- [ ] All should render with server
- [ ] Verified! ✓

---

## Part 5: Error Handling

### Test Server Unavailable
- [ ] Stop server: Close terminal running `python run_server.py`
- [ ] Reload webpage with chemistry formula
- [ ] Should fall back gracefully
- [ ] Check console: `❌ MoleculeViewer fetch failed`
- [ ] Fallback to CodeCogs (URL based rendering)
- [ ] Formula still displays (just via CodeCogs)
- [ ] Verified! ✓

### Test Invalid SMILES
- [ ] (Shouldn't happen with proper conversion)
- [ ] If it does: Server returns error as SVG text
- [ ] Page still functional, not broken
- [ ] Verified! ✓

---

## Part 6: Switching Engines

### Switch Back to CodeCogs
- [ ] Popup: Rendering Engine → select "CodeCogs"
- [ ] MoleculeViewer Options **hidden**
- [ ] Go to webpage with chemistry
- [ ] Reload page
- [ ] Structure renders using CodeCogs (check URL in console)
- [ ] Works! ✓

### Switch Back to MoleculeViewer
- [ ] Popup: Rendering Engine → select "🧪 MoleculeViewer"
- [ ] MoleculeViewer Options **visible** again
- [ ] Go to webpage with chemistry
- [ ] Reload page
- [ ] Structure renders using server
- [ ] `🔬 Using MoleculeViewer server` in console
- [ ] Works! ✓

---

## Part 7: Performance

### Test Load Time
- [ ] Page with 5 structures
- [ ] Measure time from page load to all structures visible
- [ ] Should be reasonable (< 10 seconds)
- [ ] Structures load progressively (lazy-loading)
- [ ] Verified! ✓

### Test Memory Usage
- [ ] Open DevTools → Memory tab
- [ ] Load page with 20+ structures
- [ ] Check memory doesn't spike excessively
- [ ] Tab remains responsive
- [ ] Verified! ✓

---

## Part 8: Bonus Tests

### Test Dark Mode
- [ ] Browser dark mode ON
- [ ] Load structure
- [ ] SVG renders with light colors on dark background
- [ ] Verified! ✓

### Test Mobile View
- [ ] DevTools → Toggle device toolbar
- [ ] Select mobile device (iPhone, etc.)
- [ ] Load chemistry formula
- [ ] Structure responsive on mobile
- [ ] Options still accessible in popup
- [ ] Verified! ✓

### Test Storage Migration
- [ ] Delete chrome storage: DevTools → Application → Storage → Clear
- [ ] Reload extension popup
- [ ] All settings should use defaults
- [ ] Toggle options again, all should work
- [ ] Verified! ✓

---

## Summary

### All Tests Passed? ✅
- [ ] UI tests passed
- [ ] Rendering tests passed
- [ ] Option tests passed
- [ ] Complex structure tests passed
- [ ] Error handling tests passed
- [ ] Engine switching tests passed
- [ ] Performance tests passed
- [ ] Bonus tests passed

**Status: READY FOR PRODUCTION** 🎉

### Issues Found?
Document them here:
```
Issue #1: [Description]
- Severity: [Low/Medium/High]
- Steps to reproduce: [Steps]
- Expected vs Actual: [Expected] vs [Actual]

Issue #2: ...
```

---

## Performance Metrics to Track

| Metric | Target | Actual |
|--------|--------|--------|
| Page load time | < 5s | ? |
| Avg structure render | < 1s | ? |
| Max concurrent loads | 3 | ? |
| Lazy-loading delay | < 300ms | ? |
| Memory per structure | < 500KB | ? |
| Console errors | 0 | ? |

---

**Last Updated:** [Current Date]
**Tested By:** [Your Name]
**Status:** ✅ READY
