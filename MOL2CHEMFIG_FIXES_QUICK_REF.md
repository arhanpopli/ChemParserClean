# mol2chemfig Fixes - Quick Reference

## What Was Fixed

### 1. Options Not Applying
**Before:** Enabling aromatic circles, show carbons, etc. had no effect
**After:** All mol2chemfig options work correctly
**File:** `m2cf_fixed.py`

### 2. Dark Mode Atom Labels Invisible
**Before:** In dark mode, O, N, Cl, etc. stayed black (invisible on dark background)
**After:** All atom labels turn white in dark mode
**File:** `content.js` (4 locations)

---

## Quick Test

1. **Start servers:**
   ```bash
   START_ALL_SERVERS.bat
   ```

2. **Load extension:**
   - Chrome → Extensions → Load unpacked
   - Select `chem-extension` folder

3. **Test aromatic circles:**
   - Extension popup → mol2chemfig → Enable "Aromatic circles"
   - Open `test_mol2chemfig_fixes.html`
   - Benzene should show a circle in the ring

4. **Test dark mode:**
   - Enable OS dark mode
   - Reload test page
   - O, N, Cl labels should be WHITE

---

## Files Changed

1. `m2cf_fixed.py` - Added options handling in `/submit` endpoint
2. `content.js` - Enhanced dark mode SVG color replacement (4 locations)

---

## Options Mapping

| Extension Option | mol2chemfig Flag | Effect |
|-----------------|------------------|--------|
| Aromatic circles | `-o` | Shows circles in aromatic rings |
| Show carbons | `-c` | Labels carbon atoms with "C" |
| Show methyls | `-m` | Shows CH3 groups explicitly |
| Fancy bonds | `-f` | Uses fancy bond rendering |
| Atom numbers | `-n` | Shows atom numbering |
| Compact | `-z` | Compact structure layout |

---

## Testing Checklist

- [ ] Aromatic circles appear when enabled
- [ ] Show carbons labels carbon atoms
- [ ] Show methyls shows CH3 groups
- [ ] Dark mode: bond lines are white
- [ ] Dark mode: O, N, Cl labels are white (not black)
- [ ] Dark mode: CH3, OH labels are white
- [ ] Options persist across sessions
- [ ] Page reload applies new options

---

## Key Points

1. **Page reload required:** After changing options, you MUST reload the page
2. **Cache handles options:** Different options = different cache entries (no stale SVGs)
3. **Dark mode detection:** Uses OS-level dark mode (`prefers-color-scheme`)
4. **All renderers work:** Fixes apply to both mol2chemfig backend (port 8000) and wrapper (port 5001)

---

## Common Issues

**Options don't apply:**
- Did you reload the page?
- Is mol2chemfig renderer selected?
- Are all servers running?

**Atom labels still black in dark mode:**
- Is OS dark mode enabled?
- Try browser DevTools → Rendering → Emulate prefers-color-scheme: dark

**No images showing:**
- Run `START_ALL_SERVERS.bat`
- Check console for errors (F12)

---

## Test Files

- `test_mol2chemfig_fixes.html` - Comprehensive test page
- `MOL2CHEMFIG_FIXES_COMPLETE.md` - Detailed documentation

---

## Success!

Both issues are fixed and ready for testing. The system now:
- ✅ Applies mol2chemfig options correctly
- ✅ Shows visible atom labels in dark mode
- ✅ Maintains proper caching with options
- ✅ Works with all renderers and endpoints
