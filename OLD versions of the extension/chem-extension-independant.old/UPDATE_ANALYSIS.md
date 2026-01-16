# ChemTex Independent Version - Update Analysis

## ⚠️ CRITICAL: EXTENSIVE DIVERGENCE DETECTED

After detailed code review, the two versions have **significantly diverged** over the past month.
This is NOT a simple update - it requires substantial refactoring.

## Architecture Comparison

### Server Version (`chem-extension-server/`)
- **Rendering**: Server-side via ChemistryLaTeX Server (Vercel) using mol2chemfig
- **SVG Generation**: Server generates SVGs with marker-based theming
- **chem: Parsing**: Strict `type=value` only (`chem:mol=benzene:`)
- **Features**: IUPAC support, direct ID lookups, stricter flag parsing
- **File Sizes**: content.js = 231KB, popup.js = 70KB, popup.html = 53KB

### Independent Version (`chem-extension-independant/`)
- **Rendering**: Client-side via SmilesDrawer (bundled JS)
- **SVG Generation**: SmilesDrawer draws directly to canvas/SVG
- **chem: Parsing**: Uses integrated-search.js for compound lookup
- **Features**: Older feature set, still supports `chem:name:` format
- **File Sizes**: content.js = 295KB, popup.js = 59KB, popup.html = 46KB

## Key Missing Features in Independent Version

### 1. Stricter chem: Tag Parsing ❌
**Server has:** Only accepts `chem:type=value:` format
**Independent has:** Still accepts generic `chem:name:` format

### 2. IUPAC Support ❌
**Server has:** `chem:iupac=2-methylpropan-1-ol:` using OPSIN
**Independent has:** No IUPAC support

### 3. Direct ID Lookups ❌
**Server has:** `chem:pbdid=4RHV:`, `chem:cid=2244:`, `chem:codid=1234567:`
**Independent has:** Limited direct ID support in integrated-search.js

### 4. Marker-Based Theming ❌
**Server has:** Server renders gray marker SVGs, client applies theme colors
**Independent has:** SmilesDrawer has built-in themes (different approach)

### 5. Flag Parsing Fixes ❌
**Server has:** Fixed legacy flag parsing that misread hyphens as flags
**Independent has:** May still have the old buggy flag parsing

### 6. Updated Popup UI ❌
**Server has:** Updated popup.html (53KB) with new settings
**Independent has:** Older popup.html (46KB)

### 7. Updated USAGE.md ❌
**Server has:** Strict `type=value` documentation
**Independent has:** Old documentation with `chem:aspirin:` example

---

## RECOMMENDED APPROACH

### Option A: Full Sync (RECOMMENDED - 6+ hours)
1. Copy entire content.js structure from server version
2. Replace server API calls with SmilesDrawer rendering
3. Keep integrated-search.js for compound lookups
4. Sync popup.html and popup.js

### Option B: Selective Updates (2-3 hours)
1. Update only the chem: parsing logic
2. Update USAGE.md
3. Add IUPAC and direct ID support
4. Skip theming changes (SmilesDrawer has its own)

### Option C: Minimal Critical Fixes (1 hour)
1. Fix flag parsing bug (hyphens misread as flags)
2. Update USAGE.md
3. Update popup Discord link

---

## Immediate Actions (Starting Now)

### Phase 1: Critical Bug Fixes
1. ✅ Update USAGE.md with new Discord link and stricter syntax
2. Fix the flag parsing in content.js or integrated-search.js

### Phase 2: Feature Parity
1. Add IUPAC support
2. Add direct ID lookups (pbdid, cid, codid)
3. Update popup.html/popup.js

---

## Files to Modify
1. `content.js` - Major updates (parsing, flag handling)
2. `integrated-search.js` - Add IUPAC support
3. `popup.html` - UI sync
4. `popup.js` - Settings sync
5. `background.js` - Minor sync
6. `USAGE.md` - Documentation update
7. `manifest.json` - Version bump

## Estimated Total Effort: 6-8 hours for full sync
