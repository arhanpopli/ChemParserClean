# ChemParser Smart Caching & Instant Settings Plan

## Problem Summary
1. **Settings not persisting** - When extension reloads, settings reset to defaults
2. **Page reload required** - Changing settings requires full page reload
3. **No caching** - Same molecules re-fetched and re-rendered on every page

## Solution Design

### 1. Fix Settings Persistence
**Issue:** `chrome.storage.sync.get` uses defaults if key missing, but extension reload may clear state
**Fix:** Settings ARE persisting in `chrome.storage.sync` - the issue is the content.js reload listener

### 2. Instant Settings Apply (No Page Reload)

#### Architecture:
```
popup.js  →  chrome.runtime.sendMessage  →  background.js  →  chrome.tabs.sendMessage  →  content.js
                                                                                              ↓
                                                                                    Re-render all images
```

#### Flow:
1. User toggles setting in popup
2. Popup saves to storage AND sends message to background
3. Background broadcasts to all tabs
4. Content.js receives message, updates settings object, re-renders all `.chemfig-diagram` images

#### Implementation:
- Add `reRenderAllImages()` function to content.js
- For each rendered image, store original SMILES/moleculeData in dataset
- On settings change, iterate images and re-render with new options

### 3. Smart Caching System

#### Two-Level Cache:

**Level 1: SMILES Cache (IndexedDB)**
- Key: Molecule name (normalized lowercase)
- Value: { smiles, source, cid, pdbid, codid, timestamp }
- Cross-page, persistent
- Purpose: Skip API lookup if we already know the SMILES

**Level 2: Rendered Image Cache (IndexedDB)**
- Key: SMILES + options hash (showCarbons, aromaticRings, etc.)
- Value: { svgDataUrl, timestamp }
- Cross-page, persistent
- Purpose: Skip SmilesDrawer rendering if exact same options

#### Cache Strategy:
```
chem:histamine:
    ↓
[Check SMILES Cache]
    ↓ hit                    ↓ miss
    use cached SMILES    →   Call IntegratedSearch API
                              ↓
                         Save to SMILES Cache
    ↓
[Check Image Cache with options hash]
    ↓ hit                    ↓ miss
    use cached SVG       →   Render with SmilesDrawer
                              ↓
                         Save to Image Cache
    ↓
Display image
```

#### Options Hash:
```javascript
function getOptionsHash(options) {
  return JSON.stringify({
    showCarbons: options.showCarbons,
    showAromaticRings: options.showAromaticRings,
    showHydrogens: options.showHydrogens,
    terminalCarbons: options.terminalCarbons,
    atomNumbering: options.atomNumbering,
    solidBondColors: options.solidBondColors
  });
}
```

#### Cache Invalidation:
- SMILES cache: Never expires (molecule names don't change)
- Image cache: Invalidate when options change OR manual clear
- Max entries: 1000 SMILES, 500 images (LRU eviction)

### 4. Database Schema (IndexedDB)

```javascript
// Database: ChemParserCache
// Store 1: smilesCache
{
  name: "histamine",           // Primary key (lowercase)
  smiles: "NCCc1c[nH]cn1",
  source: "pubchem",           // or "rcsb", "cod"
  cid: 774,                    // PubChem CID if available
  pdbid: null,                 // RCSB PDB ID if biomolecule
  codid: null,                 // COD ID if mineral
  type: "compound",            // "compound", "biomolecule", "mineral"
  timestamp: 1701936000000
}

// Store 2: imageCache
{
  key: "NCCc1c[nH]cn1_abc123", // SMILES + options hash
  smiles: "NCCc1c[nH]cn1",
  optionsHash: "abc123",
  svgDataUrl: "data:image/svg+xml;...",
  width: 300,
  height: 240,
  timestamp: 1701936000000
}
```

### 5. Implementation Order

1. **Phase 1: Instant Settings** (no cache yet)
   - Add message passing for settings changes
   - Add re-render function to content.js
   - Store moleculeData on each image for re-rendering

2. **Phase 2: SMILES Cache**
   - Create IndexedDB helper module
   - Cache SMILES lookups from IntegratedSearch
   - Check cache before API calls

3. **Phase 3: Image Cache**
   - Cache rendered SVGs with options hash
   - Check cache before SmilesDrawer rendering
   - LRU eviction for memory management

## Files to Modify

1. `popup.js` - Send messages on settings change
2. `background.js` - Relay messages to content scripts
3. `content.js` - Handle messages, re-render, caching
4. NEW: `cache.js` - IndexedDB helper module
