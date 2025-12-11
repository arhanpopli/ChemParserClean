# Dynamic Image Reload System

## Overview
The extension now intelligently reloads only affected images when settings change, without requiring a full page reload. It uses lazy loading to prioritize visible images and shows visual indicators for pending updates.

## How It Works

### 1. **Intelligent Type Detection**
Images are categorized by type using the `data-compound-type` attribute:
- **`compound`**: Small molecules rendered with SmilesDrawer (benzene, histamine, phenol, etc.)
- **`biomolecule`**: Proteins and biological macromolecules (rhinovirus, hemoglobin, etc.)
- **`mineral`**: Crystal structures (NaCl, quartz, etc.)

### 2. **Settings-to-Type Mapping**
The `getAffectedImageTypes()` function determines which image types need reloading based on changed settings:

**SmilesDrawer Settings** (only affect `compound` type):
- `sdShowCarbons`, `sdShowMethyls`, `sdAromaticRings`
- `sdShowHydrogens`, `sdAtomNumbers`
- `sdFlipHorizontal`, `sdFlipVertical`
- `sdTheme`, `sdRotate`, `sdBondThickness`, `sdGradientColors`

**Protein/Biomolecule Settings** (only affect `biomolecule` and `mineral` types):
- `proteinRemoveWhiteBg`
- `viewer3DStyle`, `viewer3DAutoRotate`, `viewer3DBgColor`, `viewer3DSize`

**Universal Settings** (affect ALL types):
- `showTags`, `performanceMode`

### 3. **Lazy Loading Mode (Performance Mode ON)**

When lazy loading is enabled:

1. **Visible images** (in viewport) are re-rendered **immediately**
2. **Off-screen images** are:
   - Marked with `data-needs-re-render="true"`
   - Given a **pulsing red dot indicator** (top-left corner)
   - Monitored by an IntersectionObserver

3. When you **scroll** to a pending image:
   - It's automatically re-rendered
   - Red dot is removed
   - Pending count updates

4. **Global indicator** (top-left of page):
   - Shows "X pending" with a spinning clock icon
   - Click it to force-render ALL pending images

### 4. **Immediate Mode (Performance Mode OFF)**

When lazy loading is disabled:
- ALL affected images are re-rendered immediately
- Page may feel slower, but everything updates at once
- No pending indicators needed

## Example Scenarios

### Scenario 1: Change SmilesDrawer Settings
**Page has:** benzene, histamine (compounds), rhinovirus (protein), NaCl (mineral)

**You're viewing:** benzene, histamine (bottom of page)

**You toggle:** Show Carbons + Flip Horizontal

**What happens:**
- ✅ Benzene and histamine re-render **immediately** (visible)
- ❌ Rhinovirus and NaCl **ignored** (not affected by SmilesDrawer settings)
- No page reload needed!

### Scenario 2: Change Protein Settings
**Page has:** benzene, histamine (compounds), rhinovirus (protein), NaCl (mineral)

**You're viewing:** benzene, histamine

**You toggle:** Remove Background (protein setting)

**What happens:**
- ❌ Benzene and histamine **ignored** (not affected by protein settings)
- 🔴 Rhinovirus gets a **red dot** (pending, off-screen)
- 🔴 NaCl gets a **red dot** (pending, off-screen)
- When you scroll up: rhinovirus and NaCl re-render automatically

### Scenario 3: Performance Mode OFF
**You toggle:** Show Methyls (with lazy loading OFF)

**What happens:**
- ALL compound images re-render immediately
- Page might lag briefly during rendering
- No red dots or pending indicators
- Everything updates at once

## Visual Indicators

### Red Dot (Pulsing)
- **Location**: Top-left corner of each pending image
- **Meaning**: Image needs re-rendering with new settings
- **Removal**: Automatic when image scrolls into view

### Global Pending Indicator
- **Location**: Top-left of page (fixed position)
- **Shows**: "X pending" with spinning clock
- **Click Action**: Force-render all pending images
- **Auto-hides**: When no pending images remain

## Technical Details

### Cache Integration
The cache system (`generateImageCacheKey`) already includes all rendering options:
```javascript
`${smiles}__c1_m0_a1_h0_n0_light_300_240`
```

When settings change:
- New cache key is generated
- Old cached version is ignored
- New image is rendered and cached
- Future renders use new cached version

### IntersectionObserver
```javascript
rootMargin: '100px'  // Pre-load 100px before entering viewport
```

Images are re-rendered 100px before they become visible, ensuring smooth scrolling.

## Benefits

1. **No Page Reloads**: Settings apply instantly without losing scroll position
2. **Smart Filtering**: Only affected image types are regenerated
3. **Performance**: Lazy loading prevents rendering off-screen images
4. **User Feedback**: Visual indicators show what's pending
5. **Cache Efficient**: Leverages existing cache system with setting-based keys
6. **Scroll-Aware**: Images update as you scroll to them

## Edge Cases Handled

- **Missing elements**: Null checks prevent crashes
- **No SMILES data**: Images without SMILES are skipped gracefully
- **SmilesDrawer unavailable**: Fallback prevents errors
- **Mixed content**: Compounds, proteins, and minerals coexist peacefully
- **Rapid setting changes**: Pending settings are tracked and applied correctly
