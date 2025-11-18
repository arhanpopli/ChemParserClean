# 🧪 Chem Extension - Chemfig Support via mol2chemfig

## ✅ What's New

The Chem extension now supports **chemfig notation** powered by the mol2chemfig backend! This allows you to render LaTeX chemfig diagrams directly in ChatGPT or any webpage.

---

## 🎯 Supported Patterns

### Pattern 1: `\chemfig{...}` (Standard LaTeX)
**Syntax:**
```
\chemfig{CH_3-CH_2-OH}
\chemfig{C-C-C}:30
```

**Examples:**
- `\chemfig{CH_3-CH_2-OH}` → Renders ethanol structure
- `\chemfig{C-C-C}` → Renders propane structure  
- `\chemfig{*6(-=-=(-OH)-=)}:90` → Renders phenol rotated 90°

**How it works:**
1. Extension detects `\chemfig{...}` syntax
2. Converts chemfig code to SMILES using built-in converter
3. Sends SMILES to MoleculeViewer server
4. Renders as inline SVG image

---

### Pattern 2: `chem:\chemfig{...}` (Explicit Chemistry Notation)
**Syntax:**
```
chem:\chemfig{CH_3-CH_2-OH}
chem:\chemfig{C=C}:45
```

**Examples:**
- `chem:\chemfig{CH_3-CH_2-OH}` → Renders ethanol
- `chem:\chemfig{*6(=-=-=-)}` → Renders benzene ring
- `chem:\chemfig{C(-[1]NH_2)-[7]C(=O)OH}` → Renders amino acid

**Why use this?**
- More explicit - clearly marks it as chemistry content
- Works alongside existing `chem:benzene:` notation
- Future-proof syntax

---

### Pattern 3: `chem:text:` (Existing - Still Works!)
**Syntax:**
```
chem:benzene:
chem:CCO:
chem:2-methylpropane:
```

This pattern still works as before - converts chemical names or SMILES to structures.

---

## 🔧 How Chemfig Conversion Works

### Step 1: Chemfig Detection
```
Input: \chemfig{CH_3-CH_2-OH}
Detected: chemfig code
```

### Step 2: Chemfig → SMILES Conversion
```javascript
// Built-in converter in extension
chemfigToSmiles("CH_3-CH_2-OH")
  → "CCO"  // Ethanol SMILES
```

**Conversion Rules:**
- `CH_3-CH_2-OH` → `CCO` (ethanol)
- `C-C-C` → `CCC` (propane)
- `*6(=-=-=-)` → `c1ccccc1` (benzene)
- `C=C` → `C=C` (ethene)
- Removes angle brackets: `C-[1]C-[7]C` → `CCC`
- Handles branches: `C(-OH)` → `C(O)`

### Step 3: SMILES → SVG Rendering
```
MoleculeViewer Server:
  Input: CCO
  Output: <svg>...ethanol structure...</svg>
```

### Step 4: Display
```
Extension injects SVG inline with:
- Rotation support (`:30`, `:90`, etc.)
- Dark mode detection
- Proper sizing (300x200px)
- Lazy loading for performance
```

---

## 📋 Complete Example

### In ChatGPT:
```
Me: What is the structure of \chemfig{CH_3-CH_2-OH}?

Extension transforms to:
Me: What is the structure of [ETHANOL SVG IMAGE]?
```

### Mixed Usage:
```
Comparing structures:
- chem:benzene: ← Chemical name
- chem:c1ccccc1: ← SMILES notation
- \chemfig{*6(=-=-=-)} ← Chemfig code
- chem:\chemfig{*6(=-=-=-)} ← Explicit chemfig

All four render benzene! ✅
```

---

## ⚙️ Settings & Configuration

### Enable/Disable Chemfig Support
```javascript
// In extension settings
settings.renderChemfig = true;  // Enable (default)
settings.renderChemfig = false; // Disable
```

### Developer Mode
```javascript
settings.devMode = true;
// Shows raw chemfig text instead of rendering
// Useful for debugging
```

### Console Commands
```javascript
// Test chemfig conversion
window.chemRendererDebug.testChemfig()

// View conversion logs
window.chemRendererLogs()
```

---

## 🔬 Chemfig → SMILES Examples

| Chemfig | SMILES | Name |
|---------|--------|------|
| `CH_3-CH_2-OH` | `CCO` | Ethanol |
| `C-C-C` | `CCC` | Propane |
| `C=C` | `C=C` | Ethene |
| `*6(=-=-=-)` | `c1ccccc1` | Benzene |
| `C(-OH)-C` | `C(O)C` | Isopropanol (simplified) |
| `C(=O)C` | `CC=O` | Acetone |
| `C-[1]C-[7]C` | `CCC` | Propane (zigzag) |

---

## 🎨 Visual Comparison

### Before (Without Chemfig Support):
```
User types: \chemfig{CH_3-CH_2-OH}
Display: \chemfig{CH_3-CH_2-OH} (just text)
```

### After (With Chemfig Support):
```
User types: \chemfig{CH_3-CH_2-OH}
Display: [ETHANOL SVG IMAGE]
```

---

## 🚀 Performance

- ✅ Conversion happens **locally** (no external API call for conversion)
- ✅ Rendering uses **cached mol2chemfig backend** (Docker on port 8000)
- ✅ **Lazy loading** - only renders visible structures
- ✅ **Batch processing** - converts all chemfig on page at once
- ✅ **Memory efficient** - uses base64 encoding for data transfer

---

## 🐛 Troubleshooting

### Issue: Chemfig not rendering
**Solution:**
1. Check settings: `settings.renderChemfig` should be `true`
2. Open console: Look for "🧪 Applying Pattern 3: \\chemfig{...}"
3. Check mol2chemfig backend: Docker container on port 8000 should be running

### Issue: Wrong structure rendered
**Cause:** Chemfig syntax might be too complex for automatic conversion

**Solution:** Use SMILES directly instead:
```
Instead of: \chemfig{complex-structure}
Use: chem:YOUR_SMILES:
```

### Issue: Rotation not working
**Example:** `\chemfig{C-C-C}:90` not rotating

**Check:**
1. Rotation syntax correct (`:` followed by angle)
2. Console logs show "Rotation: 90°"
3. CSS transform applied to image

---

## 📝 Summary

### Three Ways to Render Chemistry:

1. **Chemical Names/SMILES:**
   ```
   chem:benzene:
   chem:CCO:
   ```

2. **Standard Chemfig:**
   ```
   \chemfig{CH_3-CH_2-OH}
   ```

3. **Explicit Chemfig:**
   ```
   chem:\chemfig{CH_3-CH_2-OH}
   ```

### All Three Use:
- ✅ MoleculeViewer Server (Flask on port 5000)
- ✅ mol2chemfig Backend (Docker on port 8000)
- ✅ SVG rendering with caching
- ✅ Worldwide shareable cache links

---

## 🎉 Result

**You can now use chemfig notation directly in ChatGPT!**

Just type `\chemfig{YOUR_STRUCTURE}` and watch it transform into a beautiful molecular diagram! 🧪

---

**Created:** November 6, 2025  
**Version:** Chem Extension v3.0+ with mol2chemfig support  
**Status:** Production Ready ✅
