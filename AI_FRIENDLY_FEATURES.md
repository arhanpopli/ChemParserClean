# AI-Friendly ChemTex Features

## Overview
This document describes the new AI-friendly features added to ChemTex that make it easier for AI assistants like ChatGPT to generate chemical structure tags with custom names and direct ID lookups.

## New Syntax Features

### 1. Custom Named Structures
AI can now specify custom display names for molecules using the format:
```
chem:<name>type=value:
```

**Examples:**
- `chem:Ethanolsmiles=CCO:` - Displays "Ethanol" as the tag name
- `chem:Aspirinmol=aspirin:` - Displays "Aspirin" for the PubChem lookup
- `chem:Hemoglobinbiomol=hemoglobin:` - Displays "Hemoglobin" for biomolecule
- `chem:Quartzminer=quartz:` - Displays "Quartz" for mineral

### 2. Direct ID Lookups
AI can now directly specify database IDs instead of names for faster, more accurate lookups:

#### PDB ID (Protein Data Bank) - For Biomolecules
```
chem:pbdid=4RHV:
chem:Insulinpbdid=3I40:
```
- **Note:** pbdid does NOT support visual flags (biomolecules have different rendering)
- Use for proteins, DNA, RNA, and other biomolecules

#### CID (PubChem Compound ID) - For Compounds
```
chem:cid=1234:
chem:Aspirincid=2244:
```
- Supports all visual flags (+c, +n, +o, etc.)
- Faster than name lookups
- Avoids ambiguity

#### COD ID (Crystallography Open Database) - For Minerals
```
chem:codid=1234567:
chem:Quartzcodid=9000869:
```
- Supports visual flags
- Direct crystal structure lookup

### 3. Flag Syntax (Already Existed, Now Documented)
Flags use the `+` symbol for enable and `-` for disable:

**Visual Flags:**
- `+c` / `-c` - Show/hide carbons
- `+o` / `-o` - Show/hide aromatic circles
- `+n` / `-n` - Show/hide atom numbers
- `+h` / `-h` - Add/remove hydrogens
- `+p` / `-p` - Flip horizontal
- `+q` / `-q` - Flip vertical
- `+d` - Use default settings as base (important!)

**Examples:**
```
chem:mol=benzene+c+o:        # Show carbons and aromatic circles
chem:smiles=CCO+d+n:          # Use defaults, add atom numbers
chem:Aspirincid=2244+c+n-o:   # Show carbons and numbers, hide aromatic circles
```

**Flag Behavior:**
- **Without +d**: Flags completely override settings (all features off except those specified)
- **With +d**: Flags modify the user's default settings

## AI Flag Control Setting

### Purpose
The "Enable AI Flag Control" toggle allows users to control whether flags in chem tags override their settings.

### Location
Settings popup → Main Controls section

### Behavior
- **When ENABLED (default):** Flags in chem tags (like `+c`, `+n`) override user settings
- **When DISABLED:** All flags are ignored, only user's popup settings are used
- Updates all molecules **instantly** when toggled (no page reload needed)

### Use Cases
- **Enabled:** Good for AI-generated content where AI has specific visual preferences
- **Disabled:** User wants consistent styling across all molecules regardless of flags

## Tag Visibility Enhancement

### Behavior
The "Show Tags" toggle now has enhanced behavior:

- **When ENABLED:** Molecule name tags are always visible (opacity 0.7, becomes 1.0 on hover)
- **When DISABLED:** Tags are hidden by default (opacity 0), but appear on hover (opacity 0.7)

This allows users to keep their interface clean while still being able to see molecule names when needed.

## Complete Syntax Reference

### Standard Syntax (No Custom Name)
```
chem:smiles=CCO:              # Direct SMILES
chem:mol=benzene:             # PubChem compound lookup
chem:biomol=insulin:          # Biomolecule search
chem:mineral=quartz:          # Mineral search
chem:pbdid=4RHV:              # Direct PDB ID
chem:cid=1234:                # Direct PubChem CID
chem:codid=1234567:           # Direct COD ID
```

### Named Syntax (Custom Display Name)
```
chem:Ethanolsmiles=CCO:
chem:Aspirinmol=aspirin:
chem:Hemoglobinbiomol=hemoglobin:
chem:Quartzmineral=quartz:
chem:Insulinpbdid=3I40:
chem:Caffeinecid=2519:
chem:Quartzcodid=9000869:
```

### With Flags
```
chem:mol=benzene+c+o:                    # Standard with flags
chem:Benzenesmiles=c1ccccc1+d+c+n:       # Named with flags (use defaults)
chem:Aspirincid=2244+c+n-o:              # Direct CID with flags
chem:Quartzcodid=9000869+d+p:            # COD ID with defaults + flip
```

## Implementation Details

### Files Modified
1. **content.js**
   - Updated `parseChemFlags()` to support new syntax
   - Added direct ID handling in conversion logic
   - Implemented AI Flag Control check in rendering
   - Enhanced tag visibility CSS logic
   - Added `enableAIFlagControl` to universal settings for instant updates

2. **popup.html**
   - Updated "Show Tags" description to reflect new hover behavior
   - Corrected "Enable AI Flag Control" description syntax

3. **popup.js**
   - Already had `enableAIFlagControl` toggle implementation
   - Broadcasts changes for instant updates

### Caching Support
The caching system automatically works with all new ID types:
- CID lookups are cached by CID
- PBDID lookups are cached by PDB ID
- CODID lookups are cached by COD ID
- Named structures use their display name for cache keys

### Performance
- Direct ID lookups are faster than name searches
- Instant updates when toggling AI Flag Control
- Tag visibility changes apply immediately via CSS

## AI Assistant Guidelines

When generating ChemTex tags for users:

1. **Use descriptive names** when the chemical name is important for context
   ```
   chem:Aspirinmol=aspirin:
   ```

2. **Use direct IDs** when you know the database ID for accuracy
   ```
   chem:Aspirincid=2244:
   ```

3. **Add flags thoughtfully:**
   - Use `+d` when building on user preferences
   - Use specific flags without `+d` for precise control
   - Remember: pbdid doesn't support visual flags

4. **For biomolecules**, prefer pbdid when you know the PDB ID
   ```
   chem:Insulinpbdid=3I40:
   ```

5. **For common compounds**, CIDs provide the fastest, most accurate results
   ```
   chem:Caffeinecid=2519:
   chem:Ethanolcid=702:
   ```

## Examples for Common Use Cases

### Pharmaceutical Compounds
```
chem:Aspirinmol=acetylsalicylic acid+c+n:
chem:Ibuprofencid=3672+d+c:
```

### Proteins and Biomolecules
```
chem:Insulinpbdid=3I40:
chem:Hemoglobinpbdid=4HHB:
```

### Organic Chemistry
```
chem:Benzenesmiles=c1ccccc1+o:
chem:Ethanolsmiles=CCO+d+c:
```

### Minerals
```
chem:Quartzcodid=9000869:
chem:Diamondmineral=diamond:
```

### Educational Content (with flags)
```
chem:Methanesmiles=C+c+n:           # Show carbons and numbers for teaching
chem:Benzenesmiles=c1ccccc1+c+o+n:  # Full detail for aromatic systems
```
