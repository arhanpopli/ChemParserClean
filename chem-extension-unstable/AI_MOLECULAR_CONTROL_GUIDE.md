# AI Molecular Control - Advanced Flags System

## Overview
The Chemistry Renderer extension supports advanced inline flags that allow AI assistants (like ChatGPT) to control individual molecule rendering options without requiring users to change their default settings.

## Basic Syntax
```
chem:molecule_name/flags:
```

## Flag Categories

### 1. **Default Settings Flag (`d`)**
- **Flag**: `d`
- **Usage**: `chem:phenol/d:`
- **Behavior**:
  - **When AI Molecular Control is ENABLED**: Only applies user's default settings if `d` is present
  - **When AI Molecular Control is DISABLED**: Always uses user's default settings (ignores `d`)
- **Purpose**: Allows AI to explicitly request user's default settings

### 2. **Display Flags (mol2chemfig only)**
These flags control what elements are shown in the molecular structure:

| Flag | Meaning | Example |
|------|---------|---------|
| `+c` | Show carbon atoms | `chem:benzene/+c:` |
| `+o` | Show aromatic circles | `chem:benzene/+o:` |
| `+m` | Show methyl groups (CH₃) | `chem:ethanol/+m:` |
| `+n` | Show atom numbers | `chem:caffeine/+n:` |

**Example**: `chem:histamine/+c+o+m+n:` - Shows carbons, aromatic circles, methyls, and atom numbers

### 3. **Override Flags (with defaults)**
When using `d` (defaults), you can selectively disable specific options:

| Flag | Meaning | Example |
|------|---------|---------|
| `d-c` | Use defaults BUT hide carbons | `chem:phenol/d-c:` |
| `d-o` | Use defaults BUT hide aromatic circles | `chem:benzene/d-o:` |
| `d-m` | Use defaults BUT hide methyls | `chem:ethanol/d-m:` |
| `d-n` | Use defaults BUT hide atom numbers | `chem:caffeine/d-n:` |

**Example**: `chem:phenol/d-c:` - Uses all default settings except hides carbon atoms

### 4. **Transformation Flags**
These work for all rendering engines (mol2chemfig, MoleculeViewer, PubChem):

| Flag | Meaning | Values | Example |
|------|---------|--------|---------|
| `-p` | Flip horizontal | N/A | `chem:water/-p:` |
| `-q` | Flip vertical | N/A | `chem:water/-q:` |
| `-i` or `+inv` | Invert colors | N/A | `chem:benzene/-i:` |
| `+s###` | Set size | 1-500 (percentage) | `chem:caffeine/+s200:` (2x size) |
| `+r###` | Rotate | 0-360 (degrees) | `chem:benzene/+r60:` (60° rotation) |

**Examples**:
- `chem:water/+s150+r45:` - 1.5x size, rotated 45°
- `chem:benzene/-p-q:` - Flipped horizontally and vertically
- `chem:ethanol/+s80+inv:` - 80% size with inverted colors

### 5. **3D Viewer Flag**
| Flag | Meaning | Notes |
|------|---------|-------|
| `+3d` | Show 3D interactive viewer | **Do not combine with other flags** - they won't work in 3D mode |

**Example**: `chem:caffeine/+3d:` - Shows interactive 3D model

⚠️ **Important**: When using `+3d`, do NOT include display flags like `+c`, `+o`, etc. They only work for 2D renderings.

### 6. **Engine Selection**
| Flag | Meaning | Notes |
|------|---------|-------|
| `+pubchem` | Force PubChem renderer | Uses PubChem's official images |

**Example**: `chem:aspirin/+pubchem:` - Uses PubChem instead of default renderer

## AI Molecular Control Setting

### When ENABLED (in extension popup):
- AI **must** include `d` flag to use user's default settings
- Without `d`, the extension uses minimal/clean rendering
- Example: `chem:benzene:` → Clean rendering, `chem:benzene/d:` → User's defaults

### When DISABLED (default):
- Extension **always** uses user's default settings
- The `d` flag is ignored
- Example: `chem:benzene:` → Always uses user's defaults

## Complete Examples

### Example 1: Custom Display Options
```
chem:histamine/+c+o+m+n:
```
Shows histamine with carbons, aromatic circles, methyls, and atom numbers visible.

### Example 2: Using Defaults with Override
```
chem:phenol/d-c:
```
Uses user's default settings but hides carbon atoms for this specific molecule.

### Example 3: Size and Rotation
```
chem:caffeine/+s150+r90:
```
Renders caffeine at 150% size, rotated 90 degrees.

### Example 4: Flipped and Inverted
```
chem:water/-p-q+inv:
```
Renders water flipped horizontally, vertically, and with inverted colors.

### Example 5: 3D Viewer
```
chem:aspirin/+3d:
```
Shows aspirin in an interactive 3D viewer (do not add other flags).

### Example 6: Complex Combination
```
chem:benzene/+c+o+s120+r45-p:
```
Shows benzene with carbons and aromatic circles visible, at 120% size, rotated 45°, and flipped horizontally.

## Rendering Engine Compatibility

| Flag Type | mol2chemfig | MoleculeViewer | PubChem |
|-----------|-------------|----------------|---------|
| Display (`+c`, `+o`, `+m`, `+n`) | ✅ | ❌ | ❌ |
| Size (`+s###`) | ✅ | ✅ | ✅ |
| Rotation (`+r###`) | ✅ | ✅ | ✅ |
| Flip (`-p`, `-q`) | ✅ | ✅ | ✅ |
| Invert (`-i`, `+inv`) | ✅ | ✅ | ✅ |
| 3D (`+3d`) | ✅ | ✅ | ✅ |
| PubChem (`+pubchem`) | ✅ | ✅ | ✅ |

## Best Practices for AI Assistants

1. **Use `+3d` alone**: Don't combine with display flags
   - ✅ `chem:caffeine/+3d:`
   - ❌ `chem:caffeine/+3d+c+o:`

2. **Respect AI Molecular Control setting**:
   - If enabled: Use `d` to apply user defaults
   - If disabled: User defaults always apply

3. **Combine flags logically**:
   - ✅ `chem:benzene/+c+o+s150:` (display + size)
   - ✅ `chem:water/d-c+r90:` (defaults override + rotation)

4. **Size values**:
   - `+s50` = 50% size (half)
   - `+s100` = 100% size (normal)
   - `+s200` = 200% size (double)

5. **Rotation values**:
   - `+r0` = No rotation
   - `+r90` = 90° clockwise
   - `+r180` = Upside down
   - `+r270` = 90° counter-clockwise

## Common Use Cases

### Show detailed structure for teaching
```
chem:benzene/+c+o+n+s150:
```

### Clean, minimal structure
```
chem:benzene:
```
(Without AI Molecular Control, or with it disabled)

### User's preferred style
```
chem:benzene/d:
```
(With AI Molecular Control enabled)

### Interactive 3D exploration
```
chem:caffeine/+3d:
```

### Rotated for better viewing angle
```
chem:cholesterol/+r45+s120:
```

### Inverted for dark backgrounds
```
chem:aspirin/+inv:
```

## Notes for Implementation

- Flags can appear in **any order** after the `/`
- Multiple flags are concatenated: `+c+o+m` or `+c+o+s150+r90`
- The `:` at the end is **required** to close the chem notation
- Flags are case-sensitive (use lowercase)
- Invalid flags are silently ignored
- The extension handles all rendering and transformations automatically
