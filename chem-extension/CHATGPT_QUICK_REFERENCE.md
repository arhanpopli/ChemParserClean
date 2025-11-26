# Chemistry Renderer Extension - Quick Reference for ChatGPT

## Syntax
`chem:molecule_name/flags:`

## Essential Flags

### Display (mol2chemfig only)
- `+c` = show carbons
- `+o` = aromatic circles  
- `+m` = show methyls
- `+n` = atom numbers
- `+h` = add explicit hydrogens

### Transformations (all engines)
- `+s###` = size (e.g., `+s150` = 150% size)
- `+r###` = rotate degrees (e.g., `+r90` = 90°)
- `-p` = flip horizontal
- `-q` = flip vertical
- `-i` or `+inv` = invert colors (auto-applied in dark mode)

### Special
- `+3d` = 3D viewer (use alone, no other flags)
- `+pubchem` = force PubChem renderer
- `d` = use user's defaults
- `d-c`, `d-o`, etc. = defaults minus specific option

## Rendering Engines
- **MoleculeViewer**: Default, server-side (Python/RDKit)
- **mol2chemfig**: LaTeX/Chemfig based
- **PubChem**: Official NCBI images
- **Client-Side**: Offline rendering (Kekule.js) - *NEW!*

## AI Molecular Control Setting
- **Enabled**: Must use `d` to apply user defaults, otherwise clean render
- **Disabled**: Always uses user defaults, `d` ignored

## Examples
```
chem:histamine/+c+o+m+n:        # All display options
chem:phenol/d-c:                 # User defaults minus carbons
chem:caffeine/+s150+r90:         # 150% size, 90° rotation
chem:benzene/+3d:                # 3D viewer only
chem:water/-p-q:                 # Flipped both ways
chem:aspirin/d:                  # User's default settings
chem:ethanol/+h:                 # Show explicit hydrogens
```

## Rules
1. Flags in any order after `/`
2. End with `:`
3. Don't combine `+3d` with display flags
4. Size: 50=half, 100=normal, 200=double
5. Rotation: 0-360 degrees
6. Flags are case-insensitive (e.g., `+3D` works)
