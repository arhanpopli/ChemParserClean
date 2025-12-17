# ChemTex Quick Reference Card

## Basic Syntax
```
chem:type=value:
chem:<name>type=value:
chem:type=value+flags:
chem:<name>type=value+flags:
```

## Types

| Type | Database | Example | Flags? |
|------|----------|---------|--------|
| `smiles` | Direct SMILES | `chem:smiles=CCO:` | ✅ Yes |
| `mol` | PubChem name | `chem:mol=benzene:` | ✅ Yes |
| `biomol` | RCSB PDB name | `chem:biomol=insulin:` | ❌ No |
| `mineral` | COD name | `chem:mineral=quartz:` | ✅ Yes |
| `pbdid` | PDB ID | `chem:pbdid=3I40:` | ❌ No |
| `cid` | PubChem CID | `chem:cid=2244:` | ✅ Yes |
| `codid` | COD ID | `chem:codid=9000869:` | ✅ Yes |

## Flags

| Flag | Effect | Example |
|------|--------|---------|
| `+c` | Show carbons | `chem:mol=benzene+c:` |
| `-c` | Hide carbons | `chem:mol=benzene-c:` |
| `+o` | Show aromatic circles | `chem:mol=benzene+o:` |
| `-o` | Hide aromatic circles | `chem:mol=benzene-o:` |
| `+n` | Show atom numbers | `chem:smiles=CCO+n:` |
| `-n` | Hide atom numbers | `chem:smiles=CCO-n:` |
| `+h` | Add hydrogens | `chem:smiles=C+h:` |
| `-h` | Remove hydrogens | `chem:smiles=C-h:` |
| `+p` | Flip horizontal | `chem:mol=benzene+p:` |
| `+q` | Flip vertical | `chem:mol=benzene+q:` |
| `+d` | Use defaults + flags | `chem:mol=benzene+d+c:` |

## Flag Behavior

### Without +d (Override Mode)
```
chem:mol=benzene+c+o:
```
- Starts with ALL features OFF
- Only shows carbons and aromatic circles
- Ignores user settings

### With +d (Additive Mode)
```
chem:mol=benzene+d+c+o:
```
- Starts with user's settings
- Adds carbons and aromatic circles
- Builds on user preferences

## Named Structures

```
chem:Ethanolsmiles=CCO:              → Tag shows "Ethanol"
chem:Aspirincid=2244:                → Tag shows "Aspirin"
chem:Insulinpbdid=3I40:              → Tag shows "Insulin"
chem:Benzenesmiles=c1ccccc1+c+o:     → Tag shows "Benzene" with flags
```

## Common Examples

### Simple Molecules
```
chem:smiles=CCO:                     # Ethanol (direct SMILES)
chem:Ethanolsmiles=CCO:              # Ethanol (with name)
chem:mol=ethanol:                    # Ethanol (PubChem lookup)
chem:cid=702:                        # Ethanol (direct CID)
```

### Teaching Mode (Show Everything)
```
chem:Methanesmiles=C+c+n+h:          # Methane with all features
chem:Benzenesmiles=c1ccccc1+c+o+n:   # Benzene with details
```

### Professional Mode (Clean)
```
chem:mol=benzene+d:                  # Use user's settings
chem:Aspirincid=2244+d:              # Aspirin with defaults
```

### Biomolecules (No Flags)
```
chem:pbdid=3I40:                     # Insulin
chem:Insulinpbdid=3I40:              # Insulin (with name)
chem:biomol=hemoglobin:              # Hemoglobin (search)
```

### Minerals
```
chem:codid=9000869:                  # Quartz
chem:Quartzcodid=9000869:            # Quartz (with name)
chem:mineral=diamond:                # Diamond (search)
```

## Settings

### Enable AI Flag Control
- **Location:** Popup → Main Controls
- **ON (default):** Flags work normally
- **OFF:** All flags ignored, use only popup settings

### Show Tags
- **Location:** Popup → Main Controls
- **ON:** Tags always visible
- **OFF:** Tags hidden, show on hover

## Tips for AI Assistants

### ✅ DO
- Use descriptive names: `chem:Aspirinmol=aspirin:`
- Use direct IDs when known: `chem:Aspirincid=2244:`
- Use +d for consistency: `chem:mol=benzene+d+c:`
- Provide clean names: `chem:Ethanolsmiles=CCO:`

### ❌ DON'T
- Use flags with pbdid: `chem:pbdid=3I40+c:` ← Won't work
- Mix different flag styles: `chem:mol=benzene+c-o/d:` ← Wrong
- Forget the colon: `chem:smiles=CCO` ← Won't parse
- Use spaces in names: `chem:Ethyl Alcoholsmiles=CCO:` ← Use "EthylAlcohol"

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Molecule not rendering | Check syntax, verify ID is valid |
| Flags not working | Check "Enable AI Flag Control" is ON |
| Wrong molecule shown | Use direct ID (cid/pbdid/codid) instead of name |
| Tags not visible | Check "Show Tags" setting or hover over molecule |
| Changes not applying | Toggle setting off/on, or use Reload All button |

## Advanced Combinations

```
# Named + Direct ID + Flags
chem:Aspirincid=2244+d+c+n:

# Teaching example with everything visible
chem:Ethanoismiles=CCO+c+n+h+o:

# Custom name + PubChem lookup + specific styling
chem:BenzoicAcidmol=benzoic acid+d+c+o:

# Multiple molecules with consistent styling
chem:Methanosmiles=C+d+c+n+h:
chem:Ethanolsmiles=CCO+d+c+n+h:
chem:Propanolsmiles=CCCO+d+c+n+h:
```

---

**Version:** 1.0 (December 16, 2025)
**Documentation:** See AI_FRIENDLY_FEATURES.md for full details
**Test Page:** AI_FEATURES_TEST.html
