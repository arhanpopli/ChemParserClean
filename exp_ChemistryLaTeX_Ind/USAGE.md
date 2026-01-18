# exp_ChemistryLaTeX_Ind - Client-Side Chemistry Renderer

## Usage Guide (Independent Version v7.0)

This is the **client-side version** of exp_ChemistryLaTeX_Ind that requires **NO server dependency**. All rendering happens locally using SmilesDrawer, and API queries go directly to PubChem, RCSB, COD, and OPSIN.

## Syntax

All chemistry tags use the strict `type=value` format:

```
chem:TYPE=VALUE:
```

### Supported Types

| Type | Description | Example |
|------|-------------|---------|
| `mol` | Compound by name (queries PubChem) | `chem:mol=benzene:` |
| `smiles` | Direct SMILES input | `chem:smiles=CCO:` |
| `biomol` | Biomolecule (queries RCSB) | `chem:biomol=hemoglobin:` |
| `mineral` | Mineral (queries COD) | `chem:mineral=quartz:` |
| `iupac` | IUPAC name (queries OPSIN) | `chem:iupac=2-methylpropan-1-ol:` |
| `pbdid` | Direct PDB ID | `chem:pbdid=4RHV:` |
| `cid` | Direct PubChem CID | `chem:cid=2244:` |
| `codid` | Direct COD ID | `chem:codid=1234567:` |

### Examples

```
chem:mol=aspirin:         → Renders aspirin from PubChem
chem:smiles=C1=CC=CC=C1:  → Renders benzene from direct SMILES
chem:biomol=insulin:      → Shows insulin structure from RCSB
chem:mineral=calcite:     → Shows calcite crystal from COD
chem:iupac=ethanol:       → Converts "ethanol" via OPSIN, renders SMILES
```

## Rendering Flags

Add flags after the compound value with `+` or `-`:

| Flag | Effect |
|------|--------|
| `+c` | Show carbon labels |
| `+o` | Show aromatic circles |
| `+n` | Show atom numbers |
| `+h` | Show hydrogens |
| `+m` | Show methyl groups |
| `+i` | Show implicit hydrogens (H₂, etc.) |
| `+p` | Flip horizontal |
| `+q` | Flip vertical |
| `+3d` | Show 3D viewer instead of 2D |

Example: `chem:mol=caffeine+c+o:` → Caffeine with carbons and aromatic circles

## Important Notes

1. **No `chem:name:` format** - You MUST use explicit type: `chem:mol=caffeine:` not `chem:caffeine:`
2. **Flags don't work with biomol/mineral** - Only `mol`, `smiles`, and `iupac` support rendering flags
3. **MineralNames.js fallback** - If COD is down, the extension uses a local mineral database
4. **SmilesDrawer rendering** - 2D structures are rendered client-side, no server needed
5. **Direct API queries** - PubChem, RCSB, COD, OPSIN are queried directly

## API Sources

- **PubChem** - Compound names → SMILES (100M+ compounds)
- **RCSB PDB** - Biomolecule names → PDB ID (250K+ structures)
- **COD** - Mineral names → COD ID (500K+ crystals)
- **OPSIN** - IUPAC names → SMILES

## Support

- **Discord**: [https://discord.gg/jRrAFKu9hf](https://discord.gg/jRrAFKu9hf)

## Version

ChemTex v7.0 (Independent Client-Side)