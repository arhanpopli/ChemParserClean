*ChemTex Extension Guide*
*Syntax:* chem:TYPE=VALUE+flag-flag:
*Types:*
smiles: Direct render (no API call)
mol: Compound (PubChem lookup)
biomol: Protein (RCSB PDB lookup)
mineral: Mineral (COD crystal lookup)
[none]: Auto-detect (slower)

**Examples:**
chem:ethanolsmiles=CCO: (Named SMILES - "ethanol" appears as tag, CCO is the SMILES)
chem:mol=benzene: (PubChem lookup)
chem:biomol=rhinovirus: (RCSB PDB lookup)
chem:mineral=quartz: (COD crystal lookup)
chem:aspirin: (Auto-detect - slower)

**Named SMILES Syntax:**
Format: chem:<displayname>smiles=<SMILES>:
Example: chem:Ethanolsmiles=CCO: → displays "Ethanol" as the tag, renders CCO structure
Example: chem:Cyclohexanesmiles=C1CCCCC1: → displays "Cyclohexane", renders cyclohexane structure
This is useful when ChatGPT generates SMILES directly - specify the compound name before "smiles="


Flags (+ enable, - disable):
You can use these flags to show and explain stuff, like resonance using aromatic rings, Or using atom numbers or show carbons to explain nomenclature etc, They are your tools
c=Show carbons
n=Atom numbers
o=Aromatic rings
h=Explicit H
m=Show methyls
i=Implicit H labels
p=Flip Horiz
q=Flip Vert
d=Use Defaults Base
s= size, usecase s150

Combinations:
chem:mol=benzene+c+n-o: (Overrides ALL settings)
chem:mol=benzene+d+c-n: (Base defaults + overrides)
chem:smiles=C1CCCCC1: (cyclohexane)
Notes
- Always use chem: prefix and end with :
- `+d` ensures user's popup settings apply first, Without using +d overrides users default settings

- Without `+d`, flags strictly define the view (good for teaching)
- 3D mode falls back to 2D if 3D unavailable
- Never EVER wrap `chem:text:` in inline backticks NEVER USE THEM (pill-style inline code).
- When writing long reactions, You can use codeblocks with multiple compounds, so the user can scroll horizontally to see more molecules