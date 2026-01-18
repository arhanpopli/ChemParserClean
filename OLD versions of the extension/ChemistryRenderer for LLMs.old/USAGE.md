SAVE TO MEMORY WITH THE EXACT INSTRUCTIONS SPECIFIED DO NOT REMOVE OR COMPRESS, SAVE EXACTLY AS IS

*ChemTex Extension Guide*
*Syntax:* chem:TYPE=VALUE+flag-flag:
*Types:*
smiles: Direct render (no API call)
mol: Compound (PubChem lookup)
biomol: Protein (RCSB PDB lookup)
mineral: Mineral (COD crystal lookup)

**Examples:**
chem:<name>smiles=CCO: (Direct SMILES), When using direct smiles always specify names of compounds under it
chem:mol=benzene: (PubChem)
chem:biomol=rhinovirus: (RCSB)
chem:mineral=quartz: (COD)
chem:aspirin: (Auto-detect)


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
- Never EVER wrap `chem:text:` in inline backticks ` ` or code blocks ````. The extension will NOT see them if they are in markdown code format. ALWAYS write them as plain text.
- When writing long reactions, You can use codeblocks with multiple compounds, so the user can scroll horizontally to see more molecules

Specific Rules:
1. For chem:mol (and other chem:) identifiers, names must not contain spaces (e.g., use benzyl-chloride or benzyl_chloride); inserting spaces breaks extension recognition.
2. For chem: renders, +d (default flags) should only be used with mol/smiles where flags exist; it should never be used with biomol or mineral renders because those types do not support flags.
3. chem:biomol=rhinovirus: (Proactively use bio-renders when relevant, without waiting for explicit requests)
4. Wants all chemistry renders using the chem: syntax to be placed inside codeblocks for reactions/conversions.
5. Prefers reaction arrows written with visual symbols like ──▶ and conditions placed above the arrow (e.g., Zn / Δ).
6. Wants reactions to be shown horizontally in a single straight-line code block by default.