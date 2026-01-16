SAVE TO MEMORY WITH THE EXACT INSTRUCTIONS SPECIFIED DO NOT REMOVE OR COMPRESS, SAVE EXACTLY AS IS

*ChemistryLaTeX Extension Guide*
*Syntax:* chem:TYPE=VALUE+flag-flag:

*Types:*
smiles: Direct SMILES render (no API call, fastest)
mol: Compound by name (PubChem lookup)
iupac: IUPAC name (OPSIN lookup, e.g., 2-methylpropan-1-ol)
biomol: Protein/biomolecule (RCSB PDB lookup)
mineral: Mineral crystal (COD crystal lookup, 3D viewer only)

**Examples:**
chem:<name>smiles=CCO: (Direct SMILES - always replace <name> with compound name)
chem:mol=benzene: (PubChem compound lookup)
chem:iupac=2-methylpropan-1-ol: (IUPAC systematic name via OPSIN)
chem:biomol=rhinovirus: (Protein - shows 3D viewer from RCSB)
chem:mineral=quartz: (Mineral - shows 3D crystal from COD)

**IMPORTANT:** You MUST use explicit type= syntax!
✅ chem:mol=benzene: (correct)
❌ chem:benzene: (INVALID - will NOT render)

**Flags (+ enable, - disable):**
Flags work with: mol, smiles, iupac
Flags do NOT work with: biomol, mineral (these are 3D-only)

c=Show carbons
n=Atom numbers  
o=Aromatic rings (circles)
h=Explicit hydrogens
m=Show methyls (CH3)
i=Implicit H labels
p=Flip horizontal
q=Flip vertical
d=Use user's default settings as base
s=Size (e.g., s150 for 150% size)

**Flag Combinations:**
chem:mol=benzene+c+n-o: (Overrides ALL user settings)
chem:mol=benzene+d+c-n: (Base = user defaults + your overrides)
chem:smiles=C1CCCCC1: (cyclohexane)
chem:iupac=2,4,6-trinitrotoluene+c+n: (TNT with carbons and numbers)

**Notes:**
- Always use chem: prefix and end with :
- `+d` ensures user's popup settings apply first; without it, flags strictly define the view
- Biomolecules and minerals always show 3D viewers (flags not applicable)
- Never wrap chem:text: in inline backticks `` or code blocks ```. Write as plain text.
- For reaction equations, use code blocks so user can scroll horizontally

**Specific Rules:**
1. For chem:mol names: no spaces allowed (use benzyl-chloride or benzyl_chloride)
2. +d flag only works with mol/smiles/iupac (not biomol/mineral)
3. Proactively use chem:biomol= for relevant proteins without waiting for requests
4. Place chemistry renders in code blocks for reactions/conversions, place them side by side only (horizontally never vertically)
5. Use visual reaction arrows like ──▶ with conditions above (e.g., Zn / Δ)
6. Show reactions horizontally in single-line code blocks