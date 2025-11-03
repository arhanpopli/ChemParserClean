"""
Test to see what color/attributes dummy atoms have vs regular atoms
"""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# Test with a molecule that has a real functional group (OH) and a methyl
mol = Chem.MolFromSmiles('CCO')  # Ethanol - has OH and CH3
AllChem.Compute2DCoords(mol)

# Draw without modification (normal)
print("1. NORMAL (no modifications):")
drawer1 = rdMolDraw2D.MolDraw2DSVG(400, 300)
drawer1.DrawMolecule(mol)
drawer1.FinishDrawing()
svg1 = drawer1.GetDrawingText()

with open('test_normal_ethanol.svg', 'w', encoding='utf-8') as f:
    f.write(svg1)
print("   Saved: test_normal_ethanol.svg")

# Draw with methyl as dummy atom
print("\n2. WITH DUMMY ATOM (methyl replaced with E):")
mol2 = Chem.RWMol(mol)
# Replace the methyl carbon (atom 0) with dummy
atom0 = mol2.GetAtomWithIdx(0)
if atom0.GetSymbol() == 'C' and atom0.GetTotalNumHs() == 3:
    atom0.SetAtomicNum(0)  # Dummy
    atom0.SetProp('atomLabel', 'E')
    atom0.SetNoImplicit(True)

AllChem.Compute2DCoords(mol2)
drawer2 = rdMolDraw2D.MolDraw2DSVG(400, 300)
drawer2.DrawMolecule(mol2)
drawer2.FinishDrawing()
svg2 = drawer2.GetDrawingText()

with open('test_dummy_ethanol.svg', 'w', encoding='utf-8') as f:
    f.write(svg2)
print("   Saved: test_dummy_ethanol.svg")

print("\n3. Comparing atom-0 in both SVGs:")
print("   Looking for color differences...")

# Search for atom-0 in normal
import re
normal_atom0 = re.search(r"<path class='atom-0'[^>]*fill='([^']+)'", svg1)
dummy_atom0 = re.search(r"<path class='atom-0'[^>]*fill='([^']+)'", svg2)
normal_atom2 = re.search(r"<path class='atom-2'[^>]*fill='([^']+)'", svg1)

if normal_atom0:
    print(f"   Normal atom-0 (CH3): fill='{normal_atom0.group(1)}'")
else:
    print("   Normal atom-0: NOT FOUND (no label)")

if dummy_atom0:
    print(f"   Dummy atom-0 (E):    fill='{dummy_atom0.group(1)}'")
    
if normal_atom2:
    print(f"   Normal atom-2 (OH):  fill='{normal_atom2.group(1)}'")

print("\n✓ Check the SVG files to see the difference!")
