from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

# Test if atomLabel works
mol = Chem.MolFromSmiles('CCCC')
AllChem.Compute2DCoords(mol)

# Try setting atomLabel
for i, atom in enumerate(mol.GetAtoms()):
    if atom.GetSymbol() == 'C':
        atom.SetProp('atomLabel', 'C')
        print(f'Atom {i}: Set label to C, has {atom.GetTotalNumHs()} hydrogens')

# Draw
drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
svg = drawer.GetDrawingText()

# Check if C appears in SVG
if '>C<' in svg or 'C</text>' in svg:
    print('\n✓ Carbon labels are visible in SVG!')
else:
    print('\n✗ Carbon labels NOT visible')
    print('\nChecking SVG content...')
    # Look for text elements
    import re
    text_matches = re.findall(r'<text[^>]*>(.*?)</text>', svg)
    print(f'Text elements found: {text_matches[:10]}')

with open('test_atomlabel.svg', 'w') as f:
    f.write(svg)
print('Saved to test_atomlabel.svg')
