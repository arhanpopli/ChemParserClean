"""
Test using atomLabels dictionary to show carbons and methyls
"""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def test_atom_labels():
    """Test atomLabels option"""
    
    mol = Chem.MolFromSmiles('CCC')  # Propane
    AllChem.Compute2DCoords(mol)
    
    print("Testing atomLabels for Propane (CCC):")
    print("=" * 70)
    
    # Check hydrogen counts
    print("\nAtom details:")
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        num_h = atom.GetTotalNumHs()
        print(f"  Atom {idx}: {symbol} with {num_h} hydrogens")
    
    # Test 1: Set custom labels for methyls
    print("\n1. Setting custom atomLabels for CH3 groups:")
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    
    # Set labels for atoms with 3 hydrogens (CH3)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 3:
            opts.atomLabels[idx] = 'CH₃'
            print(f"   Set atom {idx} label to 'CH₃'")
        elif atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 2:
            opts.atomLabels[idx] = 'CH₂'
            print(f"   Set atom {idx} label to 'CH₂'")
    
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"\n   Text elements in SVG: {svg.count('<text')}")
    print(f"   Contains 'CH': {'CH' in svg}")
    
    with open('test_propane_atomlabels.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    print("   ✓ Saved: test_propane_atomlabels.svg")
    
    # Test 2: Try with toluene (methyl on benzene)
    print("\n2. Testing with Toluene (Cc1ccccc1):")
    mol2 = Chem.MolFromSmiles('Cc1ccccc1')
    AllChem.Compute2DCoords(mol2)
    
    drawer2 = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts2 = drawer2.drawOptions()
    
    for atom in mol2.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 3:
            opts2.atomLabels[idx] = 'CH₃'
            print(f"   Set atom {idx} (methyl) label to 'CH₃'")
    
    drawer2.DrawMolecule(mol2)
    drawer2.FinishDrawing()
    svg2 = drawer2.GetDrawingText()
    print(f"\n   Text elements in SVG: {svg2.count('<text')}")
    print(f"   Contains 'CH': {'CH' in svg2}")
    
    with open('test_toluene_atomlabels.svg', 'w', encoding='utf-8') as f:
        f.write(svg2)
    print("   ✓ Saved: test_toluene_atomlabels.svg")
    
    # Test 3: Show all carbons
    print("\n3. Testing show ALL carbons with labels:")
    mol3 = Chem.MolFromSmiles('CCCCCC')  # Hexane
    AllChem.Compute2DCoords(mol3)
    
    drawer3 = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts3 = drawer3.drawOptions()
    
    for atom in mol3.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetSymbol() == 'C':
            num_h = atom.GetTotalNumHs()
            if num_h == 3:
                opts3.atomLabels[idx] = 'CH₃'
            elif num_h == 2:
                opts3.atomLabels[idx] = 'CH₂'
            elif num_h == 1:
                opts3.atomLabels[idx] = 'CH'
            else:
                opts3.atomLabels[idx] = 'C'
            print(f"   Atom {idx}: {opts3.atomLabels[idx]}")
    
    drawer3.DrawMolecule(mol3)
    drawer3.FinishDrawing()
    svg3 = drawer3.GetDrawingText()
    print(f"\n   Text elements in SVG: {svg3.count('<text')}")
    
    with open('test_hexane_all_labels.svg', 'w', encoding='utf-8') as f:
        f.write(svg3)
    print("   ✓ Saved: test_hexane_all_labels.svg")
    
    print("\n" + "=" * 70)
    print("✓ Test complete! Check SVG files.")

if __name__ == '__main__':
    test_atom_labels()
