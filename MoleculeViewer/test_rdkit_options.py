"""
Explore RDKit drawing options for showing carbons and methyls
"""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def test_all_carbon_options():
    """Test various RDKit options for showing carbons"""
    
    mol = Chem.MolFromSmiles('CCC')  # Propane
    AllChem.Compute2DCoords(mol)
    
    print("Testing RDKit drawing options for Propane (CCC):")
    print("=" * 70)
    
    # Test 1: Default
    print("\n1. DEFAULT (no options):")
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"   Text elements: {svg.count('<text')}")
    
    # Test 2: explicitMethyl
    print("\n2. explicitMethyl = True:")
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    opts.explicitMethyl = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"   Text elements: {svg.count('<text')}")
    
    # Test 3: addAtomIndices
    print("\n3. addAtomIndices = True:")
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    opts.addAtomIndices = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"   Text elements: {svg.count('<text')}")
    with open('test_atom_indices.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    
    # Test 4: atomLabelDeuteriumTritium
    print("\n4. atomLabelDeuteriumTritium = True:")
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    opts.atomLabelDeuteriumTritium = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"   Text elements: {svg.count('<text')}")
    
    # Test 5: Inspect all available options
    print("\n5. ALL AVAILABLE OPTIONS in MolDrawOptions:")
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    
    # Get all attributes
    relevant_attrs = [attr for attr in dir(opts) if not attr.startswith('_')]
    
    carbon_related = []
    for attr in relevant_attrs:
        attr_lower = attr.lower()
        if any(keyword in attr_lower for keyword in ['carbon', 'methyl', 'atom', 'label', 'explicit', 'show']):
            try:
                value = getattr(opts, attr)
                if not callable(value):
                    carbon_related.append((attr, value))
            except:
                pass
    
    print("\n   Carbon/Atom/Label related options:")
    for name, value in sorted(carbon_related):
        print(f"      {name:40s} = {value}")
    
    # Test 6: Try modifying the molecule to have explicit Hs
    print("\n6. Add explicit hydrogens:")
    mol_with_h = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol_with_h)
    
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    drawer.DrawMolecule(mol_with_h)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"   Text elements: {svg.count('<text')}")
    with open('test_explicit_h.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    
    # Test 7: Try dummy atoms
    print("\n7. Try using atomic labels directly:")
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            # Try to force the atom to display
            atom.SetProp('atomLabel', 'CH3' if atom.GetTotalNumHs() == 3 else 'C')
    
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    print(f"   Text elements: {svg.count('<text')}")
    with open('test_atom_labels.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    
    print("\n" + "=" * 70)
    print("✓ Tests complete!")

if __name__ == '__main__':
    test_all_carbon_options()
