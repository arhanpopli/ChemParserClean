"""
Test RDKit's explicitMethyl option to see how it renders methyls
"""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def test_explicit_methyl():
    """Test explicitMethyl option with various molecules"""
    
    molecules = [
        ('CCC', 'Propane (2 methyls)'),
        ('CC(C)C', 'Isobutane (3 methyls)'),
        ('Cc1ccccc1', 'Toluene (1 methyl on benzene)'),
        ('CCCCCC', 'Hexane (2 methyls)')
    ]
    
    for smiles, name in molecules:
        print(f"\n{name}:")
        print("=" * 60)
        
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        
        # Test WITHOUT explicitMethyl
        drawer1 = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer1.DrawMolecule(mol)
        drawer1.FinishDrawing()
        svg1 = drawer1.GetDrawingText()
        
        # Test WITH explicitMethyl
        drawer2 = rdMolDraw2D.MolDraw2DSVG(400, 300)
        opts = drawer2.drawOptions()
        opts.explicitMethyl = True
        drawer2.DrawMolecule(mol)
        drawer2.FinishDrawing()
        svg2 = drawer2.GetDrawingText()
        
        # Check for CH3 in SVG
        has_ch3_without = 'CH' in svg1 or '>C<' in svg1
        has_ch3_with = 'CH' in svg2 or '>C<' in svg2
        
        print(f"  WITHOUT explicitMethyl: Contains 'CH' or 'C': {has_ch3_without}")
        print(f"  WITH explicitMethyl:    Contains 'CH' or 'C': {has_ch3_with}")
        
        # Count text elements
        text_count_without = svg1.count('<text')
        text_count_with = svg2.count('<text')
        
        print(f"  Text elements WITHOUT: {text_count_without}")
        print(f"  Text elements WITH:    {text_count_with}")
        
        # Save both versions for visual inspection
        with open(f'test_{name.split()[0].lower()}_without_methyl.svg', 'w', encoding='utf-8') as f:
            f.write(svg1)
        with open(f'test_{name.split()[0].lower()}_with_methyl.svg', 'w', encoding='utf-8') as f:
            f.write(svg2)
        
        print(f"  ✓ Saved: test_{name.split()[0].lower()}_without_methyl.svg")
        print(f"  ✓ Saved: test_{name.split()[0].lower()}_with_methyl.svg")

if __name__ == '__main__':
    test_explicit_methyl()
    print("\n" + "=" * 60)
    print("✓ Test complete! Check the SVG files to see the difference.")
