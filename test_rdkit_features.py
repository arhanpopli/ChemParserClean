"""Test RDKit's aromatic circles and show carbons features"""
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os

def test_rdkit_features():
    """Generate SVGs with different RDKit options"""
    
    # Test molecule: benzene (aromatic)
    smiles = "c1ccccc1"
    mol = Chem.MolFromSmiles(smiles)
    
    output_dir = "rdkit_test_output"
    os.makedirs(output_dir, exist_ok=True)
    
    tests = [
        {
            "name": "1_default",
            "options": {},
            "description": "Default (no carbons, no circles)"
        },
        {
            "name": "2_show_carbons",
            "options": {"addAtomIndices": False, "explicitMethyl": True},
            "description": "Show carbons (explicitMethyl)"
        },
        {
            "name": "3_aromatic_circles",
            "options": {"useBWAtomPalette": False},
            "description": "Aromatic circles attempt"
        },
        {
            "name": "4_both",
            "options": {"explicitMethyl": True, "useBWAtomPalette": False},
            "description": "Both features"
        }
    ]
    
    print("🧪 Testing RDKit Features\n")
    print("Molecule: Benzene (c1ccccc1)\n")
    
    for test in tests:
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 250)
        
        # Apply options
        opts = drawer.drawOptions()
        for key, value in test["options"].items():
            if hasattr(opts, key):
                setattr(opts, key, value)
        
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        output_file = os.path.join(output_dir, f"{test['name']}.svg")
        with open(output_file, 'w') as f:
            f.write(svg)
        
        print(f"✅ {test['name']}: {test['description']}")
        print(f"   Saved to: {output_file}")
    
    # Test with a more complex molecule
    print("\n" + "="*50)
    print("Testing with Toluene (methylbenzene)")
    print("="*50 + "\n")
    
    smiles2 = "Cc1ccccc1"
    mol2 = Chem.MolFromSmiles(smiles2)
    
    tests2 = [
        {
            "name": "5_toluene_default",
            "options": {},
            "description": "Toluene - Default"
        },
        {
            "name": "6_toluene_show_carbons",
            "options": {"explicitMethyl": True},
            "description": "Toluene - Show carbons"
        }
    ]
    
    for test in tests2:
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 250)
        opts = drawer.drawOptions()
        for key, value in test["options"].items():
            if hasattr(opts, key):
                setattr(opts, key, value)
        
        drawer.DrawMolecule(mol2)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        output_file = os.path.join(output_dir, f"{test['name']}.svg")
        with open(output_file, 'w') as f:
            f.write(svg)
        
        print(f"✅ {test['name']}: {test['description']}")
        print(f"   Saved to: {output_file}")
    
    print(f"\n📁 All SVGs saved to: {os.path.abspath(output_dir)}")
    print("\n🔍 Check the SVGs to see if RDKit supports:")
    print("   1. Aromatic circles (circles inside benzene rings)")
    print("   2. Show carbons (explicit C labels)")

if __name__ == "__main__":
    test_rdkit_features()
