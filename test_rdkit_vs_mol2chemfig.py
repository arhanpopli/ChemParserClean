"""Generate comparison SVGs: RDKit vs mol2chemfig features"""
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

def generate_rdkit_svg(smiles, filename, show_carbons=False, aromatic_circles=False):
    """Generate SVG with RDKit"""
    mol = Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    opts = drawer.drawOptions()
    
    # Show carbons
    if show_carbons:
        opts.explicitMethyl = True
        # Also try these for showing all carbons
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                atom.SetProp('atomLabel', 'C')
    
    # Aromatic circles - RDKit uses Kekule vs aromatic rendering
    if aromatic_circles:
        # Don't kekulize - keep aromatic representation
        # RDKit will draw circles for aromatic rings
        pass
    else:
        # Kekulize to show explicit double bonds
        Chem.Kekulize(mol, clearAromaticFlags=True)
    
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    
    with open(filename, 'w') as f:
        f.write(svg)
    
    return filename

# Test molecules
tests = [
    ("Benzene", "c1ccccc1"),
    ("Toluene", "Cc1ccccc1"),
    ("Naphthalene", "c1ccc2ccccc2c1"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
]

print("🧪 RDKit Feature Test\n")
print("="*60)

for name, smiles in tests:
    print(f"\n{name} ({smiles}):")
    
    # 1. Default (Kekule - double bonds)
    f1 = generate_rdkit_svg(smiles, f"rdkit_test_output/{name}_default.svg", 
                            show_carbons=False, aromatic_circles=False)
    print(f"  ✅ Default (Kekule): {f1}")
    
    # 2. Aromatic circles
    f2 = generate_rdkit_svg(smiles, f"rdkit_test_output/{name}_circles.svg",
                            show_carbons=False, aromatic_circles=True)
    print(f"  ✅ Aromatic circles: {f2}")
    
    # 3. Show carbons
    f3 = generate_rdkit_svg(smiles, f"rdkit_test_output/{name}_carbons.svg",
                            show_carbons=True, aromatic_circles=False)
    print(f"  ✅ Show carbons: {f3}")
    
    # 4. Both
    f4 = generate_rdkit_svg(smiles, f"rdkit_test_output/{name}_both.svg",
                            show_carbons=True, aromatic_circles=True)
    print(f"  ✅ Both features: {f4}")

print("\n" + "="*60)
print("\n📁 Check rdkit_test_output/ folder for SVGs")
print("\n🔍 Compare with mol2chemfig to see if RDKit can replace it!")
