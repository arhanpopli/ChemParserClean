"""
Test script to generate SVGs using the NATIVE mol2chemfig implementation.
Tests features: Aromatic Circles, Show Carbons, Show Methyls, Hydrogens.
"""
import os
import native_mol2chemfig

# Create output directory
OUTPUT_DIR = "native_test_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def generate_test_svg(name, smiles, options):
    print(f"\n🧪 Testing: {name}")
    print(f"   SMILES: {smiles}")
    print(f"   Options: {options}")
    
    svg_content = native_mol2chemfig.run_mol2chemfig(smiles, options)
    
    if svg_content:
        filename = os.path.join(OUTPUT_DIR, f"{name}.svg")
        with open(filename, "w", encoding="utf-8") as f:
            f.write(svg_content)
        print(f"✅ Saved: {filename}")
        return True
    else:
        print(f"❌ Failed to generate {name}")
        return False

# Define test cases
tests = [
    {
        "name": "1_benzene_aromatic",
        "smiles": "c1ccccc1",
        "options": {
            "m2cfAromaticCircles": True,
            "m2cfShowCarbons": False
        }
    },
    {
        "name": "2_toluene_methyls",
        "smiles": "Cc1ccccc1",
        "options": {
            "m2cfShowMethyls": True,
            "m2cfAromaticCircles": True
        }
    },
    {
        "name": "3_ethanol_carbons_hydrogens",
        "smiles": "CCO",
        "options": {
            "m2cfShowCarbons": True,
            "m2cfHydrogensMode": "add"  # Explicit hydrogens
        }
    },
    {
        "name": "4_caffeine_all_features",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "options": {
            "m2cfAromaticCircles": True,
            "m2cfShowMethyls": True,
            "m2cfShowCarbons": True,
            "m2cfFancyBonds": True
        }
    }
]

print("🚀 Starting Native Mol2ChemFig Tests...")
print(f"📂 Output Directory: {os.path.abspath(OUTPUT_DIR)}")

for test in tests:
    generate_test_svg(test["name"], test["smiles"], test["options"])

print("\n✨ Tests Completed!")
