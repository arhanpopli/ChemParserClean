"""
Extract chemical structures from ChemDoodle MOL files
Convert to SMILES and create a lookup database for MoleculeViewer
"""

# NOTE: small automated edit added 2025-11-03 — comment only, no behavior change

import os
import re
from pathlib import Path

# Try to import RDKit
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("⚠️  RDKit not available - will parse MOL files manually")

chemdoodle_path = r"C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\ChemDoodleWeb-11.0.0"
data_folder = os.path.join(chemdoodle_path, "data", "molecules")

print("=" * 70)
print("CHEMDOODLE MOL FILE EXTRACTION")
print("=" * 70)

# Dictionary to store extracted data
extracted_compounds = {}

# Find all MOL files
mol_files = []
if os.path.exists(data_folder):
    for root, dirs, files in os.walk(data_folder):
        for file in files:
            if file.lower().endswith('.mol'):
                filepath = os.path.join(root, file)
                mol_files.append(filepath)

print(f"\n✅ Found {len(mol_files)} MOL files")

# Process each MOL file
print("\nProcessing MOL files:")
print("-" * 70)

for mol_file in mol_files:
    filename = os.path.basename(mol_file)
    compound_name = filename.replace('.mol', '').replace('_', ' ').lower()
    
    try:
        if HAS_RDKIT:
            # Use RDKit to read MOL and convert to SMILES
            mol = Chem.MolFromMolFile(mol_file)
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                extracted_compounds[compound_name] = smiles
                print(f"✅ {filename:30s} → {smiles}")
            else:
                print(f"⚠️  {filename:30s} (invalid MOL format)")
        else:
            # Manual parsing: extract chemical name from MOL comment field
            with open(mol_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 3:
                    # Line 3 (index 3) often contains compound name
                    mol_name = lines[3].strip()
                    if mol_name:
                        compound_name = mol_name.lower()
                        extracted_compounds[compound_name] = f"[MOL from {filename}]"
                        print(f"ℹ️  {filename:30s} → {compound_name}")
    except Exception as e:
        print(f"❌ {filename:30s} (Error: {str(e)[:40]})")

print("\n" + "=" * 70)
print(f"Successfully extracted: {len(extracted_compounds)} compounds")
print("=" * 70)

# Display extracted compounds
if extracted_compounds:
    print("\n📊 EXTRACTED COMPOUNDS:\n")
    for name, smiles in sorted(extracted_compounds.items()):
        print(f"  '{name}': '{smiles}',")

# Generate Python code for MoleculeViewer integration
print("\n" + "=" * 70)
print("PYTHON CODE FOR MOLECULEVIEWER INTEGRATION")
print("=" * 70)

output_code = """
# ChemDoodle extracted compounds (add to chemistry.py)
CHEMDOODLE_COMPOUNDS = {
"""

for name, smiles in sorted(extracted_compounds.items()):
    output_code += f"    '{name}': '{smiles}',\n"

output_code += "}\n"

print("\n" + output_code)

# Save to file
output_file = r"C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\chemdoodle_compounds.py"
try:
    with open(output_file, 'w') as f:
        f.write("# ChemDoodle compound database extracted from MOL files\n")
        f.write("# Auto-generated - do not edit manually\n\n")
        f.write("CHEMDOODLE_COMPOUNDS = {\n")
        for name, smiles in sorted(extracted_compounds.items()):
            f.write(f"    '{name}': '{smiles}',\n")
        f.write("}\n")
    print(f"\n✅ Saved to: {output_file}")
except Exception as e:
    print(f"\n❌ Error saving file: {e}")

print("\n" + "=" * 70)
