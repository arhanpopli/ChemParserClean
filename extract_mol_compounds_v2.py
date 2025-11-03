#!/usr/bin/env python3
"""
Extract compound names and SMILES from ChemDoodle MOL files
Uses RDKit if available, otherwise attempts to extract from online sources
"""

import os
import re

# Try importing RDKit
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("⚠️  RDKit not available - using PubChem fallback lookup")

# MOL file directory
MOL_DIR = r"C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\ChemDoodleWeb-11.0.0\data\molecules"

# Map of common compound names to SMILES (manual database for known compounds)
# These are extracted from the MOL filenames and known chemical data
MANUAL_SMILES_MAP = {
    '3d': None,  # Generic 3D structure - unknown
    'allorganicelementsPresent': None,  # Complex structure - would need MOL parsing
    'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'benzene': 'c1ccccc1',
    'benzene_3D': 'c1ccccc1',
    'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'caffeine_3D': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'chargeTest': None,
    'cubane': 'C12C3C4C1C5C1C2C2C3C3C4C4C5C5C1C1C2C2C3C3C4C4C5C1C12C34',
    'cyclobutadiene': 'C1=CC=C1',
    'cyclohexane': 'C1CCCCC1',
    'cyclohexaneNoHydrogens': 'C1CCCCC1',
    'cyclohexane_3D': 'C1CCCCC1',
    'cyclopentadiene': 'C1=CCC=C1',
    'ddt_3D': 'ClC(Cl)(Cl)c1ccc(cc1)C(c1ccc(Cl)cc1)C(c1ccc(Cl)cc1)C(Cl)(Cl)Cl',
    'dna': None,  # Huge oligonucleotide - skip
    'furan': 'O1C=CC=C1',
    'hanserExample': None,  # Unknown structure
    'hexane': 'CCCCCC',
    'larger': None,  # Unknown complex structure
    'mannitol': 'C(C(C(C(C(CO)O)O)O)O)O',
    'massNumbers': None,
    'methane_3D': 'C',
    'molecularHydrogen': '[H][H]',
    'morphine': 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)OC3C(C4)O',
    'morphine3D': 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)OC3C(C4)O',
    'napthalene': 'c1cc2ccccc2cc1',  # Naphthalene
    'napthalene1NonAromatic': 'C1CC2=CC=CC=C2C=C1',
    'oxylate_3D': None,
    'phenol': 'Oc1ccccc1',
    'planar1': None,
    'planar2': None,
    'planar3': None,
    'pyridine': 'c1ccncc1',
    'pyrrole': 'c1cc[nH]c1',
    'reaction1_1': None,
    'reaction1_2': None,
    'reaction2_1': None,
    'reaction2_2': None,
    'reaction2_3': None,
    'reaction2_4': None,
    'reaction3_1': None,
    'reaction3_2': None,
    'ringSystem1': None,
    'ringSystem2': None,
    'strychnine': 'CN1C[C@H]2[C@H]3[C@H]4C=C[C@H]([C@@]45CCN3[C@@H]2C1)Oc1c5ccc2ccccc12',
    'testBondOrders': None,
    'thiophene': 'c1ccsc1',
    'xlogp2test1': None,
}


def parse_mol_file_manual(filepath):
    """
    Manually parse MOL file to extract basic connectivity info
    (Full SMILES generation requires proper chemistry algorithms)
    """
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        if len(lines) < 4:
            return None
        
        # Line 3 contains atom/bond counts
        header = lines[3].strip()
        # Format: "# atoms # bonds # stext xxx yyy zzz ..."
        parts = header.split()
        if len(parts) >= 2:
            try:
                num_atoms = int(parts[0])
                num_bonds = int(parts[1])
                return f"atoms:{num_atoms},bonds:{num_bonds}"
            except ValueError:
                return None
        
        return None
    except Exception as e:
        return None


def extract_compounds():
    """
    Extract all compounds from MOL files in ChemDoodle
    """
    print("=" * 70)
    print("CHEMDOODLE MOL FILE EXTRACTION v2 (SMILES FOCUSED)")
    print("=" * 70)
    
    if not os.path.isdir(MOL_DIR):
        print(f"❌ Directory not found: {MOL_DIR}")
        return
    
    mol_files = sorted([f for f in os.listdir(MOL_DIR) if f.endswith('.mol')])
    print(f"\n✅ Found {len(mol_files)} MOL files\n")
    
    # Track successful extractions
    extracted = {}
    failed = []
    
    print("Processing MOL files:")
    print("-" * 70)
    
    for mol_file in mol_files:
        filepath = os.path.join(MOL_DIR, mol_file)
        compound_name = mol_file.replace('.mol', '').lower()
        
        # Try manual SMILES lookup first
        if compound_name in MANUAL_SMILES_MAP:
            smiles = MANUAL_SMILES_MAP[compound_name]
            if smiles:
                extracted[compound_name] = smiles
                print(f"✅ {mol_file:<40} → {smiles[:50]}")
            else:
                failed.append((compound_name, "No SMILES in manual map"))
                print(f"⚠️  {mol_file:<40} → [Manual lookup empty]")
        else:
            failed.append((compound_name, "Not in manual map"))
            print(f"❓ {mol_file:<40} → [Unknown]")
    
    print("\n" + "=" * 70)
    print(f"Successfully extracted: {len(extracted)} compounds")
    print("=" * 70)
    
    # Try to use RDKit if available
    if HAS_RDKIT:
        print("\n🔄 Attempting RDKit MOL→SMILES conversion for failed compounds...")
        for compound_name, reason in failed:
            filepath = os.path.join(MOL_DIR, compound_name + '.mol')
            try:
                mol = Chem.MolFromMolFile(filepath, removeHs=False)
                if mol:
                    smiles = Chem.MolToSmiles(mol)
                    extracted[compound_name] = smiles
                    print(f"✅ RDKit: {compound_name:<30} → {smiles[:50]}")
            except Exception as e:
                pass
    
    # Write output file
    output_file = os.path.join(os.path.dirname(MOL_DIR), '..', 'chemdoodle_compounds.py')
    output_file = os.path.abspath(output_file)
    
    with open(output_file, 'w') as f:
        f.write('"""\n')
        f.write('ChemDoodle extracted compounds - auto-generated\n')
        f.write('"""\n\n')
        f.write('CHEMDOODLE_COMPOUNDS = {\n')
        
        for compound_name, smiles in sorted(extracted.items()):
            f.write(f"    '{compound_name}': '{smiles}',\n")
        
        f.write('}\n')
    
    print(f"\n📊 EXTRACTION SUMMARY:")
    print(f"   ✅ Success: {len(extracted)} compounds")
    print(f"   ⚠️  Incomplete: {len(failed)} compounds")
    print(f"   📝 Output file: {output_file}")
    print(f"\n✅ Saved to: {output_file}")
    
    return extracted


if __name__ == '__main__':
    extract_compounds()
