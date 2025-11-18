"""
Test script for nomenclature parser configuration.
Tests different parser modes with various compound names.
"""

import sys
sys.path.insert(0, '.')
from app import config
from app.chemistry import nomenclature_to_smiles

print("="*70)
print("NOMENCLATURE PARSER CONFIGURATION TEST")
print("="*70)

# Test compounds covering different scenarios
test_compounds = [
    # Common names (ChemDoodle/Fallback/PubChem)
    ("aspirin", "Common name"),
    ("caffeine", "Common name"),
    ("glucose", "Common name"),
    
    # IUPAC names (OPSIN)
    ("2-methylpropane", "IUPAC systematic"),
    ("ethanol", "IUPAC/Common"),
    ("benzene", "IUPAC/Common"),
    
    # Complex/Obscure (PubChem needed)
    ("adrenaline", "Trade name"),
    ("paracetamol", "Drug name"),
]

print(f"\nCurrent Configuration:")
print(f"  Parser Mode: {config.NOMENCLATURE_PARSER}")
print(f"  Timeout: {config.PARSER_TIMEOUT}s")
if config.NOMENCLATURE_PARSER == 'auto':
    print(f"  Enabled Parsers:")
    print(f"    ChemDoodle: {config.ENABLE_CHEMDOODLE}")
    print(f"    OPSIN: {config.ENABLE_OPSIN}")
    print(f"    Fallback: {config.ENABLE_FALLBACK}")
    print(f"    PubChem: {config.ENABLE_PUBCHEM}")
print()

print("-"*70)
print("Testing Compounds:")
print("-"*70)

success_count = 0
for compound, category in test_compounds:
    print(f"\n{compound} ({category})")
    error, smiles, source = nomenclature_to_smiles(compound)
    
    if error:
        print(f"  ✗ Error: {error}")
    else:
        print(f"  ✓ Success!")
        print(f"    SMILES: {smiles}")
        print(f"    Source: {source}")
        success_count += 1

print("\n" + "="*70)
print(f"Results: {success_count}/{len(test_compounds)} compounds successfully parsed")
print("="*70)

# Show configuration recommendations
print("\nConfiguration Recommendations:")
print("-"*70)

if config.NOMENCLATURE_PARSER == 'auto':
    if success_count == len(test_compounds):
        print("✓ Current configuration is working well!")
    elif success_count < len(test_compounds) / 2:
        print("⚠ Many failures. Consider:")
        print("  - Enable all parsers (ENABLE_CHEMDOODLE/OPSIN/FALLBACK/PUBCHEM = True)")
        print("  - Check internet connection for PubChem")
        print("  - Install Java for OPSIN")
elif config.NOMENCLATURE_PARSER == 'chemdoodle':
    print(f"Using ChemDoodle only. Coverage: {success_count}/{len(test_compounds)}")
    if success_count < len(test_compounds):
        print("  → Try 'auto' mode for better coverage")
elif config.NOMENCLATURE_PARSER == 'opsin':
    print(f"Using OPSIN only. Coverage: {success_count}/{len(test_compounds)}")
    if success_count < len(test_compounds):
        print("  → OPSIN works best with IUPAC names")
        print("  → Try 'auto' mode for common names")
elif config.NOMENCLATURE_PARSER == 'pubchem':
    print(f"Using PubChem only. Coverage: {success_count}/{len(test_compounds)}")
    print("  → Good coverage but requires internet")
elif config.NOMENCLATURE_PARSER == 'fallback':
    print(f"Using Fallback only. Coverage: {success_count}/{len(test_compounds)}")
    if success_count < len(test_compounds):
        print("  → Fallback dictionary is very limited")
        print("  → Try 'auto' mode for better coverage")

print("\n" + "="*70)
print("To change configuration, edit: app/config.py")
print("Then restart the server for changes to take effect.")
print("="*70)
