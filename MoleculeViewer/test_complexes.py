from app.chemistry import nomenclature_to_smiles

test_names = [
    "Potassium hexacyanoferrate(II)",
    "Hexaamminecobalt(III) chloride",
    "Sodium pentacyanonitrosylferrate(II)",
    "Tetrachloronickelate(II) ion",
    "Diaquatetrachlorochromate(III) ion",
    "Bis(ethylenediamine)nickel(II)",
    "Pentaamminechlororuthenium(III) ion",
]

print("\nTesting Coordination Complex Conversion:\n")
for name in test_names:
    error, smiles = nomenclature_to_smiles(name)
    status = "✓" if not error else "✗"
    print(f"{status} {name}")
    if smiles:
        print(f"  → {smiles}\n")
    else:
        print(f"  Error: {error}\n")
