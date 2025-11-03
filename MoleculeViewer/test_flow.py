from app.chemistry import smiles_to_svg, nomenclature_to_smiles

# Test nomenclature
name_err, name_smiles, name_source = nomenclature_to_smiles('acetone')
print(f"Nomenclature 'acetone':")
print(f"  Error: {name_err}")
print(f"  SMILES: {name_smiles}")
print(f"  Source: {name_source}")

# Test SVG generation
svg_err, svg_data = smiles_to_svg(name_smiles if name_smiles else 'CC(=O)C')
print(f"\nSVG Generation:")
print(f"  Error: {svg_err}")
print(f"  SVG Length: {len(svg_data) if svg_data else 'None'}")
if svg_data:
    print(f"  First 200 chars: {svg_data[:200]}")
