"""
Test the updated show_carbons and show_methyls implementation
"""
import sys
sys.path.insert(0, '.')

from app.chemistry import smiles_to_svg

def test_updated_implementation():
    """Test the RDKit atomLabels approach"""
    
    test_cases = [
        ('CCC', 'Propane', {'show_methyls': True}),
        ('CCC', 'Propane', {'show_carbons': True}),
        ('Cc1ccccc1', 'Toluene', {'show_methyls': True}),
        ('CCCCCC', 'Hexane', {'show_carbons': True}),
        ('CC(C)C', 'Isobutane', {'show_methyls': True}),
    ]
    
    print("Testing updated show_carbons and show_methyls implementation:")
    print("=" * 70)
    
    for smiles, name, opts in test_cases:
        print(f"\n{name} ({smiles}) with {opts}:")
        
        error, svg = smiles_to_svg(smiles, 400, 300, opts)
        
        if error:
            print(f"  ❌ Error: {error}")
        else:
            # Count atom path elements (RDKit renders text as paths)
            atom_paths = svg.count("class='atom-")
            bond_paths = svg.count("class='bond-")
            
            print(f"  ✓ SVG generated")
            print(f"    Atom paths: {atom_paths}")
            print(f"    Bond paths: {bond_paths}")
            
            # Save for visual inspection
            filename = f"test_updated_{name.lower()}_{list(opts.keys())[0]}.svg"
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(svg)
            print(f"    Saved: {filename}")
    
    print("\n" + "=" * 70)
    print("✓ All tests complete! Check the SVG files in a browser.")

if __name__ == '__main__':
    test_updated_implementation()
