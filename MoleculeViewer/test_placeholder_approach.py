"""
Test the placeholder replacement approach for showing methyls at bond endpoints
"""
import sys
sys.path.insert(0, '.')

from app.chemistry import smiles_to_svg

def test_placeholder_methyls():
    """Test methyls appearing at bond endpoints like NO2"""
    
    test_cases = [
        ('CCC', 'Propane', {'show_methyls': True}),
        ('Cc1ccccc1', 'Toluene', {'show_methyls': True}),
        ('CC(C)C', 'Isobutane', {'show_methyls': True}),
        ('CCCCCC', 'Hexane', {'show_carbons': True}),
    ]
    
    print("Testing placeholder approach for methyl display:")
    print("=" * 70)
    
    for smiles, name, opts in test_cases:
        print(f"\n{name} ({smiles}) with {opts}:")
        
        error, svg = smiles_to_svg(smiles, 400, 300, opts)
        
        if error:
            print(f"  ❌ Error: {error}")
            import traceback
            traceback.print_exc()
        else:
            # Check for our labels in the SVG
            has_ch3 = 'CH₃' in svg or 'CH3' in svg
            has_ch2 = 'CH₂' in svg or 'CH2' in svg
            has_ch = '>CH<' in svg
            has_c = '>C<' in svg
            
            print(f"  ✓ SVG generated")
            print(f"    Contains CH₃: {has_ch3}")
            print(f"    Contains CH₂: {has_ch2}")
            print(f"    Contains CH: {has_ch}")
            print(f"    Contains C: {has_c}")
            
            # Save for visual inspection
            filename = f"test_placeholder_{name.lower()}_{list(opts.keys())[0]}.svg"
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(svg)
            print(f"    Saved: {filename}")
    
    print("\n" + "=" * 70)
    print("✓ Tests complete! Open SVG files in browser to verify positioning.")

if __name__ == '__main__':
    test_placeholder_methyls()
