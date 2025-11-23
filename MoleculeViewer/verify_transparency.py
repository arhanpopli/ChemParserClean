
import sys
import os

# Add the current directory to path so we can import app
sys.path.insert(0, os.getcwd())

from app.chemistry import smiles_to_svg

def test_generation():
    smiles = "c1ccccc1" # Benzene
    error, svg = smiles_to_svg(smiles, width=300, height=300)
    
    if error:
        print(f"Error: {error}")
        return

    print("SVG Generated. Length:", len(svg))
    
    # Check for white background indicators
    if 'fill="#FFFFFF"' in svg or 'fill="white"' in svg or 'fill:#FFFFFF' in svg:
        print("FAIL: White fill found in SVG")
    else:
        print("PASS: No explicit white fill found")
        
    if '<rect' in svg:
        print("WARNING: <rect> tag found. Checking if it's a background...")
        # Simple check if it covers the whole image
        if 'width="300.0"' in svg and 'height="300.0"' in svg:
             print("FAIL: Full size rect found, likely background.")
    
    with open("debug_benzene.svg", "w") as f:
        f.write(svg)
    print("Saved to debug_benzene.svg")

if __name__ == "__main__":
    test_generation()
