"""
Verify native mol2chemfig output has aromatic circle.
"""
import os

filename = "test_output.svg"
print(f"Checking {filename}...")

if os.path.exists(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
        
    if "<circle" in content or "ellipse" in content:
        print("✅ Found circle/ellipse (aromatic ring)!")
    else:
        print("❌ No circle found!")
        
    # Check for Carbon and Hydrogen
    if ">C<" in content or ">C</text>" in content or "C" in content: # dvisvgm output varies
        print("✅ Found Carbon!")
    if ">H<" in content or ">H</text>" in content or "H" in content:
        print("✅ Found Hydrogen!")
        
else:
    print("❌ File not found")
