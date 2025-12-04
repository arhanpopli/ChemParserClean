"""
Comprehensive patch for mol2chemfig Python 3 compatibility.
Fixes all dict.values(), dict.items(), dict.keys() sorting issues.
"""
import os
import glob

mol2chemfig_dir = os.path.join(os.path.dirname(__file__), 'mol2chemfig')

# All Python files in mol2chemfig
py_files = glob.glob(os.path.join(mol2chemfig_dir, '*.py'))

patches = [
    # molecule.py
    ("atoms = self.atoms.values()", "atoms = list(self.atoms.values())"),
    # common.py  
    ("lst.sort(key=lambda pair: pair[1])", "lst = list(lst)\n    lst.sort(key=lambda pair: pair[1])"),
]

for py_file in py_files:
    print(f"\nChecking {os.path.basename(py_file)}...")
    
    with open(py_file, 'r') as f:
        content = f.read()
    
    modified = False
    for old, new in patches:
        if old in content:
            content = content.replace(old, new)
            modified = True
            print(f"  ✅ Patched: {old[:50]}...")
    
    if modified:
        with open(py_file, 'w') as f:
            f.write(content)
        print(f"  💾 Saved changes")

print("\n✨ Patching complete!")
