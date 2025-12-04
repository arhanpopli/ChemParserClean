"""
One-time patch for mol2chemfig Python 3 compatibility.
Fixes dict.values() sorting issue in molecule.py
"""
import os

file_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig', 'molecule.py')

print(f"Patching {file_path}...")

with open(file_path, 'r') as f:
    content = f.read()

# Fix line 314-315
old_code = "            atoms = self.atoms.values()\n            atoms.sort(key=lambda atom: len(atom.neighbors))"
new_code = "            atoms = list(self.atoms.values())\n            atoms.sort(key=lambda atom: len(atom.neighbors))"

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(file_path, 'w') as f:
        f.write(content)
    print("✅ Patched successfully!")
else:
    print("⚠️ Code not found - might already be patched or file changed")
