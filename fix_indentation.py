"""
Fix the indentation issue in common.py
"""
import os

file_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig', 'common.py')

print(f"Fixing {file_path}...")

with open(file_path, 'r') as f:
    lines = f.readlines()

# Fix lines 106-108 (0-indexed: 105-107)
# Should be:
# line 106: lst = self._d.items()
# line 107: lst = list(lst)
# line 108: lst.sort(key=lambda pair: pair[1])

if len(lines) >= 108:
    # Find the problematic line
    for i, line in enumerate(lines):
        if 'lst = list(lst)' in line and 'lst.sort' in lines[i+1]:
            # Fix indentation on next line
            lines[i+1] = '        ' + lines[i+1].lstrip()
            print(f"✅ Fixed indentation on line {i+2}")
            break

with open(file_path, 'w') as f:
    f.writelines(lines)

print("✅ Done!")
