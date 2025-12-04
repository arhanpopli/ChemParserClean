"""
Fix print statement in common.py for Python 3
"""
import os

file_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig', 'common.py')

print(f"Fixing print statement in {file_path}...")

with open(file_path, 'r') as f:
    content = f.read()

# Fix Python 2 print statement
old = 'print >> _debug_file, " ".join(args)'
new = 'print(" ".join(args), file=_debug_file)'

if old in content:
    content = content.replace(old, new)
    with open(file_path, 'w') as f:
        f.write(content)
    print("✅ Fixed print statement")
else:
    print("⚠️ Print statement not found or already fixed")
