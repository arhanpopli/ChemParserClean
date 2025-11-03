import os
import json
import re
from pathlib import Path

chemdoodle_path = r"C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\ChemDoodleWeb-11.0.0"
data_folder = os.path.join(chemdoodle_path, "data")

print("=" * 60)
print("CHEMDOODLE DATA EXTRACTION")
print("=" * 60)

# List all files in data folder
if os.path.exists(data_folder):
    print(f"\n📁 Contents of {data_folder}:")
    for root, dirs, files in os.walk(data_folder):
        level = root.replace(data_folder, '').count(os.sep)
        indent = ' ' * 2 * level
        print(f'{indent}{os.path.basename(root)}/')
        subindent = ' ' * 2 * (level + 1)
        for file in files[:10]:
            print(f'{subindent}{file}')
        if len(files) > 10:
            print(f'{subindent}... and {len(files) - 10} more files')
else:
    print(f"\n❌ Data folder not found: {data_folder}")

# Look for MOL files
print("\n" + "=" * 60)
print("SEARCHING FOR MOL FILES")
print("=" * 60)

mol_files = []
for root, dirs, files in os.walk(data_folder):
    for file in files:
        if file.lower().endswith('.mol'):
            filepath = os.path.join(root, file)
            mol_files.append(filepath)

print(f"Found {len(mol_files)} MOL files")
if mol_files:
    print("Sample MOL files:")
    for f in mol_files[:5]:
        print(f"  - {os.path.basename(f)}")

# Look for JSON files
print("\n" + "=" * 60)
print("CHECKING FOR JSON DATA FILES")
print("=" * 60)

json_files = []
for root, dirs, files in os.walk(data_folder):
    for file in files:
        if file.lower().endswith('.json'):
            filepath = os.path.join(root, file)
            json_files.append(filepath)
            print(f"✅ JSON: {file}")

if not json_files:
    print("⚠️  No JSON files found")

# Check main folder
print("\n" + "=" * 60)
print("MAIN CHEMDOODLE FOLDER")
print("=" * 60)

if os.path.exists(chemdoodle_path):
    main_contents = os.listdir(chemdoodle_path)
    print(f"Total items: {len(main_contents)}")
    for item in sorted(main_contents)[:30]:
        item_path = os.path.join(chemdoodle_path, item)
        if os.path.isdir(item_path):
            print(f"📁 {item}/")
        else:
            print(f"📄 {item}")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"MOL files: {len(mol_files)}")
print(f"JSON files: {len(json_files)}")
print(f"ChemDoodle is primarily a VISUALIZATION library")
print("=" * 60)
