import requests
import json
from rdkit import Chem

# Create MOL block for phenol with explicit H
mol = Chem.MolFromSmiles("OC1=CC=CC=C1")
mol_h = Chem.AddHs(mol)
mol_block = Chem.MolToMolBlock(mol_h)

url = "http://localhost:5001/m2cf/submit"
data = {
    "textAreaData": mol_block, 
    "selections": [],
    "h2": "keep" 
}

print(f"Sending MOL block to {url}")

try:
    response = requests.post(url, json=data)
    print(f"Status: {response.status_code}")
    if response.status_code != 200:
        print(f"Error Response: {response.text}")
    result = response.json()
    print(f"Chem Format: {result.get('chem_format')}")
    
    if result.get('svglink'):
        svg_url = f"http://localhost:5001{result.get('svglink')}"
        svg_resp = requests.get(svg_url)
        content = svg_resp.text
        h_count = content.count(">H<") + content.count(">H</")
        print(f"Potential H atoms found in SVG: {h_count}")
        
except Exception as e:
    print(f"Error: {e}")
