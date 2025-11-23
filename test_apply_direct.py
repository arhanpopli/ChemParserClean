import requests
import json

url = "http://localhost:5001/m2cf/apply" # Note: calling apply directly
explicit_smiles = "[H][O][c]1[c]([H])[c]([H])[c]([H])[c]([H])[c]1[H]"
data = {
    "chem_data": explicit_smiles,
    "chem_format": "smiles",
    "chemfig": "",
    "selections": [],
    "h2": "keep" 
}

print(f"Sending explicit SMILES to {url}")

try:
    response = requests.post(url, json=data)
    print(f"Status: {response.status_code}")
    result = response.json()
    print(f"Chem Data (Response): {result.get('chem_data')}")
    
    if result.get('svglink'):
        svg_url = f"http://localhost:5001{result.get('svglink')}"
        svg_resp = requests.get(svg_url)
        content = svg_resp.text
        h_count = content.count(">H<") + content.count(">H</")
        print(f"Potential H atoms found in SVG: {h_count}")
        
except Exception as e:
    print(f"Error: {e}")
