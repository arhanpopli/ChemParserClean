import requests
import json

base_url = "http://localhost:5001"
implicit_smiles = "OC1=CC=CC=C1"
explicit_smiles = "[H][O][c]1[c]([H])[c]([H])[c]([H])[c]([H])[c]1[H]"

# Step 1: Submit implicit
print("Step 1: Submitting implicit SMILES...")
submit_data = {
    "textAreaData": implicit_smiles,
    "h2": "keep"
}
resp1 = requests.post(f"{base_url}/m2cf/submit", json=submit_data)
print(f"Submit Status: {resp1.status_code}")
res1 = resp1.json()

if resp1.status_code != 200:
    print("Submit failed")
    exit(1)

chemfig = res1.get('chemfig')
print(f"Got chemfig: {chemfig[:20]}...")

# Step 2: Apply explicit
print("Step 2: Applying explicit SMILES...")
apply_data = {
    "chem_data": explicit_smiles,
    "chem_format": "smiles",
    # "chemfig": chemfig, # Try without chemfig
    "selections": [],
    "h2": "keep" # We already added H manually
}

resp2 = requests.post(f"{base_url}/m2cf/apply", json=apply_data)
print(f"Apply Status: {resp2.status_code}")
res2 = resp2.json()
print(f"Chem Data (Response): {res2.get('chem_data')}")

if res2.get('svglink'):
    svg_url = f"{base_url}{res2.get('svglink')}"
    svg_resp = requests.get(svg_url)
    content = svg_resp.text
    h_count = content.count(">H<") + content.count(">H</")
    print(f"Potential H atoms found in SVG: {h_count}")
