import requests
import json

url = "http://localhost:5001/m2cf/submit"
# Exact payload from extension logs
data = {
    "textAreaData": "OC1=CC=CC=C1", 
    "selections": [],
    "h2": "on"
}

print(f"Sending request to {url} with data: {data}")

try:
    response = requests.post(url, json=data)
    print(f"Status: {response.status_code}")
    result = response.json()
    print(f"SVG Link: {result.get('svglink')}")
    print(f"Chem Data (Response): {result.get('chem_data')}")
    
    # Check if we can download the SVG and see if it has hydrogens
    if result.get('svglink'):
        svg_url = f"http://localhost:5001{result.get('svglink')}"
        svg_resp = requests.get(svg_url)
        content = svg_resp.text
        # Count occurrences of ">H<" or similar text nodes in SVG which might indicate explicit H
        # This is a rough check
        h_count = content.count(">H<") + content.count(">H</") 
        print(f"SVG Content Length: {len(content)}")
        print(f"Potential H atoms found in SVG: {h_count}")
        
except Exception as e:
    print(f"Error: {e}")
