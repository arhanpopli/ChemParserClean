
import requests
import json
import sys

def check_endpoint():
    url = "http://localhost:5000/img/smiles?smiles=c1ccccc1&json=true"
    try:
        print(f"Requesting {url}...")
        response = requests.get(url)
        print(f"Status Code: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            svg = data.get('svg', '')
            print(f"SVG Length: {len(svg)}")
            
            # Check for white background
            if 'fill="#FFFFFF"' in svg or 'fill="white"' in svg or 'fill:#FFFFFF' in svg:
                print("FAIL: White fill found in SVG from server")
            else:
                print("PASS: No explicit white fill found in SVG from server")
                
            if '<rect' in svg:
                print("WARNING: <rect> tag found. Checking content...")
                print(svg[:500]) # Print first 500 chars to see the rect
            
            with open("server_response.svg", "w") as f:
                f.write(svg)
            print("Saved server response to server_response.svg")
        else:
            print("Error: Server returned non-200 status")
            print(response.text)
            
    except Exception as e:
        print(f"Error connecting to server: {e}")

if __name__ == "__main__":
    check_endpoint()
