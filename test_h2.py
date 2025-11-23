import requests
import json

url = "http://localhost:5001/m2cf/submit"
data = {
    "textAreaData": "c1ccccc1",
    "h2": "on"
}

try:
    response = requests.post(url, json=data)
    print(f"Status: {response.status_code}")
    print(f"Response: {response.text[:200]}...")
    result = response.json()
    print(f"SVG Link: {result.get('svglink')}")
except Exception as e:
    print(f"Error: {e}")
