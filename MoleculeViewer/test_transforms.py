#!/usr/bin/env python
import urllib.request, json

# Test without rotation
data1 = json.dumps({
    'smiles': 'c1ccccc1',
    'width': 300,
    'height': 300,
    'options': {}
}).encode()
req1 = urllib.request.Request('http://localhost:5000/api/smiles-to-svg', data=data1, headers={'Content-Type': 'application/json'})
resp1 = urllib.request.urlopen(req1)
result1 = json.loads(resp1.read().decode())
svg1 = result1['svg']

# Test with 90 degree rotation
data2 = json.dumps({
    'smiles': 'c1ccccc1',
    'width': 300,
    'height': 300,
    'options': {'rotate': 90}
}).encode()
req2 = urllib.request.Request('http://localhost:5000/api/smiles-to-svg', data=data2, headers={'Content-Type': 'application/json'})
resp2 = urllib.request.urlopen(req2)
result2 = json.loads(resp2.read().decode())
svg2 = result2['svg']

# Test with flip X
data3 = json.dumps({
    'smiles': 'c1ccccc1',
    'width': 300,
    'height': 300,
    'options': {'flip': 1}
}).encode()
req3 = urllib.request.Request('http://localhost:5000/api/smiles-to-svg', data=data3, headers={'Content-Type': 'application/json'})
resp3 = urllib.request.urlopen(req3)
result3 = json.loads(resp3.read().decode())
svg3 = result3['svg']

print("=== Test Results ===")
print(f"\n1. Normal SVG:")
print(f"   Contains transform: {'transform' in svg1}")
print(f"   First 200 chars: {svg1[:200]}")

print(f"\n2. With 90° rotation:")
print(f"   Contains transform: {'transform' in svg2}")
print(f"   Contains 'rotate(90': {'rotate(90)' in svg2}")
idx = svg2.find('transform')
if idx != -1:
    print(f"   Transform snippet: ...{svg2[max(0,idx-20):idx+60]}...")

print(f"\n3. With flip X:")
print(f"   Contains transform: {'transform' in svg3}")
print(f"   Contains 'scaleX': {'scaleX' in svg3}")
idx = svg3.find('transform')
if idx != -1:
    print(f"   Transform snippet: ...{svg3[max(0,idx-20):idx+60]}...")
