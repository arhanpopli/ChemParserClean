#!/usr/bin/env python
import urllib.request, json

data = json.dumps({
    'smiles': 'c1ccccc1',
    'width': 600,
    'height': 500,
    'options': {
        'rotate': 90,
        'fancy_bonds': True,
        'aromatic': True,
        'keep_hydrogens': False
    }
}).encode()

req = urllib.request.Request(
    'http://localhost:5000/api/smiles-to-svg',
    data=data,
    headers={'Content-Type': 'application/json'}
)

try:
    resp = urllib.request.urlopen(req)
    result = json.loads(resp.read().decode())
    print('✓ Success!')
    print(f'  Error: {result.get("error")}')
    print(f'  Has SVG: {bool(result.get("svg"))}')
    print(f'  SVG length: {len(result.get("svg", ""))} characters')
    print(f'  SMILES: {result.get("smiles")}')
    if 'transform' in result.get('svg', ''):
        print('  ✓ Transform found in SVG!')
    else:
        print('  Note: No transform in SVG (may be expected)')
except urllib.error.HTTPError as e:
    try:
        result = json.loads(e.read().decode())
        print(f'✗ HTTP {e.code}: {result}')
    except:
        print(f'✗ HTTP {e.code}: {e.reason}')
except Exception as e:
    print(f'✗ Error: {e}')
