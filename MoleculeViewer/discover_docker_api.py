#!/usr/bin/env python
"""Discover Docker mol2chemfig API parameters"""

import urllib.request, json

print("Testing Docker mol2chemfig Backend (port 8000)...\n")

# Try the main endpoint with extensive parameter testing
endpoints = [
    ('/api/smiles-to-svg', 'SMILES endpoint'),
    ('/api/name-to-svg', 'Name/nomenclature endpoint'),
    ('/api/mol-to-svg', 'MOL file endpoint'),
    ('/api/inchi-to-svg', 'InChI endpoint'),
    ('/api/query', 'Query/info endpoint'),
]

for endpoint_path, description in endpoints:
    print(f"Testing: {endpoint_path} ({description})")
    
    test_payload = {
        'smiles': 'c1ccccc1',
        'name': 'benzene',
        'mol': '',
        'inchi': '',
        'width': 400,
        'height': 400,
        # Visualization options from Docker version
        'compact': 0,
        'compact_view': 0,
        'fancy_bonds': True,
        'aromatic': True,
        'show_carbon': False,
        'show_methyl': False,
        'keep_hydrogens': False,
        'atom_numbers': False,
        'flip': 0,
        'rotate': 0,
        'indentation': 'keep',
        'wedge': True,
        'dash': True,
        '3d': True,
    }
    
    try:
        data = json.dumps(test_payload).encode()
        req = urllib.request.Request(
            f'http://localhost:8000{endpoint_path}',
            data=data,
            headers={'Content-Type': 'application/json'}
        )
        resp = urllib.request.urlopen(req)
        result = json.loads(resp.read().decode())
        
        print(f"  ✓ FOUND")
        print(f"    Response keys: {list(result.keys())}")
        if 'error' in result:
            print(f"    Error: {result['error']}")
        
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"  ✗ Not found (404)")
        elif e.code == 400:
            try:
                error_data = json.loads(e.read().decode())
                print(f"  ~ Bad request (400): {error_data}")
            except:
                print(f"  ~ Bad request (400)")
        else:
            print(f"  ~ HTTP {e.code}")
    except Exception as e:
        print(f"  ✗ Error: {type(e).__name__}: {str(e)[:80]}")
    
    print()
