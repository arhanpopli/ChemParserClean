#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test inorganic compound and 3D support"""

import urllib.request, json
import sys

# Fix encoding on Windows
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def test_compound(name, compound_type="nomenclature"):
    """Test a compound conversion"""
    print(f"\n{'='*60}")
    print(f"Testing: {name}")
    print(f"Type: {compound_type}")
    print(f"{'='*60}")
    
    data = json.dumps({
        'nomenclature': name if compound_type == 'nomenclature' else None,
        'smiles': name if compound_type == 'smiles' else None,
        'width': 400,
        'height': 400,
        'options': {
            'wedge_dash': True,
            '3d_mode': True
        }
    }).encode()
    
    endpoint = '/api/nomenclature-to-svg' if compound_type == 'nomenclature' else '/api/smiles-to-svg'
    
    req = urllib.request.Request(
        f'http://localhost:5000{endpoint}',
        data=data,
        headers={'Content-Type': 'application/json'}
    )
    
    try:
        resp = urllib.request.urlopen(req)
        result = json.loads(resp.read().decode())
        
        if result.get('error'):
            print(f"[X] Error: {result['error']}")
            return False
        
        print(f"[OK] Success!")
        print(f"  SMILES: {result.get('smiles', 'N/A')[:50]}")
        if 'info' in result:
            print(f"  Molecular weight: {result['info'].get('molecular_weight', 'N/A')}")
            print(f"  Formula: {result['info'].get('formula', 'N/A')}")
        print(f"  SVG size: {len(result.get('svg', ''))} bytes")
        
        return True
        
    except urllib.error.HTTPError as e:
        try:
            result = json.loads(e.read().decode())
            print(f"[X] HTTP {e.code}: {result.get('error', e.reason)}")
        except:
            print(f"[X] HTTP {e.code}: {e.reason}")
        return False
    except Exception as e:
        print(f"[X] Error: {e}")
        return False

print("INORGANIC COMPOUND SUPPORT TEST")
print("="*60)

# Test organic compounds (baseline)
test_compound('benzene')
test_compound('caffeine')

# Test inorganic coordination complexes
test_compound('potassium hexacyanoferrate(ii)')
test_compound('hexacyanoferrate(ii)')
test_compound('hexamminecobalt(iii) chloride')
test_compound('tetraaquadichlorochromium(iii) chloride')
test_compound('ferrocyanide')

# Test 3D with stereochemistry
test_compound('[C@H](F)(Cl)Br', 'smiles')  # Chiral center

print("\n" + "="*60)
print("TEST COMPLETE")
print("="*60)
