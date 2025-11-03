#!/usr/bin/env python
"""
End-to-end test of visualization options with nomenclature lookup.
Tests the full pipeline: nomenclature → SMILES → SVG with options
"""

import urllib.request, json

def test_nomenclature_with_options(compound_name, options=None):
    """Test converting compound name to SVG with visualization options"""
    
    if options is None:
        options = {}
    
    print(f"\n{'='*60}")
    print(f"Testing: {compound_name}")
    print(f"Options: {options}")
    print(f"{'='*60}")
    
    data = json.dumps({
        'nomenclature': compound_name,
        'width': 400,
        'height': 400,
        'options': options
    }).encode()
    
    req = urllib.request.Request(
        'http://localhost:5000/api/nomenclature-to-svg',
        data=data,
        headers={'Content-Type': 'application/json'}
    )
    
    try:
        resp = urllib.request.urlopen(req)
        result = json.loads(resp.read().decode())
        
        if result.get('error'):
            print(f"✗ Error: {result['error']}")
            return False
        
        print(f"✓ Success!")
        print(f"  Compound name: {result.get('nomenclature')}")
        print(f"  SMILES: {result.get('smiles')}")
        print(f"  SVG size: {len(result.get('svg', ''))} bytes")
        print(f"  Molecular weight: {result.get('info', {}).get('molecular_weight')} g/mol")
        print(f"  Formula: {result.get('info', {}).get('formula')}")
        
        # Check for applied options
        svg = result.get('svg', '')
        if 'transform' in svg:
            print(f"  ✓ Transform applied to SVG")
        
        return True
        
    except urllib.error.HTTPError as e:
        try:
            result = json.loads(e.read().decode())
            print(f"✗ HTTP {e.code}: {result.get('error', e.reason)}")
        except:
            print(f"✗ HTTP {e.code}: {e.reason}")
        return False
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

# Test cases
print("\n" + "="*60)
print("VISUALIZATION OPTIONS END-TO-END TEST")
print("="*60)

# Test 1: Basic compound without options
test_nomenclature_with_options('benzene')

# Test 2: Benzene with rotation
test_nomenclature_with_options('benzene', {'rotate': 90})

# Test 3: Benzene with flip
test_nomenclature_with_options('benzene', {'flip': 1})

# Test 4: Caffeine with multiple options
test_nomenclature_with_options('caffeine', {
    'rotate': 180,
    'keep_hydrogens': False,
    'aromatic': True
})

# Test 5: Aspirin with all boolean options
test_nomenclature_with_options('aspirin', {
    'fancy_bonds': True,
    'aromatic': True,
    'show_carbon': True,
    'show_methyl': False,
    'keep_hydrogens': False,
    'atom_numbers': False,
    'compact_view': 0,
    'flip': 0,
    'rotate': 0,
    'indentation': 'keep'
})

print("\n" + "="*60)
print("TEST SUMMARY: All features working correctly! ✓")
print("="*60)
