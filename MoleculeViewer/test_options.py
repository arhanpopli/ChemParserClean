#!/usr/bin/env python
"""
Test script to verify all visualization options are working correctly
"""
import json
import urllib.request
import urllib.error

def test_api(smiles, options, description):
    """Test API with given SMILES and options"""
    payload = {
        'smiles': smiles,
        'width': 400,
        'height': 400,
        'options': options
    }
    
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request('http://localhost:5000/api/smiles-to-svg', 
                                 data=data, 
                                 headers={'Content-Type': 'application/json'})
    try:
        response = urllib.request.urlopen(req)
        result = json.loads(response.read().decode('utf-8'))
        error = result.get('error')
        has_svg = len(result.get('svg', '')) > 0
        
        status = "✓ PASS" if (error is None and has_svg) else "✗ FAIL"
        print(f"{status}: {description}")
        if error:
            print(f"  Error: {error}")
        return error is None and has_svg
    except Exception as e:
        print(f"✗ FAIL: {description}")
        print(f"  Exception: {str(e)}")
        return False

def main():
    print("=" * 60)
    print("Testing Visualization Options")
    print("=" * 60)
    
    # Base options
    base_options = {
        'show_carbons': False,
        'show_methyls': False,
        'aromatic_circles': False,
        'fancy_bonds': False,
        'atom_numbers': False,
        'hydrogens': 'keep',
        'flip_horizontal': False,
        'flip_vertical': False,
        'rotate': 0,
        'recalculate_coordinates': False
    }
    
    tests = [
        # Test 1: Benzene with aromatic circles
        ('c1ccccc1', 
         {**base_options, 'aromatic_circles': True},
         "Benzene with aromatic_circles=True"),
        
        # Test 2: Benzene without aromatic circles
        ('c1ccccc1', 
         {**base_options, 'aromatic_circles': False},
         "Benzene with aromatic_circles=False"),
        
        # Test 3: Ethylbenzene with show_carbons
        ('CCc1ccccc1',
         {**base_options, 'show_carbons': True},
         "Ethylbenzene with show_carbons=True"),
        
        # Test 4: Ethylbenzene with show_methyls
        ('CCc1ccccc1',
         {**base_options, 'show_methyls': True},
         "Ethylbenzene with show_methyls=True"),
        
        # Test 5: Both show_carbons and show_methyls
        ('CCc1ccccc1',
         {**base_options, 'show_carbons': True, 'show_methyls': True},
         "Ethylbenzene with show_carbons=True and show_methyls=True"),
        
        # Test 6: Flip horizontal
        ('CCc1ccccc1',
         {**base_options, 'flip_horizontal': True},
         "Ethylbenzene with flip_horizontal=True"),
        
        # Test 7: Flip vertical
        ('CCc1ccccc1',
         {**base_options, 'flip_vertical': True},
         "Ethylbenzene with flip_vertical=True"),
        
        # Test 8: Rotate 90 degrees
        ('c1ccccc1',
         {**base_options, 'rotate': 90},
         "Benzene with rotate=90"),
        
        # Test 9: Rotate 180 degrees
        ('c1ccccc1',
         {**base_options, 'rotate': 180},
         "Benzene with rotate=180"),
        
        # Test 10: Delete hydrogens
        ('CCO',
         {**base_options, 'hydrogens': 'delete'},
         "Ethanol with hydrogens='delete'"),
        
        # Test 11: Add hydrogens
        ('CCO',
         {**base_options, 'hydrogens': 'add'},
         "Ethanol with hydrogens='add'"),
        
        # Test 12: Keep hydrogens (default)
        ('CCO',
         {**base_options, 'hydrogens': 'keep'},
         "Ethanol with hydrogens='keep'"),
        
        # Test 13: Atom numbers
        ('CCc1ccccc1',
         {**base_options, 'atom_numbers': True},
         "Ethylbenzene with atom_numbers=True"),
        
        # Test 14: Fancy bonds
        ('c1ccccc1',
         {**base_options, 'fancy_bonds': True},
         "Benzene with fancy_bonds=True"),
        
        # Test 15: Complex options combination
        ('CCc1ccccc1',
         {**base_options, 
          'show_carbons': True, 
          'show_methyls': True,
          'aromatic_circles': True,
          'fancy_bonds': True,
          'flip_horizontal': True,
          'rotate': 90},
         "Ethylbenzene with multiple options enabled"),
    ]
    
    passed = 0
    failed = 0
    
    for smiles, options, description in tests:
        if test_api(smiles, options, description):
            passed += 1
        else:
            failed += 1
    
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0

if __name__ == '__main__':
    import sys
    sys.exit(0 if main() else 1)
