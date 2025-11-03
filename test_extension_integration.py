#!/usr/bin/env python3
"""
Test script for Chrome Extension ↔ MoleculeViewer Integration
Verifies that the /api/render-smiles endpoint works correctly
"""

import requests
import json
import sys
from pathlib import Path

# Configuration
SERVER_URL = 'http://localhost:5000'
RENDER_ENDPOINT = f'{SERVER_URL}/api/render-smiles'
HEALTH_ENDPOINT = f'{SERVER_URL}/health'

# Test cases: (SMILES, description)
TEST_CASES = [
    ('CCO', 'Ethanol (simple alcohol)'),
    ('CCC', 'Propane (simple chain)'),
    ('CC(C)C', 'Isobutane (branched chain)'),
    ('c1ccccc1', 'Benzene (aromatic ring)'),
    ('CC(=O)C', 'Acetone (ketone)'),
    ('CC(=O)O', 'Acetic acid (carboxylic acid)'),
    ('C1CCCCC1', 'Cyclohexane (cyclic)'),
    ('CCN', 'Ethylamine (amine)'),
    ('CC(=O)N', 'Acetamide (amide)'),
]

def test_health():
    """Test if server is running"""
    print('\n' + '='*60)
    print('🔍 HEALTH CHECK')
    print('='*60)
    try:
        response = requests.get(HEALTH_ENDPOINT, timeout=5)
        if response.status_code == 200:
            print(f'✅ Server is running at {SERVER_URL}')
            print(f'   Response: {response.json()}')
            return True
        else:
            print(f'❌ Server returned status {response.status_code}')
            return False
    except requests.exceptions.ConnectionError:
        print(f'❌ Cannot connect to server at {SERVER_URL}')
        print(f'   Make sure MoleculeViewer is running:')
        print(f'   cd MoleculeViewer && python run_server.py')
        return False
    except Exception as e:
        print(f'❌ Error: {str(e)}')
        return False

def test_render_smiles():
    """Test /api/render-smiles endpoint"""
    print('\n' + '='*60)
    print('🧪 TESTING /api/render-smiles ENDPOINT')
    print('='*60)
    
    passed = 0
    failed = 0
    
    for smiles, description in TEST_CASES:
        try:
            print(f'\n  📝 {description}')
            print(f'     SMILES: {smiles}')
            
            # Make request with parameters
            params = {
                'smiles': smiles,
                'width': '300',
                'height': '200',
                'dark': 'false'
            }
            
            response = requests.get(RENDER_ENDPOINT, params=params, timeout=10)
            
            # Check response
            if response.status_code == 200:
                # Check if it's SVG
                if 'svg' in response.headers.get('Content-Type', '').lower():
                    svg_size = len(response.content)
                    print(f'     ✅ SUCCESS - SVG rendered ({svg_size} bytes)')
                    print(f'     Content-Type: {response.headers.get("Content-Type")}')
                    print(f'     Cache-Control: {response.headers.get("Cache-Control")}')
                    passed += 1
                else:
                    print(f'     ⚠️  Wrong content type: {response.headers.get("Content-Type")}')
                    failed += 1
            else:
                print(f'     ❌ FAILED - Status {response.status_code}')
                print(f'     Response: {response.text[:100]}')
                failed += 1
                
        except requests.exceptions.Timeout:
            print(f'     ❌ TIMEOUT - Request took too long')
            failed += 1
        except Exception as e:
            print(f'     ❌ ERROR - {str(e)}')
            failed += 1
    
    # Summary
    print('\n' + '-'*60)
    print(f'📊 Results: {passed} passed, {failed} failed out of {len(TEST_CASES)} tests')
    print('-'*60)
    
    return failed == 0

def test_dark_mode():
    """Test dark mode rendering"""
    print('\n' + '='*60)
    print('🌙 TESTING DARK MODE')
    print('='*60)
    
    try:
        print(f'\n  Testing CCO (ethanol) in dark mode...')
        
        # Light mode
        response_light = requests.get(RENDER_ENDPOINT, params={
            'smiles': 'CCO',
            'dark': 'false'
        }, timeout=10)
        
        # Dark mode
        response_dark = requests.get(RENDER_ENDPOINT, params={
            'smiles': 'CCO',
            'dark': 'true'
        }, timeout=10)
        
        if response_light.status_code == 200 and response_dark.status_code == 200:
            light_size = len(response_light.content)
            dark_size = len(response_dark.content)
            print(f'  ✅ Light mode: {light_size} bytes')
            print(f'  ✅ Dark mode: {dark_size} bytes')
            return True
        else:
            print(f'  ❌ One or both requests failed')
            return False
            
    except Exception as e:
        print(f'  ❌ ERROR - {str(e)}')
        return False

def test_cors_headers():
    """Test CORS headers"""
    print('\n' + '='*60)
    print('🔐 TESTING CORS HEADERS')
    print('='*60)
    
    try:
        response = requests.get(RENDER_ENDPOINT, params={'smiles': 'CCO'}, timeout=10)
        
        cors_header = response.headers.get('Access-Control-Allow-Origin')
        
        if cors_header:
            print(f'\n  ✅ CORS header found: {cors_header}')
            print(f'     Allows: {response.headers.get("Access-Control-Allow-Methods", "N/A")}')
            return True
        else:
            print(f'\n  ⚠️  No CORS header found')
            print(f'     Chrome extension may have trouble accessing the endpoint')
            return False
            
    except Exception as e:
        print(f'  ❌ ERROR - {str(e)}')
        return False

def test_error_handling():
    """Test error handling for invalid SMILES"""
    print('\n' + '='*60)
    print('⚠️  TESTING ERROR HANDLING')
    print('='*60)
    
    invalid_smiles = ['INVALID!!!', 'xxx', '']
    
    all_handled = True
    for smiles in invalid_smiles:
        try:
            print(f'\n  Testing with invalid SMILES: "{smiles}"')
            
            response = requests.get(RENDER_ENDPOINT, params={'smiles': smiles}, timeout=10)
            
            # Should return SVG even on error (with error message in SVG)
            if response.status_code in [200, 400, 500]:
                if 'svg' in response.headers.get('Content-Type', '').lower():
                    print(f'  ✅ Returns SVG with error message (Status {response.status_code})')
                else:
                    print(f'  ⚠️  Returns status {response.status_code} but not SVG')
                    all_handled = False
            else:
                print(f'  ⚠️  Unexpected status code: {response.status_code}')
                all_handled = False
                
        except Exception as e:
            print(f'  ❌ ERROR - {str(e)}')
            all_handled = False
    
    return all_handled

def main():
    """Run all tests"""
    print('\n' + '='*60)
    print('🧪 CHROME EXTENSION ↔ MOLECULEVIEWER INTEGRATION TEST')
    print('='*60)
    
    results = []
    
    # Test health
    if not test_health():
        print('\n❌ Server not running. Cannot continue.')
        sys.exit(1)
    
    # Run all tests
    results.append(('Render SMILES', test_render_smiles()))
    results.append(('Dark Mode', test_dark_mode()))
    results.append(('CORS Headers', test_cors_headers()))
    results.append(('Error Handling', test_error_handling()))
    
    # Final summary
    print('\n' + '='*60)
    print('📋 FINAL RESULTS')
    print('='*60)
    
    for test_name, passed in results:
        status = '✅ PASS' if passed else '❌ FAIL'
        print(f'  {status} - {test_name}')
    
    all_passed = all(result[1] for result in results)
    
    print('\n' + '='*60)
    if all_passed:
        print('🎉 ALL TESTS PASSED!')
        print('✅ Integration is working correctly')
        print('✅ Chrome extension can use MoleculeViewer server')
        print('\nNext steps:')
        print('  1. Load the chrome-extension/ folder in Chrome')
        print('  2. Select "🧪 MoleculeViewer Server" in rendering engine')
        print('  3. Open ChatGPT or any webpage with chemistry notation')
        print('  4. Watch structures render beautifully! 🧪✨')
    else:
        print('❌ SOME TESTS FAILED')
        print('Check the errors above and review the documentation')
        print('See: CHROME_EXTENSION_INTEGRATION.md')
    print('='*60 + '\n')
    
    return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())
