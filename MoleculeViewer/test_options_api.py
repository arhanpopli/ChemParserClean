#!/usr/bin/env python
"""
Verify that options are being passed and cached correctly
"""
import requests
import json
import time
import sys

def test_api():
    print("\n" + "="*60)
    print(" MoleculeViewer Options Test")
    print("="*60 + "\n")
    
    tests = [
        {
            "name": "Benzene - No Options",
            "url": "http://localhost:5000/img/smiles?smiles=c1ccccc1",
            "expect_in_url": []
        },
        {
            "name": "Benzene - Aromatic Circles",
            "url": "http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true",
            "expect_in_url": ["aromatic_circles"]
        },
        {
            "name": "Benzene - Show Carbons",
            "url": "http://localhost:5000/img/smiles?smiles=c1ccccc1&show_carbons=true",
            "expect_in_url": ["show_carbons"]
        },
        {
            "name": "Benzene - Aromatic + Carbons",
            "url": "http://localhost:5000/img/smiles?smiles=c1ccccc1&aromatic_circles=true&show_carbons=true",
            "expect_in_url": ["aromatic_circles", "show_carbons"]
        },
        {
            "name": "Propane - Show Methyls",
            "url": "http://localhost:5000/img/smiles?smiles=CCC&show_methyls=true",
            "expect_in_url": ["show_methyls"]
        },
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            print(f"TEST: {test['name']}")
            print(f"  URL: {test['url'][:80]}...")
            
            response = requests.get(test['url'], timeout=10)
            response.raise_for_status()
            
            data = response.json()
            
            if not data.get('success'):
                print(f"  ❌ FAILED: API returned success=false")
                print(f"  Error: {data.get('error')}")
                failed += 1
                continue
            
            cache_url = data.get('cache_url', '')
            print(f"  Cache URL: {cache_url}")
            
            # Check if expected options are in the URL
            all_found = True
            for expected in test['expect_in_url']:
                if expected in cache_url:
                    print(f"  ✓ Found '{expected}' in cache URL")
                else:
                    print(f"  ✗ Missing '{expected}' in cache URL")
                    all_found = False
            
            if all_found:
                print(f"  ✅ PASSED\n")
                passed += 1
            else:
                print(f"  ❌ FAILED\n")
                failed += 1
                
        except requests.exceptions.ConnectionError:
            print(f"  ❌ FAILED: Cannot connect to server (is it running?)\n")
            failed += 1
        except Exception as e:
            print(f"  ❌ FAILED: {str(e)}\n")
            failed += 1
    
    print("="*60)
    print(f"Results: {passed} passed, {failed} failed")
    print("="*60 + "\n")
    
    return failed == 0

if __name__ == '__main__':
    # Wait for server to be ready
    print("Waiting for server to be ready...")
    for i in range(10):
        try:
            requests.head('http://localhost:5000')
            print("Server is ready!\n")
            break
        except:
            if i == 9:
                print("ERROR: Cannot reach server at localhost:5000")
                sys.exit(1)
            time.sleep(1)
    
    success = test_api()
    sys.exit(0 if success else 1)
