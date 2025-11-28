#!/usr/bin/env python3
"""Quick test script for MolView API endpoints"""
import requests
import json

BASE_URL = "http://localhost:5003"

def test_api():
    print("=" * 50)
    print("Testing MolView API Endpoints")
    print("=" * 50)
    
    # Test 1: Root page
    print("\n1. Testing root page (/)...")
    try:
        r = requests.get(f"{BASE_URL}/", timeout=5)
        print(f"   Status: {r.status_code}")
        print(f"   Content type: {r.headers.get('content-type', 'unknown')}")
    except Exception as e:
        print(f"   ERROR: {e}")
    
    # Test 2: API detect-type
    print("\n2. Testing /api/detect-type/aspirin...")
    try:
        r = requests.get(f"{BASE_URL}/api/detect-type/aspirin", timeout=10)
        print(f"   Status: {r.status_code}")
        if r.status_code == 200:
            print(f"   Response: {r.json()}")
    except Exception as e:
        print(f"   ERROR: {e}")
    
    # Test 3: API search
    print("\n3. Testing /api/search?q=aspirin...")
    try:
        r = requests.get(f"{BASE_URL}/api/search?q=aspirin", timeout=15)
        print(f"   Status: {r.status_code}")
        if r.status_code == 200:
            data = r.json()
            print(f"   Response: {json.dumps(data, indent=2)[:500]}...")
        else:
            print(f"   Response body: {r.text[:200]}")
    except Exception as e:
        print(f"   ERROR: {e}")
    
    # Test 4: API search for protein
    print("\n4. Testing /api/search?q=insulin...")
    try:
        r = requests.get(f"{BASE_URL}/api/search?q=insulin", timeout=15)
        print(f"   Status: {r.status_code}")
        if r.status_code == 200:
            print(f"   Response: {r.json()}")
    except Exception as e:
        print(f"   ERROR: {e}")
    
    # Test 5: MolView with query parameter
    print("\n5. Testing /?q=aspirin (how MolView handles search)...")
    try:
        r = requests.get(f"{BASE_URL}/?q=aspirin", timeout=5)
        print(f"   Status: {r.status_code}")
        print(f"   Content type: {r.headers.get('content-type', 'unknown')}")
        print(f"   Response length: {len(r.text)} bytes")
    except Exception as e:
        print(f"   ERROR: {e}")
    
    print("\n" + "=" * 50)
    print("Tests complete!")
    print("=" * 50)

if __name__ == "__main__":
    test_api()
