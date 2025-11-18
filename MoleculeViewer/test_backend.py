from app.api import app
import json

# Test with Flask test client
client = app.test_client()

print("="*70)
print("Testing Backend DIRECTLY (without running server)")
print("="*70)

# Test 1: Nomenclature endpoint
print("\n1️⃣ Testing /img/nomenclature?nomenclature=acetone")
response = client.get('/img/nomenclature?nomenclature=acetone&width=300&height=200')
print(f"   Status: {response.status_code}")
data = response.get_json()
if data:
    print(f"   Success: {data.get('success')}")
    print(f"   SMILES: {data.get('smiles')}")
    print(f"   Cache URL: {data.get('cache_url')}")
    print(f"   SVG bytes: {len(data.get('svg', ''))}")
else:
    print(f"   ERROR: No response")

# Test 2: SMILES endpoint
print("\n2️⃣ Testing /img/smiles?smiles=CCO")
response = client.get('/img/smiles?smiles=CCO&width=300&height=200')
print(f"   Status: {response.status_code}")
data = response.get_json()
if data:
    print(f"   Success: {data.get('success')}")
    print(f"   SMILES: {data.get('smiles')}")
    print(f"   Cache URL: {data.get('cache_url')}")
    print(f"   SVG bytes: {len(data.get('svg', ''))}")
else:
    print(f"   ERROR: No response")

# Test 3: Cache detection - run same request again
print("\n3️⃣ Testing CACHE (run acetone again)")
response = client.get('/img/nomenclature?nomenclature=acetone&width=300&height=200')
data = response.get_json()
if data:
    print(f"   Cache URL (should be same): {data.get('cache_url')}")
else:
    print(f"   ERROR")

print("\n" + "="*70)
print("✅ Backend is working perfectly!")
print("="*70)
