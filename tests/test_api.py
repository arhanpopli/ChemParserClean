"""
Autonomous test suite for MoleculeViewer API
Run: python test_api.py
"""
import requests
import json
import sys
from datetime import datetime

BASE_URL = "http://localhost:5000"
VERBOSE = True

class TestRunner:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.errors = []
    
    def test(self, name, condition, expected_msg=""):
        """Run a single test"""
        try:
            if condition:
                self.passed += 1
                print(f"[PASS] {name}")
                return True
            else:
                self.failed += 1
                msg = f"[FAIL] {name}"
                if expected_msg:
                    msg += f" ({expected_msg})"
                print(msg)
                self.errors.append(msg)
                return False
        except Exception as e:
            self.failed += 1
            msg = f"[ERROR] {name} - {str(e)}"
            print(msg)
            self.errors.append(msg)
            return False
    
    def report(self):
        """Print summary report"""
        total = self.passed + self.failed
        print(f"\n{'='*60}")
        print(f"TEST SUMMARY: {self.passed}/{total} passed")
        print(f"{'='*60}")
        if self.errors:
            print("\nFailed tests:")
            for error in self.errors:
                print(f"  {error}")
        return self.failed == 0

def test_smiles_with_options():
    """Test SMILES endpoint with various options"""
    print("\n" + "="*60)
    print("TESTING: SMILES Endpoint with Options")
    print("="*60)
    
    runner = TestRunner()
    
    # Test 1: Basic SMILES
    try:
        r = requests.get(f"{BASE_URL}/img/smiles?smiles=c1ccccc1&json=true", timeout=10)
        runner.test("Basic SMILES request", r.status_code == 200, f"Status: {r.status_code}")
        data = r.json()
        runner.test("Response has SVG", 'svg' in data or 'image_url' in data)
    except Exception as e:
        runner.test("Basic SMILES request", False, str(e))
    
    # Test 2: Aromatic circles
    try:
        r = requests.get(f"{BASE_URL}/img/smiles?smiles=c1ccccc1&aromatic_circles=true&json=true", timeout=10)
        runner.test("Aromatic circles option accepted", r.status_code == 200)
        data = r.json()
        cache_url = data.get('cache_url', '')
        runner.test("Cache URL contains 'aromatic_circles'", 
                   'aromatic_circles' in cache_url.lower(), 
                   f"URL: {cache_url}")
    except Exception as e:
        runner.test("Aromatic circles option", False, str(e))
    
    # Test 3: Show carbons
    try:
        r = requests.get(f"{BASE_URL}/img/smiles?smiles=c1ccccc1&show_carbons=true&json=true", timeout=10)
        runner.test("Show carbons option accepted", r.status_code == 200)
        data = r.json()
        cache_url = data.get('cache_url', '')
        runner.test("Cache URL contains 'show_carbons'", 
                   'show_carbons' in cache_url.lower(), 
                   f"URL: {cache_url}")
    except Exception as e:
        runner.test("Show carbons option", False, str(e))
    
    # Test 4: Multiple options
    try:
        r = requests.get(f"{BASE_URL}/img/smiles?smiles=CCC&show_carbons=true&show_methyls=true&json=true", timeout=10)
        runner.test("Multiple options accepted", r.status_code == 200)
        data = r.json()
        cache_url = data.get('cache_url', '')
        runner.test("Cache URL contains both options", 
                   'show_carbons' in cache_url.lower() and 'show_methyls' in cache_url.lower(), 
                   f"URL: {cache_url}")
    except Exception as e:
        runner.test("Multiple options", False, str(e))
    
    # Test 5: Nomenclature with options
    try:
        r = requests.get(f"{BASE_URL}/img/nomenclature?nomenclature=benzene&aromatic_circles=true&json=true", timeout=10)
        runner.test("Nomenclature with options", r.status_code == 200)
        data = r.json()
        cache_url = data.get('cache_url', '')
        runner.test("Nomenclature cache URL contains option", 
                   'aromatic_circles' in cache_url.lower(), 
                   f"URL: {cache_url}")
    except Exception as e:
        runner.test("Nomenclature with options", False, str(e))
    
    return runner.report()

def test_cache_system():
    """Test cache functionality"""
    print("\n" + "="*60)
    print("TESTING: Cache System")
    print("="*60)
    
    runner = TestRunner()
    
    # Test: Same request twice should use cache
    try:
        r1 = requests.get(f"{BASE_URL}/img/smiles?smiles=c1ccccc1&aromatic_circles=true&json=true", timeout=10)
        time1 = r1.elapsed.total_seconds()
        
        r2 = requests.get(f"{BASE_URL}/img/smiles?smiles=c1ccccc1&aromatic_circles=true&json=true", timeout=10)
        time2 = r2.elapsed.total_seconds()
        
        runner.test("Cache working (2nd request faster)", time2 < time1 * 1.5, f"T1:{time1:.3f}s, T2:{time2:.3f}s")
        
        data1 = r1.json()
        data2 = r2.json()
        runner.test("Same request returns same cache URL", 
                   data1.get('cache_url') == data2.get('cache_url'))
    except Exception as e:
        runner.test("Cache system", False, str(e))
    
    return runner.report()

def test_connection():
    """Test server is running"""
    print("\n" + "="*60)
    print("TESTING: Server Connection")
    print("="*60)
    
    runner = TestRunner()
    
    try:
        r = requests.get(f"{BASE_URL}/", timeout=5)
        runner.test(f"Server running at {BASE_URL}", r.status_code in [200, 404])
    except requests.exceptions.ConnectionError:
        runner.test("Server running", False, f"Could not connect to {BASE_URL}")
        return False
    except Exception as e:
        runner.test("Server connection", False, str(e))
        return False
    
    return runner.report()

if __name__ == "__main__":
    print(f"\n[TEST] AUTONOMOUS API TEST SUITE")
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Run all tests
    all_pass = True
    
    if not test_connection():
        print("\n[FAIL] Server not running! Cannot continue.")
        sys.exit(1)
    
    all_pass = test_smiles_with_options() and all_pass
    all_pass = test_cache_system() and all_pass
    
    # Final report
    print("\n" + "="*60)
    if all_pass:
        print("[SUCCESS] ALL TESTS PASSED!")
        sys.exit(0)
    else:
        print("[FAIL] SOME TESTS FAILED")
        sys.exit(1)
