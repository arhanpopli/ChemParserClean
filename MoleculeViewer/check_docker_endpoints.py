#!/usr/bin/env python
import urllib.request, json

try:
    # Try root endpoint
    resp = urllib.request.urlopen('http://localhost:8000/')
    html = resp.read().decode()
    print("Root endpoint response (first 1000 chars):")
    print(html[:1000])
    print("\n" + "="*60 + "\n")
    
    # Check for available endpoints in the HTML
    if 'href' in html or 'action' in html or 'endpoint' in html:
        print("Found HTML content with potential endpoint references")
except Exception as e:
    print(f"Error accessing root: {e}")

# Try common API pattern
try:
    resp = urllib.request.urlopen('http://localhost:8000/api/')
    print("API root response:")
    print(resp.read().decode()[:500])
except Exception as e:
    print(f"Error accessing /api/: {e}")

# Check for mol2chemfig specific endpoint
try:
    resp = urllib.request.urlopen('http://localhost:8000/m2cf/submit')
    print("\n/m2cf/submit response:")
    content = resp.read().decode()
    print(content[:500])
except Exception as e:
    print(f"\nError accessing /m2cf/submit: {type(e).__name__}")
