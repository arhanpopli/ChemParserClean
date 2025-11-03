#!/usr/bin/env python
import urllib.request

try:
    resp = urllib.request.urlopen('http://localhost:8080')
    html = resp.read().decode('utf-8', errors='ignore')
    
    # Find form inputs, selects, checkboxes
    import re
    
    print("="*60)
    print("DOCKER MOL2CHEMFIG FRONTEND - HTML ANALYSIS")
    print("="*60)
    
    # Find all input elements
    inputs = re.findall(r'<input[^>]*?(?:name|id)=["\']([^"\']*)["\'][^>]*(?:type=["\']([^"\']*)["\'])?[^>]*>', html)
    if inputs:
        print("\nInput Elements Found:")
        for name, type_ in inputs[:20]:
            print(f"  - {name} ({type_ or 'text'})")
    
    # Find all select elements
    selects = re.findall(r'<select[^>]*?(?:name|id)=["\']([^"\']*)["\']', html)
    if selects:
        print("\nSelect Elements Found:")
        for name in selects[:20]:
            print(f"  - {name}")
    
    # Find all labels
    labels = re.findall(r'<label[^>]*>([^<]+)</label>', html)
    if labels:
        print("\nLabels Found:")
        for label in labels[:30]:
            clean = label.strip()
            if clean and len(clean) < 80:
                print(f"  - {clean}")
    
    # Find form data
    forms = re.findall(r'<form[^>]*>.*?</form>', html, re.DOTALL)
    if forms:
        print(f"\nForms Found: {len(forms)}")
        # Extract options/features mentioned in first 5000 chars
        snippet = html[:5000]
        if 'compact' in snippet:
            print("  ✓ 'compact' mentioned")
        if 'bond' in snippet:
            print("  ✓ 'bond' mentioned")
        if 'aromatic' in snippet:
            print("  ✓ 'aromatic' mentioned")
        if '3d' in snippet or '3D' in snippet:
            print("  ✓ '3D' mentioned")
        if 'wedge' in snippet:
            print("  ✓ 'wedge' mentioned")
        if 'dash' in snippet:
            print("  ✓ 'dash' mentioned")
    
    # Look for POST/submit endpoints
    endpoints = re.findall(r'action=["\']([^"\']*)["\']', html)
    if endpoints:
        print("\nForm Action Endpoints:")
        for endpoint in set(endpoints):
            print(f"  - {endpoint}")
    
    # Print a large section to manually inspect
    print("\n" + "="*60)
    print("HTML Sample (chars 1000-3000):")
    print("="*60)
    print(html[1000:3000])
    
except Exception as e:
    print(f"Error: {e}")
