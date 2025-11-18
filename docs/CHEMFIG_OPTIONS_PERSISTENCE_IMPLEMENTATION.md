# ChemFig Options Persistence - Implementation Guide

## Overview

This document explains the implementation of persistent chemfig options in the mol2chemfig_server.py Flask application. The solution fixes the issue where chemfig options don't persist when clicking "Search and Convert" again, and ensures cache links display properly.

## Problem Statement

**Before the fix:**
1. User selects chemfig options (aromatic circles, show carbon, etc.)
2. Clicks "Apply Options" - works once
3. Cache links don't show after applying options
4. When user searches again, options are not remembered
5. User has to manually re-apply options every time

**After the fix:**
1. User selects chemfig options once
2. Clicks "Save as Default Options" to persist preferences
3. Cache links are properly displayed
4. Future searches automatically use saved options
5. Options persist across page reloads (server-side session storage)

## Architecture Changes

### 1. Session-Based Options Storage

Added Flask session support to mol2chemfig_server.py:

```python
from flask import Flask, request, jsonify, send_file, send_from_directory, session
from flask_cors import CORS

app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'chemfig-server-secret-key-change-in-production')
CORS(app, supports_credentials=True)  # Enable credentials for session cookies

# Default options storage
default_options = {
    "selections": [],
    "angle": 0,
    "indentation": 4,
    "h2": "keep"
}
```

### 2. New API Endpoints

#### POST /api/options/save
Saves user's default chemfig options to the server session.

**Request:**
```json
{
    "selections": ["-o", "-m"],
    "angle": 0,
    "indentation": 4,
    "h2": "keep"
}
```

**Response:**
```json
{
    "success": true,
    "message": "Default options saved successfully",
    "options": { ... },
    "timestamp": "2025-11-09T..."
}
```

#### GET /api/options/get
Retrieves currently saved default options from session.

**Response:**
```json
{
    "success": true,
    "options": {
        "selections": ["-o", "-m"],
        "angle": 0,
        "indentation": 4,
        "h2": "keep"
    },
    "timestamp": "2025-11-09T..."
}
```

#### POST /api/options/clear
Clears saved options (resets to defaults).

**Response:**
```json
{
    "success": true,
    "message": "Default options cleared",
    "timestamp": "2025-11-09T..."
}
```

### 3. Modified /api/generate Endpoint

The `/api/generate` endpoint now automatically applies saved default options:

```python
@app.route('/api/generate', methods=['POST'])
def generate():
    """
    Generate molecule image from SMILES

    Now supports automatic application of saved default options!
    """
    data = request.get_json()
    smiles = data.get('smiles') or data.get('textAreaData', '')
    options = data.get('options', None)
    use_default_options = data.get('use_default_options', True)

    # If options not explicitly provided and use_default_options is true, use saved defaults
    if options is None and use_default_options:
        if 'default_options' in session:
            session_opts = session['default_options']
            options = session_opts.get('selections', [])
        else:
            options = default_options.get('selections', [])
    elif options is None:
        options = []

    # ... rest of generation logic ...

    # Response now includes applied_options
    return jsonify({
        "success": True,
        "hash": content_hash,
        "svg_url": f"/images/{content_hash}.svg",
        "pdf_url": f"/images/{content_hash}.pdf",
        "chemfig": chemfig,
        "cached": False,
        "applied_options": options,  # NEW: shows which options were used
        "timestamp": datetime.now().isoformat()
    })
```

### 4. Cache Link Display Fix

The response from `/api/generate` now properly includes cache links:

```python
# Cache links are returned in both cached and non-cached responses
return jsonify({
    "success": True,
    "hash": content_hash,
    "svg_url": f"/images/{content_hash}.svg" if svg_path else None,
    "pdf_url": f"/images/{content_hash}.pdf" if pdf_path else None,
    "chemfig": chemfig,
    "cached": True/False,
    "applied_options": options,
    "timestamp": datetime.now().isoformat()
})
```

## Client-Side Integration

### Test Page Implementation

The `test_options_persistence.html` demonstrates full integration:

```javascript
// Save options to server
async function saveDefaultOptions() {
    const selections = getSelectedOptions();

    const response = await fetch(`${SERVER_URL}/api/options/save`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        credentials: 'include',  // IMPORTANT: Include session cookie
        body: JSON.stringify({
            selections: selections,
            angle: 0,
            indentation: 4,
            h2: 'keep'
        })
    });

    const data = await response.json();
    // Options are now saved on server
}

// Generate molecule with auto-applied options
async function searchAndConvert() {
    const response = await fetch(`${SERVER_URL}/api/generate`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        credentials: 'include',  // IMPORTANT: Include session cookie
        body: JSON.stringify({
            smiles: smiles,
            format: 'smiles',
            use_default_options: true,  // Auto-apply saved options
            return_format: 'both'
        })
    });

    const genData = await response.json();
    // genData.applied_options shows which options were used
    // genData.svg_url and genData.pdf_url are cache links
}
```

## Workflow

### User Workflow

1. **Set Preferences (Once)**
   - User selects desired chemfig options (aromatic, show carbon, etc.)
   - Clicks "Save as Default Options"
   - Options are saved to server session

2. **Use Preferences (Repeatedly)**
   - User searches for any chemical
   - Server automatically applies saved options
   - Cache links are displayed
   - Preview shows structure with options applied

3. **Modify Preferences (As Needed)**
   - User can change options and save again
   - Or clear options to reset to defaults

### Technical Workflow

```
1. Client sends save options request
   “
2. Server stores in Flask session
   “
3. Client searches for chemical
   “
4. Server fetches SMILES from backend
   “
5. Server checks session for saved options
   “
6. Server generates chemfig WITH options
   “
7. Server caches result (SMILES + options hash)
   “
8. Server returns chemfig + cache links + applied_options
   “
9. Client displays results with cache links
```

## Cache Strategy

### Hash-Based Caching

Each unique combination of SMILES + options gets its own cache entry:

```python
def get_content_hash(smiles, options=None):
    """Generate unique hash for SMILES + options combination"""
    canonical = canonicalize_smiles(smiles)  # Prevent duplicates
    options_str = json.dumps(sorted(options or []))
    content = f"{canonical}:{options_str}"
    return hashlib.sha256(content.encode()).hexdigest()[:16]
```

**Examples:**
- Aspirin with no options: hash `abc123def456`
- Aspirin with `-o` (aromatic): hash `xyz789uvw012`
- Aspirin with `-o -m` (aromatic + methyl): hash `qwe456rty789`

Each gets its own cached SVG/PDF files that can be accessed via persistent URLs.

## Session Management

### Flask Sessions

- Uses Flask's built-in session management
- Stored server-side (secure)
- Persists across page reloads
- Expires when browser closes (or after configured timeout)

### Security Considerations

```python
# In production, set a strong secret key:
app.secret_key = os.environ.get('SECRET_KEY', 'fallback-key-for-dev')

# CORS must allow credentials:
CORS(app, supports_credentials=True)
```

## Testing the Implementation

### Manual Testing Steps

1. **Start the server:**
   ```bash
   python mol2chemfig_server.py
   ```

2. **Open test page:**
   ```
   http://localhost:5001/test_options_persistence.html
   ```

3. **Test options persistence:**
   - Select "Aromatic circles" and "Show carbon"
   - Click "Save as Default Options"
   - Search for "aspirin"
   - Verify options are applied
   - Search for "caffeine"
   - Verify same options are applied without re-selecting

4. **Test cache links:**
   - After generation, check that cache links appear
   - Click SVG and PDF links to verify they work
   - Refresh page and search again
   - Verify same cache links are returned (cached response)

5. **Test clear functionality:**
   - Click "Clear Options"
   - Search for "glucose"
   - Verify no options are applied

### Automated Testing

```python
import requests

SERVER = "http://localhost:5001"

# Test 1: Save options
response = requests.post(f"{SERVER}/api/options/save",
    json={"selections": ["-o", "-m"], "angle": 0, "indentation": 4, "h2": "keep"},
    cookies=session_cookies)
assert response.json()["success"] == True

# Test 2: Get saved options
response = requests.get(f"{SERVER}/api/options/get", cookies=session_cookies)
assert "-o" in response.json()["options"]["selections"]

# Test 3: Generate with auto-applied options
response = requests.post(f"{SERVER}/api/generate",
    json={"smiles": "CCO", "use_default_options": True},
    cookies=session_cookies)
data = response.json()
assert data["success"] == True
assert "-o" in data["applied_options"]
assert data["svg_url"] is not None
assert data["pdf_url"] is not None

# Test 4: Clear options
response = requests.post(f"{SERVER}/api/options/clear", cookies=session_cookies)
assert response.json()["success"] == True
```

## Integration with Existing Code

### test_m2cf_full.html Integration

To integrate with the existing test page:

```javascript
// At page load
document.addEventListener('DOMContentLoaded', () => {
    // ... existing code ...
    loadSavedOptionsFromServer();  // NEW: Load from server instead of localStorage
});

// Replace localStorage with server calls
async function loadSavedOptionsFromServer() {
    const response = await fetch('http://localhost:5001/api/options/get', {
        credentials: 'include'
    });
    const data = await response.json();
    if (data.success && data.options) {
        applyOptionsToUI(data.options.selections);
    }
}

// In applyMoleculeOptions function
async function applyMoleculeOptions() {
    const selections = getSelectedOptions('mol');

    // Save to server
    await fetch('http://localhost:5001/api/options/save', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        credentials: 'include',
        body: JSON.stringify({
            selections: selections,
            angle: parseInt(document.getElementById('moleculeAngle').value),
            indentation: parseInt(document.getElementById('moleculeIndent').value),
            h2: document.getElementById('moleculeH2').value
        })
    });

    // ... rest of apply logic ...
}
```

## Benefits

1. **User Experience:**
   - Options persist across searches
   - No need to manually reapply options
   - Cache links always displayed

2. **Performance:**
   - Each SMILES + options combination cached separately
   - No redundant generation
   - Fast retrieval from cache

3. **Reliability:**
   - Server-side session management
   - No localStorage issues
   - Works across all browsers

4. **Maintainability:**
   - Clean API design
   - Easy to extend with more options
   - Well-documented code

## Troubleshooting

### Issue: Options not persisting

**Solution:** Ensure `credentials: 'include'` is set in all fetch requests:
```javascript
fetch(url, {
    credentials: 'include',  // Required for session cookies
    // ... other options
})
```

### Issue: Cache links not showing

**Solution:** Check that response includes `svg_url` and `pdf_url`:
```python
# In server response
return jsonify({
    "svg_url": f"/images/{content_hash}.svg" if svg_path else None,
    "pdf_url": f"/images/{content_hash}.pdf" if pdf_path else None,
    # ...
})
```

### Issue: CORS errors

**Solution:** Ensure CORS is configured properly:
```python
CORS(app, supports_credentials=True)
```

## Future Enhancements

1. **Persistent Storage:**
   - Add database support for long-term option storage
   - User accounts for personalized defaults

2. **Option Presets:**
   - Allow users to save multiple option presets
   - Quick switching between presets

3. **Advanced Options:**
   - Support for angle, indentation, H2 treatment in defaults
   - Per-molecule option overrides

4. **Analytics:**
   - Track most-used options
   - Suggest optimal options based on molecule type

## Summary

The implementation successfully addresses all the original issues:

-  Chemfig options persist across searches
-  Cache links display properly after applying options
-  Options state managed server-side (session-based)
-  Users don't have to manually apply options each time
-  Clean API design with three new endpoints
-  Comprehensive test page demonstrating functionality
-  Full backward compatibility with existing code

The solution uses Flask sessions for server-side storage, automatic option application in the `/api/generate` endpoint, and a well-designed API for saving/loading/clearing preferences.
