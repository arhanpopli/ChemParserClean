"""
PubChem Integration Server
This server provides endpoints for fetching images and 3D models from PubChem.
Port: 5002
"""

from flask import Flask, jsonify, send_file, request, Response, redirect, send_from_directory
from flask_cors import CORS
import requests
import os
import hashlib
import json
from datetime import datetime
import io

app = Flask(__name__)
CORS(app)

# Static files directory for 3D viewer
STATIC_DIR = os.path.join(os.path.dirname(__file__), 'MoleculeViewer', 'pubchem', 'static')
if not os.path.exists(STATIC_DIR):
    # Fallback to current directory
    STATIC_DIR = os.path.join(os.path.dirname(__file__), 'static')
os.makedirs(STATIC_DIR, exist_ok=True)

# Cache directory for PubChem data
CACHE_DIR = os.path.join(os.path.dirname(__file__), 'pubchem-cache')
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(os.path.join(CACHE_DIR, 'images'), exist_ok=True)
os.makedirs(os.path.join(CACHE_DIR, 'sdf'), exist_ok=True)

PORT = 5002

# PubChem API endpoints
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

print('=' * 70)
print('🚀 PubChem Integration Server Starting...')
print(f'📁 Cache directory: {CACHE_DIR}')
print('=' * 70)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_cache_key(text):
    """Generate cache key from text"""
    return hashlib.md5(text.encode('utf-8')).hexdigest()


def get_compound_cid(name):
    """
    Get PubChem CID (Compound ID) from chemical name or SMILES
    Returns: CID as integer or None if not found
    """
    try:
        # Try by name first
        url = f"{PUBCHEM_BASE}/compound/name/{name}/cids/JSON"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                cids = data['IdentifierList']['CID']
                return cids[0] if cids else None

        # If name lookup fails, try SMILES
        url = f"{PUBCHEM_BASE}/compound/smiles/{name}/cids/JSON"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                cids = data['IdentifierList']['CID']
                return cids[0] if cids else None

        return None
    except Exception as e:
        print(f"❌ Error getting CID for {name}: {e}")
        return None


def fetch_pubchem_image(cid, image_size='large', record_type='2d'):
    """
    Fetch PNG image from PubChem
    Args:
        cid: Compound ID
        image_size: 'small' (100x100), 'large' (300x300), or custom like '500x500'
        record_type: '2d' or '3d'
    """
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/PNG"
        params = {}
        if image_size:
            params['image_size'] = image_size
        if record_type:
            params['record_type'] = record_type

        response = requests.get(url, params=params, timeout=15)

        if response.status_code == 200:
            return response.content
        return None
    except Exception as e:
        print(f"❌ Error fetching image for CID {cid}: {e}")
        return None


def fetch_pubchem_sdf(cid, record_type='3d'):
    """
    Fetch SDF (Structure Data File) from PubChem for 3D models
    Args:
        cid: Compound ID
        record_type: '3d' or '2d'
    """
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/record/SDF"
        params = {'record_type': record_type}

        response = requests.get(url, params=params, timeout=15)

        if response.status_code == 200:
            return response.text
        return None
    except Exception as e:
        print(f"❌ Error fetching SDF for CID {cid}: {e}")
        return None


def get_conformers(cid):
    """Get list of 3D conformer IDs for a compound"""
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/conformers/TXT"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            conformer_ids = response.text.strip().split('\n')
            return [int(conf_id) for conf_id in conformer_ids if conf_id.strip()]
        return []
    except Exception as e:
        print(f"❌ Error getting conformers for CID {cid}: {e}")
        return []


# ============================================================
# ROUTES - PubChem Image API
# ============================================================

@app.route('/')
def index():
    """API information page"""
    return jsonify({
        'name': 'PubChem Integration Server',
        'version': '1.0.0',
        'port': PORT,
        'endpoints': {
            'image_by_name': 'GET /pubchem/image?name=histamine&size=large',
            'image_by_cid': 'GET /pubchem/image?cid=774&size=large',
            'image_direct': 'GET /pubchem/img/histamine (returns PNG directly)',
            '3d_model': 'GET /pubchem/3d-model?name=histamine&format=sdf',
            '3d_viewer': 'GET /pubchem/3d-viewer?name=histamine',
            'compound_info': 'GET /pubchem/info?name=histamine'
        },
        'examples': {
            'image': 'http://localhost:5002/pubchem/img/histamine',
            '3d_model': 'http://localhost:5002/pubchem/3d-model?name=histamine',
            '3d_viewer': 'http://localhost:5002/pubchem/3d-viewer?name=histamine'
        }
    })


@app.route('/health')
def health():
    """Health check endpoint"""
    cache_files = len([f for f in os.listdir(os.path.join(CACHE_DIR, 'images')) if f.endswith('.png')])
    return jsonify({
        'status': 'ok',
        'uptime': (datetime.now() - app.config.get('start_time', datetime.now())).total_seconds() if 'start_time' in app.config else 0,
        'timestamp': datetime.now().isoformat(),
        'cached_cids': cache_files
    })

@app.route('/static/<path:filename>')
def serve_static(filename):
    """Serve static files (viewer HTML, JS, etc.)"""
    return send_from_directory(STATIC_DIR, filename)


@app.route('/pubchem/image', methods=['GET'])
def get_image():
    """
    Get 2D or 3D image from PubChem
    Query params:
        - name: chemical name (e.g., 'histamine')
        - cid: PubChem compound ID (optional, faster if known)
        - size: 'small', 'large', or custom like '500x500' (default: 'large')
        - type: '2d' or '3d' (default: '2d')
    """
    name = request.args.get('name', '').strip()
    cid = request.args.get('cid', '').strip()
    size = request.args.get('size', 'large')
    record_type = request.args.get('type', '2d')

    if not name and not cid:
        return jsonify({'error': 'Missing required parameter: name or cid'}), 400

    print(f"\n{'=' * 70}")
    print(f"📥 [PubChem] GET /pubchem/image")
    print(f"   Name: {name if name else 'N/A'}")
    print(f"   CID: {cid if cid else 'N/A'}")
    print(f"   Size: {size}")
    print(f"   Type: {record_type}")

    # Get CID if not provided
    if not cid:
        cid = get_compound_cid(name)
        if not cid:
            print(f"❌ Compound not found: {name}")
            return jsonify({'error': f'Compound not found: {name}'}), 404
        print(f"   ✓ Found CID: {cid}")

    # Check cache
    cache_key = get_cache_key(f"img_{cid}_{size}_{record_type}")
    cache_file = os.path.join(CACHE_DIR, 'images', f"{cache_key}.png")

    if os.path.exists(cache_file):
        print(f"✅ Served from cache: {cache_file}")
        return jsonify({
            'cid': int(cid),
            'name': name,
            'image_url': f"http://localhost:{PORT}/pubchem/img/{name}?size={size}&type={record_type}",
            'direct_url': f"{PUBCHEM_BASE}/compound/cid/{cid}/PNG?image_size={size}&record_type={record_type}",
            'cached': True
        })

    # Fetch from PubChem
    image_data = fetch_pubchem_image(cid, size, record_type)

    if not image_data:
        print(f"❌ Failed to fetch image for CID {cid}")
        return jsonify({'error': 'Failed to fetch image from PubChem'}), 500

    # Cache the image
    with open(cache_file, 'wb') as f:
        f.write(image_data)

    print(f"✅ Image fetched and cached: {cache_file}")
    print('=' * 70 + '\n')

    return jsonify({
        'cid': int(cid),
        'name': name,
        'image_url': f"http://localhost:{PORT}/pubchem/img/{name}?size={size}&type={record_type}",
        'direct_url': f"{PUBCHEM_BASE}/compound/cid/{cid}/PNG?image_size={size}&record_type={record_type}",
        'cached': True
    })


@app.route('/pubchem/img/<name>')
def get_image_direct(name):
    """
    Direct PNG image endpoint (like MoleculeViewer)
    Returns the actual PNG image data
    Usage: <img src="http://localhost:5002/pubchem/img/histamine">
    """
    size = request.args.get('size', 'large')
    record_type = request.args.get('type', '2d')

    print(f"\n{'=' * 70}")
    print(f"📥 [PubChem] GET /pubchem/img/{name}")
    print(f"   Size: {size}, Type: {record_type}")

    # Get CID
    cid = get_compound_cid(name)
    if not cid:
        print(f"❌ Compound not found: {name}")
        # Return error image
        return Response(
            create_error_png(f"Not found: {name}"),
            mimetype='image/png'
        )

    print(f"   ✓ Found CID: {cid}")

    # Check cache
    cache_key = get_cache_key(f"img_{cid}_{size}_{record_type}")
    cache_file = os.path.join(CACHE_DIR, 'images', f"{cache_key}.png")

    if os.path.exists(cache_file):
        print(f"✅ Served from cache")
        print('=' * 70 + '\n')
        return send_file(cache_file, mimetype='image/png')

    # Fetch from PubChem
    image_data = fetch_pubchem_image(cid, size, record_type)

    if not image_data:
        print(f"❌ Failed to fetch image")
        print('=' * 70 + '\n')
        return Response(
            create_error_png(f"Failed to fetch: {name}"),
            mimetype='image/png'
        )

    # Cache the image
    with open(cache_file, 'wb') as f:
        f.write(image_data)

    print(f"✅ Image fetched and cached")
    print('=' * 70 + '\n')

    return Response(image_data, mimetype='image/png')


# ============================================================
# ROUTES - PubChem 3D Model API
# ============================================================

@app.route('/pubchem/3d-model', methods=['GET'])
def get_3d_model():
    """
    Get 3D model data from PubChem
    Query params:
        - name: chemical name (e.g., 'histamine')
        - cid: PubChem compound ID (optional)
        - format: 'sdf' or 'json' (default: 'sdf')
        - conformer: conformer index (0 for default, 1-N for alternatives)
    """
    name = request.args.get('name', '').strip()
    cid = request.args.get('cid', '').strip()
    format_type = request.args.get('format', 'sdf')
    conformer_idx = int(request.args.get('conformer', 0))

    if not name and not cid:
        return jsonify({'error': 'Missing required parameter: name or cid'}), 400

    print(f"\n{'=' * 70}")
    print(f"📥 [PubChem] GET /pubchem/3d-model")
    print(f"   Name: {name if name else 'N/A'}")
    print(f"   Format: {format_type}")
    print(f"   Conformer: {conformer_idx}")

    # Get CID if not provided
    if not cid:
        cid = get_compound_cid(name)
        if not cid:
            print(f"❌ Compound not found: {name}")
            return jsonify({'error': f'Compound not found: {name}'}), 404
        print(f"   ✓ Found CID: {cid}")

    # Check cache
    cache_key = get_cache_key(f"sdf_{cid}_{conformer_idx}")
    cache_file = os.path.join(CACHE_DIR, 'sdf', f"{cache_key}.sdf")

    sdf_data = None
    if os.path.exists(cache_file):
        print(f"✅ Served from cache")
        with open(cache_file, 'r') as f:
            sdf_data = f.read()
    else:
        # Fetch from PubChem
        if conformer_idx == 0:
            # Get default conformer
            sdf_data = fetch_pubchem_sdf(cid, '3d')
        else:
            # Get specific conformer
            conformers = get_conformers(cid)
            if conformer_idx <= len(conformers):
                conformer_id = conformers[conformer_idx - 1]
                try:
                    url = f"{PUBCHEM_BASE}/conformers/{conformer_id}/SDF"
                    response = requests.get(url, timeout=15)
                    if response.status_code == 200:
                        sdf_data = response.text
                except Exception as e:
                    print(f"❌ Error fetching conformer: {e}")

        if sdf_data:
            # Cache the SDF
            with open(cache_file, 'w') as f:
                f.write(sdf_data)
            print(f"✅ 3D model fetched and cached")

    if not sdf_data:
        print(f"❌ Failed to fetch 3D model")
        print('=' * 70 + '\n')
        return jsonify({'error': 'Failed to fetch 3D model from PubChem'}), 500

    print('=' * 70 + '\n')

    if format_type == 'json':
        return jsonify({
            'cid': int(cid),
            'name': name,
            'format': 'sdf',
            'data': sdf_data,
            'viewer_url': f"http://localhost:{PORT}/pubchem/3d-viewer?cid={cid}",
            'pubchem_url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer"
        })
    else:
        # Return SDF file directly
        return Response(
            sdf_data,
            mimetype='chemical/x-mdl-sdfile',
            headers={'Content-Disposition': f'attachment; filename="{name or cid}_3d.sdf"'}
        )


@app.route('/pubchem/3d-viewer', methods=['GET'])
def get_3d_viewer():
    """
    Interactive 3D molecular viewer using PubChem data
    Query params:
        - name: chemical name
        - cid: PubChem compound ID
        - embed: 'true' to return embeddable HTML, 'false' to redirect
    """
    name = request.args.get('name', '').strip()
    cid = request.args.get('cid', '').strip()
    embed = request.args.get('embed', 'false').lower() == 'true'

    # Get CID if not provided
    if not cid:
        cid = get_compound_cid(name)
        if not cid:
            return jsonify({'error': f'Compound not found: {name}'}), 404

    pubchem_3d_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer"

    if embed:
        # Return HTML with interactive 3D viewer
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>3D Viewer - {name or cid}</title>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <style>
                * {{
                    margin: 0;
                    padding: 0;
                    box-sizing: border-box;
                }}
                body {{
                    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: #333;
                    overflow: hidden;
                }}
                .container {{
                    height: 100vh;
                    display: flex;
                    flex-direction: column;
                }}
                .header {{
                    background: rgba(255, 255, 255, 0.95);
                    padding: 15px 20px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                    z-index: 10;
                }}
                .header h1 {{
                    font-size: 24px;
                    color: #667eea;
                    margin-bottom: 5px;
                }}
                .compound-info {{
                    font-size: 14px;
                    color: #666;
                }}
                .controls {{
                    background: rgba(255, 255, 255, 0.95);
                    padding: 15px 20px;
                    display: flex;
                    gap: 10px;
                    flex-wrap: wrap;
                    box-shadow: 0 -2px 10px rgba(0,0,0,0.1);
                    z-index: 10;
                }}
                .btn {{
                    padding: 8px 16px;
                    border: none;
                    border-radius: 6px;
                    font-size: 13px;
                    font-weight: 600;
                    cursor: pointer;
                    transition: all 0.3s ease;
                    text-decoration: none;
                    display: inline-block;
                }}
                .btn-primary {{
                    background: #667eea;
                    color: white;
                }}
                .btn-primary:hover {{
                    background: #5568d3;
                    transform: translateY(-2px);
                    box-shadow: 0 4px 12px rgba(102, 126, 234, 0.4);
                }}
                .btn-secondary {{
                    background: #764ba2;
                    color: white;
                }}
                .btn-secondary:hover {{
                    background: #63408a;
                    transform: translateY(-2px);
                    box-shadow: 0 4px 12px rgba(118, 75, 162, 0.4);
                }}
                .btn-success {{
                    background: #10b981;
                    color: white;
                }}
                .btn-success:hover {{
                    background: #059669;
                }}
                .viewer-container {{
                    flex: 1;
                    background: white;
                    margin: 20px;
                    border-radius: 12px;
                    box-shadow: 0 10px 40px rgba(0,0,0,0.2);
                    overflow: hidden;
                    position: relative;
                }}
                .viewer-frame {{
                    width: 100%;
                    height: 100%;
                    border: none;
                }}
                .loading {{
                    position: absolute;
                    top: 50%;
                    left: 50%;
                    transform: translate(-50%, -50%);
                    text-align: center;
                    color: #667eea;
                }}
                .loading-spinner {{
                    width: 50px;
                    height: 50px;
                    border: 5px solid #e0e7ff;
                    border-top-color: #667eea;
                    border-radius: 50%;
                    animation: spin 1s linear infinite;
                    margin: 0 auto 20px;
                }}
                @keyframes spin {{
                    to {{ transform: rotate(360deg); }}
                }}
                .option-group {{
                    display: flex;
                    align-items: center;
                    gap: 8px;
                    padding: 5px 12px;
                    background: #f3f4f6;
                    border-radius: 6px;
                }}
                .option-group label {{
                    font-size: 13px;
                    color: #4b5563;
                    font-weight: 500;
                }}
                .option-group select {{
                    padding: 5px 10px;
                    border: 1px solid #d1d5db;
                    border-radius: 4px;
                    font-size: 13px;
                    background: white;
                    cursor: pointer;
                }}
                .checkbox-group {{
                    display: flex;
                    align-items: center;
                    gap: 6px;
                    padding: 5px 12px;
                    background: #f3f4f6;
                    border-radius: 6px;
                    cursor: pointer;
                }}
                .checkbox-group input[type="checkbox"] {{
                    cursor: pointer;
                }}
                .checkbox-group label {{
                    font-size: 13px;
                    color: #4b5563;
                    font-weight: 500;
                    cursor: pointer;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1>🔬 3D Molecular Viewer</h1>
                    <div class="compound-info">
                        <strong>{name or 'Unknown'}</strong> • PubChem CID: {cid}
                    </div>
                </div>
                
                <div class="viewer-container">
                    <div class="loading" id="loading">
                        <div class="loading-spinner"></div>
                        <p>Loading 3D structure...</p>
                    </div>
                    <iframe id="viewer" class="viewer-frame" src="{pubchem_3d_url}" style="display:none"></iframe>
                </div>
                
                <div class="controls">
                    <div class="option-group">
                        <label>Style:</label>
                        <select id="styleSelect">
                            <option value="ball-stick">Ball & Stick</option>
                            <option value="stick">Stick</option>
                            <option value="sphere">Space Filling</option>
                            <option value="wireframe">Wireframe</option>
                        </select>
                    </div>
                    
                    <div class="checkbox-group">
                        <input type="checkbox" id="showHydrogens" checked>
                        <label for="showHydrogens">Show Hydrogens</label>
                    </div>
                    
                    <div class="checkbox-group">
                        <input type="checkbox" id="rotate">
                        <label for="rotate">Auto Rotate</label>
                    </div>
                    
                    <a href="{pubchem_3d_url}" target="_blank" class="btn btn-primary">
                        🔗 Open in PubChem
                    </a>
                    
                    <a href="http://localhost:{PORT}/pubchem/3d-model?cid={cid}" class="btn btn-secondary">
                        💾 Download SDF
                    </a>
                    
                    <button class="btn btn-success" onclick="window.close()">
                        ✓ Close
                    </button>
                </div>
            </div>
            
            <script>
                // Load iframe and hide loading spinner
                const viewer = document.getElementById('viewer');
                const loading = document.getElementById('loading');
                
                viewer.onload = function() {{
                    loading.style.display = 'none';
                    viewer.style.display = 'block';
                }};
                
                // Control handlers (these would communicate with the iframe if we had access)
                // For now, they're placeholders showing the intended functionality
                document.getElementById('styleSelect').addEventListener('change', function(e) {{
                    console.log('Style changed to:', e.target.value);
                    // In a full implementation, this would send a message to the PubChem viewer
                }});
                
                document.getElementById('showHydrogens').addEventListener('change', function(e) {{
                    console.log('Show hydrogens:', e.target.checked);
                }});
                
                document.getElementById('rotate').addEventListener('change', function(e) {{
                    console.log('Auto rotate:', e.target.checked);
                }});
            </script>
        </body>
        </html>
        """
        return Response(html, mimetype='text/html')
    else:
        # Redirect to PubChem
        return redirect(pubchem_3d_url)


# ============================================================
# ROUTES - Compound Information
# ============================================================

@app.route('/pubchem/info', methods=['GET'])
def get_compound_info():
    """
    Get compound information from PubChem
    Query params:
        - name: chemical name
        - cid: PubChem compound ID
    """
    name = request.args.get('name', '').strip()
    cid = request.args.get('cid', '').strip()

    if not name and not cid:
        return jsonify({'error': 'Missing required parameter: name or cid'}), 400

    # Get CID if not provided
    if not cid:
        cid = get_compound_cid(name)
        if not cid:
            return jsonify({'error': f'Compound not found: {name}'}), 404

    # Get conformers info
    conformers = get_conformers(cid)

    return jsonify({
        'cid': int(cid),
        'name': name,
        'pubchem_url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        '3d_viewer_url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer",
        'image_url_2d': f"{PUBCHEM_BASE}/compound/cid/{cid}/PNG?image_size=large",
        'image_url_3d': f"{PUBCHEM_BASE}/compound/cid/{cid}/PNG?image_size=large&record_type=3d",
        'sdf_url': f"{PUBCHEM_BASE}/compound/cid/{cid}/record/SDF?record_type=3d",
        'conformers_count': len(conformers),
        'local_endpoints': {
            'image': f"http://localhost:{PORT}/pubchem/img/{name or cid}",
            '3d_model': f"http://localhost:{PORT}/pubchem/3d-model?cid={cid}",
            '3d_viewer': f"http://localhost:{PORT}/pubchem/3d-viewer?cid={cid}"
        }
    })


# ============================================================
# CACHE MANAGEMENT
# ============================================================

@app.route('/cache-info', methods=['GET'])
def cache_info():
    """Get cache statistics"""
    images = os.listdir(os.path.join(CACHE_DIR, 'images'))
    sdfs = os.listdir(os.path.join(CACHE_DIR, 'sdf'))

    return jsonify({
        'cache_directory': CACHE_DIR,
        'cached_images': len(images),
        'cached_3d_models': len(sdfs),
        'total_files': len(images) + len(sdfs)
    })


@app.route('/clear-cache', methods=['DELETE'])
def clear_cache():
    """Clear all cached data"""
    try:
        images_dir = os.path.join(CACHE_DIR, 'images')
        sdf_dir = os.path.join(CACHE_DIR, 'sdf')

        count = 0
        for file in os.listdir(images_dir):
            os.remove(os.path.join(images_dir, file))
            count += 1
        for file in os.listdir(sdf_dir):
            os.remove(os.path.join(sdf_dir, file))
            count += 1

        return jsonify({
            'message': f'Cleared {count} cached files',
            'status': 'success'
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def create_error_png(message):
    """Create a simple error PNG (placeholder - returns empty PNG)"""
    # This is a minimal 1x1 transparent PNG
    png_data = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x06\x00\x00\x00\x1f\x15\xc4\x89\x00\x00\x00\nIDATx\x9cc\x00\x01\x00\x00\x05\x00\x01\r\n-\xb4\x00\x00\x00\x00IEND\xaeB`\x82'
    return png_data


# ============================================================
# START SERVER
# ============================================================

if __name__ == '__main__':
    print('\n' + '=' * 70)
    print(f'✅ PubChem Server running on http://localhost:{PORT}')
    print('=' * 70)
    print('\n📍 API Endpoints:')
    print(f'   Image (direct):  http://localhost:{PORT}/pubchem/img/histamine')
    print(f'   Image (info):    http://localhost:{PORT}/pubchem/image?name=histamine')
    print(f'   3D Model:        http://localhost:{PORT}/pubchem/3d-model?name=histamine')
    print(f'   3D Viewer:       http://localhost:{PORT}/pubchem/3d-viewer?name=histamine')
    print(f'   Compound Info:   http://localhost:{PORT}/pubchem/info?name=histamine')
    print(f'\n💾 Cache Directory: {CACHE_DIR}')
    print('=' * 70 + '\n')

    app.run(host='0.0.0.0', port=PORT, debug=True)
