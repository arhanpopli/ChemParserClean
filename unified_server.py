import sys
import os
import hashlib
import json
import base64
from datetime import datetime
from pathlib import Path
from flask import Flask, request, jsonify, send_file, Response
from flask_cors import CORS
import requests

# Add MoleculeViewer to path so we can import its app
sys.path.insert(0, os.path.abspath('MoleculeViewer'))

# Import the existing MoleculeViewer app
# This gives us all the RDKit endpoints (/api/smiles-to-svg, /img/smiles, etc.)
try:
    from app.api import app
    print("✅ Imported MoleculeViewer app")
except ImportError as e:
    print(f"❌ Failed to import MoleculeViewer app: {e}")
    # Fallback if import fails
    app = Flask(__name__)
    CORS(app)

# Import native mol2chemfig
import native_mol2chemfig

# Configuration
STORAGE_DIR = Path("temp/images")
STORAGE_DIR.mkdir(parents=True, exist_ok=True)
SEARCH_SERVER = "http://localhost:8001" 

# Cache
image_cache = {}

# Helper Functions
def get_content_hash(text, selections, h2_option):
    """Generate a unique hash for the input content and options"""
    content = f"{text}|{json.dumps(selections)}|{h2_option}"
    return hashlib.md5(content.encode()).hexdigest()

def save_svg_content(svg_content, content_hash):
    """Save SVG content to disk and return the filename"""
    filename = f"{content_hash}.svg"
    filepath = STORAGE_DIR / filename
    
    # Ensure SVG has namespace
    if '<svg' in svg_content and 'xmlns=' not in svg_content:
        svg_content = svg_content.replace('<svg', '<svg xmlns="http://www.w3.org/2000/svg"')
        
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg_content)
    return filename

# =============================================================================
# MOL2CHEMFIG ENDPOINTS
# =============================================================================

@app.route('/m2cf/submit', methods=['POST'])
def m2cf_submit():
    """
    Process SMILES to chemfig + SVG using NATIVE mol2chemfig
    Replaces the Docker backend call.
    """
    try:
        data = request.get_json()
        text_area_data = data.get('textAreaData', '').strip()
        
        if not text_area_data:
            return jsonify({"error": "No input data provided"}), 400
            
        print(f"[Unified] Processing: {text_area_data[:50]}...")
        
        # Check cache
        h2_option = data.get('h2', 'keep')
        content_hash = get_content_hash(text_area_data, data.get('selections', []), h2_option)
        
        if content_hash in image_cache:
            print(f"[Unified] Serving from cache: {content_hash}")
            cached = image_cache[content_hash]
            return jsonify({
                "svglink": f"/images/{content_hash}.svg",
                "chemfig": cached.get('chemfig', ''),
                "cached": True
            })
            
        # Run Native mol2chemfig
        # Map options from request to native_mol2chemfig options
        options = {
            'm2cfHydrogensMode': h2_option,
            'm2cfShowCarbons': data.get('show_carbons', False),
            'm2cfAromaticCircles': data.get('aromatic_circles', False),
            'm2cfShowMethyls': data.get('show_methyls', False),
            'm2cfFancyBonds': data.get('fancy_bonds', False),
            'm2cfAtomNumbers': data.get('atom_numbers', False),
            'm2cfCompact': data.get('compact', False),
            'm2cfFlipHorizontal': data.get('flip_horizontal', False),
            'm2cfFlipVertical': data.get('flip_vertical', False),
            'm2cfRotate': data.get('rotate', 0)
        }
        
        svg_content = native_mol2chemfig.run_mol2chemfig(text_area_data, options)
        
        if not svg_content:
            return jsonify({"error": "Native rendering failed"}), 500
            
        # Save to cache
        filename = save_svg_content(svg_content, content_hash)
        
        # Update cache
        image_cache[content_hash] = {
            "svg": str(STORAGE_DIR / filename),
            "chemfig": "", # We don't get raw chemfig easily from the native run yet, but that's fine
            "timestamp": datetime.now().isoformat()
        }
        
        return jsonify({
            "svglink": f"/images/{filename}",
            "chemfig": "",
            "pdflink": None # We don't return PDF link in native mode yet, can add if needed
        })

    except Exception as e:
        print(f"[Unified] Error: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

@app.route('/images/<filename>', methods=['GET'])
def serve_image(filename):
    """Serve generated SVG files"""
    try:
        filepath = STORAGE_DIR / filename
        if not filepath.exists():
            return jsonify({"error": "File not found"}), 404
        return send_file(filepath, mimetype='image/svg+xml')
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/m2cf/search', methods=['GET'])
def m2cf_search():
    """
    Proxy search requests to the Universal Search API (Port 8001)
    """
    query = request.args.get('q', '').strip()
    if not query:
        return jsonify({"error": "No query provided"}), 400
        
    try:
        # Forward to Search Server (8001)
        resp = requests.get(f"{SEARCH_SERVER}/search", params={'q': query}, timeout=10)
        return Response(resp.content, resp.status_code, content_type=resp.headers['Content-Type'])
    except Exception as e:
        return jsonify({"error": f"Search failed: {e}"}), 500

@app.route('/health', methods=['GET'])
def health_check():
    return jsonify({
        "status": "running",
        "server": "unified_server",
        "port": 1000,
        "features": ["moleculeviewer", "mol2chemfig-native"]
    })

if __name__ == '__main__':
    print("🚀 Starting Unified Server on port 1000...")
    print("   • MoleculeViewer (RDKit) endpoints active")
    print("   • mol2chemfig (Native) endpoints active")
    app.run(host='0.0.0.0', port=1000, debug=True, use_reloader=False)
