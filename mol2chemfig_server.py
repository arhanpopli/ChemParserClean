"""
Mol2ChemFig Server - Flask wrapper for mol2chemfig Docker backend
Provides persistent SVG/PDF links and Chrome extension integration

Similar to MoleculeViewer but uses mol2chemfig for superior rendering
Port: 5001 (MoleculeViewer uses 5000)
"""

from flask import Flask, request, jsonify, send_file, send_from_directory, session
from flask_cors import CORS
import requests
import os
import hashlib
import json
from datetime import datetime
from pathlib import Path
import sys

# Import canonicalization utility
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from canonicalize_smiles import canonicalize_smiles

app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'chemfig-server-secret-key-change-in-production')
CORS(app, supports_credentials=True)

# Configuration
MOL2CHEMFIG_BACKEND = "http://localhost:8000"  # Docker backend
STORAGE_DIR = Path("cache") / "mol2chemfig"  # Local storage for generated files - mol2chemfig specific
STORAGE_DIR.mkdir(parents=True, exist_ok=True)

# In-memory cache for quick lookups
image_cache = {}  # key: hash -> {svg: path, pdf: path, chemfig: str, options: [], timestamp}

# Default options storage (persistent across requests)
# This will store user preferences for chemfig options
default_options = {
    "selections": [],
    "angle": 0,
    "indentation": 4,
    "h2": "keep"
}

def get_content_hash(smiles, options=None):
    """
    Generate unique hash for SMILES + options combination
    Uses canonical SMILES to prevent duplicates
    """
    # Canonicalize SMILES to ensure same molecule gets same hash
    canonical = canonicalize_smiles(smiles)
    if canonical is None:
        # If canonicalization fails, use original SMILES
        print(f"Warning: Could not canonicalize SMILES '{smiles}', using original")
        canonical = smiles

    options_str = json.dumps(sorted(options or []))
    content = f"{canonical}:{options_str}"
    return hashlib.sha256(content.encode()).hexdigest()[:16]

def save_content(content, extension, content_hash):
    """Save content to disk and return relative path"""
    filename = f"{content_hash}.{extension}"
    filepath = STORAGE_DIR / filename
    
    if extension in ['svg', 'chemfig', 'mol']:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    else:  # binary files like PDF
        with open(filepath, 'wb') as f:
            f.write(content)
    
    return str(filepath)

@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint"""
    try:
        # Check if mol2chemfig backend is accessible
        response = requests.get(f"{MOL2CHEMFIG_BACKEND}/", timeout=2)
        backend_status = "healthy" if response.status_code == 200 else "unhealthy"
    except:
        backend_status = "unreachable"
    
    return jsonify({
        "status": "running",
        "server": "mol2chemfig_server",
        "port": 5001,
        "backend": backend_status,
        "storage": str(STORAGE_DIR),
        "cached_images": len(image_cache)
    })

@app.route('/api/generate', methods=['POST'])
def generate():
    """
    Generate molecule image from SMILES

    Request:
        {
            "smiles": "CCO",
            "format": "smiles",  # or "mol"
            "options": ["-o", "-m"],  # chemfig options (optional - will use saved defaults if not provided)
            "return_format": "svg",  # "svg", "pdf", or "both"
            "use_default_options": true  # if true, merge with saved defaults
        }

    Response:
        {
            "success": true,
            "hash": "abc123...",
            "svg_url": "/images/abc123.svg",
            "pdf_url": "/images/abc123.pdf",
            "chemfig": "\\chemfig{...}",
            "cached": false,
            "applied_options": ["-o", "-m"]  # shows which options were used
        }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles') or data.get('textAreaData', '')
        chem_format = data.get('format', 'smiles')
        options = data.get('options', None)
        return_format = data.get('return_format', 'svg')
        use_default_options = data.get('use_default_options', True)

        # If options not explicitly provided and use_default_options is true, use saved defaults
        if options is None and use_default_options:
            # Get default options from session or global defaults
            if 'default_options' in session:
                session_opts = session['default_options']
                options = session_opts.get('selections', [])
            else:
                options = default_options.get('selections', [])
        elif options is None:
            options = []
        
        if not smiles:
            return jsonify({"success": False, "error": "No SMILES/MOL data provided"}), 400
        
        # Generate hash for caching
        content_hash = get_content_hash(smiles, options)
        
        # Check cache first
        if content_hash in image_cache:
            cached = image_cache[content_hash]
            return jsonify({
                "success": True,
                "hash": content_hash,
                "svg_url": f"/images/{content_hash}.svg" if cached.get('svg') else None,
                "pdf_url": f"/images/{content_hash}.pdf" if cached.get('pdf') else None,
                "chemfig": cached.get('chemfig', ''),
                "cached": True,
                "applied_options": cached.get('options', []),
                "timestamp": cached.get('timestamp')
            })
        
        # Not in cache - generate via mol2chemfig backend
        if options:
            # With options - use /apply endpoint
            payload = {
                "chem_data": smiles,
                "chem_format": chem_format,
                "selections": options,
                "angle": 0,
                "indentation": 4,
                "h2": "keep"
            }
            response = requests.post(f"{MOL2CHEMFIG_BACKEND}/m2cf/apply", 
                                    json=payload, timeout=30)
        else:
            # No options - use /submit endpoint
            payload = {"textAreaData": smiles}
            response = requests.post(f"{MOL2CHEMFIG_BACKEND}/m2cf/submit", 
                                    json=payload, timeout=30)
        
        if response.status_code != 200:
            return jsonify({"success": False, "error": "Backend request failed"}), 500
        
        result = response.json()
        
        if result.get('error'):
            return jsonify({"success": False, "error": result['error']}), 400
        
        # Save generated content
        chemfig = result.get('chemfig', '')
        svg_path = None
        pdf_path = None
        
        # Fetch and save SVG if available
        if result.get('svglink') and return_format in ['svg', 'both']:
            svg_response = requests.get(f"{MOL2CHEMFIG_BACKEND}{result['svglink']}")
            if svg_response.status_code == 200:
                svg_path = save_content(svg_response.text, 'svg', content_hash)
        
        # Fetch and save PDF if available
        if result.get('pdflink') and return_format in ['pdf', 'both']:
            pdf_response = requests.get(f"{MOL2CHEMFIG_BACKEND}{result['pdflink']}")
            if pdf_response.status_code == 200:
                pdf_path = save_content(pdf_response.content, 'pdf', content_hash)
        
        # Cache the result
        image_cache[content_hash] = {
            "svg": svg_path,
            "pdf": pdf_path,
            "chemfig": chemfig,
            "options": options,
            "timestamp": datetime.now().isoformat()
        }
        
        return jsonify({
            "success": True,
            "hash": content_hash,
            "svg_url": f"/images/{content_hash}.svg" if svg_path else None,
            "pdf_url": f"/images/{content_hash}.pdf" if pdf_path else None,
            "chemfig": chemfig,
            "cached": False,
            "applied_options": options,
            "timestamp": datetime.now().isoformat()
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/layers', methods=['POST'])
def generate_layers():
    """
    Generate layered SVGs for a molecule
    
    Request:
        {
            "smiles": "CCO",
            "format": "smiles"
        }
    
    Response:
        {
            "success": true,
            "hash": "abc123...",
            "layers": {
                "base": "/images/abc123_base.svg",
                "aromatic": "/images/abc123_aromatic.svg",
                "carbon": "/images/abc123_carbon.svg",
                ...
            }
        }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles', '')
        chem_format = data.get('format', 'smiles')
        
        if not smiles:
            return jsonify({"success": False, "error": "No SMILES data provided"}), 400
        
        # Generate base hash
        content_hash = get_content_hash(smiles, ['layers'])
        
        # Call backend /layers endpoint
        payload = {
            "chem_data": smiles,
            "chem_format": chem_format
        }
        response = requests.post(f"{MOL2CHEMFIG_BACKEND}/m2cf/layers", 
                                json=payload, timeout=30)
        
        if response.status_code != 200:
            return jsonify({"success": False, "error": "Backend request failed"}), 500
        
        result = response.json()
        
        if result.get('error'):
            return jsonify({"success": False, "error": result['error']}), 400
        
        # Save each layer
        layer_urls = {}
        layers_data = result.get('layers', {})
        
        for layer_name, layer_svg in layers_data.items():
            if layer_name.endswith('_chemfig'):
                continue  # Skip chemfig data, only save SVGs
            
            if layer_svg and isinstance(layer_svg, str) and layer_svg.startswith('<?xml'):
                layer_hash = f"{content_hash}_{layer_name}"
                layer_path = save_content(layer_svg, 'svg', layer_hash)
                layer_urls[layer_name] = f"/images/{layer_hash}.svg"
        
        return jsonify({
            "success": True,
            "hash": content_hash,
            "layers": layer_urls,
            "timestamp": datetime.now().isoformat()
        })
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/search', methods=['GET', 'POST'])
def search():
    """
    Search for molecule by common name
    
    Request:
        GET: ?name=aspirin
        POST: {"name": "aspirin"}
    
    Response:
        {
            "success": true,
            "name": "aspirin",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "svg_url": "/images/abc123.svg"
        }
    """
    try:
        if request.method == 'GET':
            search_term = request.args.get('name') or request.args.get('searchTerm', '')
        else:
            data = request.get_json()
            search_term = data.get('name') or data.get('searchTerm', '')
        
        if not search_term:
            return jsonify({"success": False, "error": "No search term provided"}), 400
        
        # Call backend search endpoint
        response = requests.get(f"{MOL2CHEMFIG_BACKEND}/m2cf/search", 
                               params={"searchTerm": search_term}, timeout=30)
        
        if response.status_code != 200:
            return jsonify({"success": False, "error": "Backend search failed"}), 500
        
        result = response.json()
        
        if result.get('error'):
            return jsonify({"success": False, "error": result['error']}), 400
        
        # Generate image for the found molecule
        smiles = result.get('smiles', '')
        if smiles:
            gen_result = generate_from_data(smiles, 'smiles', [])
            return jsonify({
                "success": True,
                "name": search_term,
                "smiles": smiles,
                "svg_url": gen_result.get('svg_url'),
                "pdf_url": gen_result.get('pdf_url'),
                "chemfig": gen_result.get('chemfig')
            })
        
        return jsonify({"success": False, "error": "No SMILES found for search term"}), 404
        
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

def generate_from_data(smiles, format_type, options):
    """Helper function to generate molecule (used internally)"""
    content_hash = get_content_hash(smiles, options)
    
    # Check cache
    if content_hash in image_cache:
        cached = image_cache[content_hash]
        return {
            "svg_url": f"/images/{content_hash}.svg" if cached.get('svg') else None,
            "pdf_url": f"/images/{content_hash}.pdf" if cached.get('pdf') else None,
            "chemfig": cached.get('chemfig', '')
        }
    
    # Generate new
    payload = {"textAreaData": smiles}
    response = requests.post(f"{MOL2CHEMFIG_BACKEND}/m2cf/submit", json=payload, timeout=30)
    result = response.json()
    
    if result.get('svglink'):
        svg_response = requests.get(f"{MOL2CHEMFIG_BACKEND}{result['svglink']}")
        svg_path = save_content(svg_response.text, 'svg', content_hash)
        
        image_cache[content_hash] = {
            "svg": svg_path,
            "chemfig": result.get('chemfig', ''),
            "timestamp": datetime.now().isoformat()
        }
        
        return {
            "svg_url": f"/images/{content_hash}.svg",
            "chemfig": result.get('chemfig', '')
        }
    
    return {}

@app.route('/images/<filename>', methods=['GET'])
def serve_image(filename):
    """Serve generated SVG/PDF files"""
    try:
        filepath = STORAGE_DIR / filename
        if not filepath.exists():
            return jsonify({"error": "File not found"}), 404
        
        # Determine MIME type
        if filename.endswith('.svg'):
            mimetype = 'image/svg+xml'
        elif filename.endswith('.pdf'):
            mimetype = 'application/pdf'
        else:
            mimetype = 'application/octet-stream'
        
        return send_file(filepath, mimetype=mimetype)
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/api/cache/clear', methods=['POST'])
def clear_cache():
    """Clear the image cache (admin endpoint)"""
    try:
        image_cache.clear()
        return jsonify({
            "success": True,
            "message": "Cache cleared",
            "timestamp": datetime.now().isoformat()
        })
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/cache/stats', methods=['GET'])
def cache_stats():
    """Get cache statistics"""
    total_files = len(list(STORAGE_DIR.glob('*')))
    total_size = sum(f.stat().st_size for f in STORAGE_DIR.glob('*'))

    return jsonify({
        "cached_entries": len(image_cache),
        "total_files": total_files,
        "storage_size_mb": round(total_size / (1024 * 1024), 2),
        "storage_path": str(STORAGE_DIR)
    })

@app.route('/api/options/save', methods=['POST'])
def save_default_options():
    """
    Save default chemfig options for the user session

    Request:
        {
            "selections": ["-o", "-m"],
            "angle": 0,
            "indentation": 4,
            "h2": "keep"
        }

    Response:
        {
            "success": true,
            "message": "Default options saved",
            "options": {...}
        }
    """
    try:
        data = request.get_json()

        # Validate the options structure
        saved_opts = {
            "selections": data.get('selections', []),
            "angle": data.get('angle', 0),
            "indentation": data.get('indentation', 4),
            "h2": data.get('h2', 'keep')
        }

        # Save to session
        session['default_options'] = saved_opts
        session.modified = True

        return jsonify({
            "success": True,
            "message": "Default options saved successfully",
            "options": saved_opts,
            "timestamp": datetime.now().isoformat()
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/options/get', methods=['GET'])
def get_default_options():
    """
    Get currently saved default chemfig options

    Response:
        {
            "success": true,
            "options": {
                "selections": ["-o", "-m"],
                "angle": 0,
                "indentation": 4,
                "h2": "keep"
            }
        }
    """
    try:
        # Get from session or return global defaults
        if 'default_options' in session:
            opts = session['default_options']
        else:
            opts = default_options

        return jsonify({
            "success": True,
            "options": opts,
            "timestamp": datetime.now().isoformat()
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/options/clear', methods=['POST'])
def clear_default_options():
    """
    Clear saved default options (reset to defaults)

    Response:
        {
            "success": true,
            "message": "Options cleared"
        }
    """
    try:
        if 'default_options' in session:
            session.pop('default_options')
            session.modified = True

        return jsonify({
            "success": True,
            "message": "Default options cleared",
            "timestamp": datetime.now().isoformat()
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/opsin', methods=['GET', 'POST'])
def opsin_conversion():
    """
    Convert nomenclature to 3D SMILES using OPSIN

    Request:
        GET: ?name=glucose
        POST: {"name": "glucose"}

    Response:
        {
            "success": true,
            "name": "glucose",
            "smiles": "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
            "source": "OPSIN"
        }
    """
    try:
        if request.method == 'GET':
            name = request.args.get('name', '')
        else:
            data = request.get_json()
            name = data.get('name', '')

        if not name:
            return jsonify({"success": False, "error": "No name provided"}), 400

        # Call OPSIN API for 3D SMILES
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        response = requests.get(opsin_url, timeout=10)

        if response.status_code != 200:
            return jsonify({"success": False, "error": "OPSIN conversion failed"}), 500

        result = response.json()

        if result.get('error') or not result.get('smiles'):
            return jsonify({"success": False, "error": "No SMILES found"}), 404

        return jsonify({
            "success": True,
            "name": name,
            "smiles": result.get('smiles'),
            "source": "OPSIN"
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/generate-3d', methods=['POST'])
def generate_3d():
    """
    Generate molecule image from nomenclature with 3D stereochemistry via OPSIN

    Request:
        {
            "name": "glucose",
            "options": ["-o", "-m"],
            "return_format": "svg"
        }

    Response:
        {
            "success": true,
            "name": "glucose",
            "smiles": "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
            "hash": "abc123...",
            "svg_url": "/images/abc123.svg",
            "pdf_url": "/images/abc123.pdf",
            "chemfig": "\\chemfig{...}",
            "source": "OPSIN"
        }
    """
    try:
        data = request.get_json()
        name = data.get('name', '')
        options = data.get('options', [])
        return_format = data.get('return_format', 'svg')

        if not name:
            return jsonify({"success": False, "error": "No name provided"}), 400

        # Step 1: Convert name to 3D SMILES via OPSIN
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        opsin_response = requests.get(opsin_url, timeout=10)

        if opsin_response.status_code != 200:
            return jsonify({"success": False, "error": "OPSIN conversion failed"}), 500

        opsin_result = opsin_response.json()

        if opsin_result.get('error') or not opsin_result.get('smiles'):
            return jsonify({"success": False, "error": "No SMILES found via OPSIN"}), 404

        smiles = opsin_result.get('smiles')

        # Step 2: Generate image using mol2chemfig
        content_hash = get_content_hash(smiles, options)

        # Check cache first
        if content_hash in image_cache:
            cached = image_cache[content_hash]
            return jsonify({
                "success": True,
                "name": name,
                "smiles": smiles,
                "hash": content_hash,
                "svg_url": f"/images/{content_hash}.svg" if cached.get('svg') else None,
                "pdf_url": f"/images/{content_hash}.pdf" if cached.get('pdf') else None,
                "chemfig": cached.get('chemfig', ''),
                "cached": True,
                "source": "OPSIN",
                "timestamp": cached.get('timestamp')
            })

        # Generate via mol2chemfig backend
        if options:
            payload = {
                "chem_data": smiles,
                "chem_format": "smiles",
                "selections": options,
                "angle": 0,
                "indentation": 4,
                "h2": "keep"
            }
            backend_response = requests.post(f"{MOL2CHEMFIG_BACKEND}/m2cf/apply",
                                           json=payload, timeout=30)
        else:
            payload = {"textAreaData": smiles}
            backend_response = requests.post(f"{MOL2CHEMFIG_BACKEND}/m2cf/submit",
                                           json=payload, timeout=30)

        if backend_response.status_code != 200:
            return jsonify({"success": False, "error": "Backend request failed"}), 500

        backend_result = backend_response.json()

        if backend_result.get('error'):
            return jsonify({"success": False, "error": backend_result['error']}), 400

        # Save generated content
        chemfig = backend_result.get('chemfig', '')
        svg_path = None
        pdf_path = None

        # Fetch and save SVG if available
        if backend_result.get('svglink') and return_format in ['svg', 'both']:
            svg_response = requests.get(f"{MOL2CHEMFIG_BACKEND}{backend_result['svglink']}")
            if svg_response.status_code == 200:
                svg_path = save_content(svg_response.text, 'svg', content_hash)

        # Fetch and save PDF if available
        if backend_result.get('pdflink') and return_format in ['pdf', 'both']:
            pdf_response = requests.get(f"{MOL2CHEMFIG_BACKEND}{backend_result['pdflink']}")
            if pdf_response.status_code == 200:
                pdf_path = save_content(pdf_response.content, 'pdf', content_hash)

        # Cache the result
        image_cache[content_hash] = {
            "svg": svg_path,
            "pdf": pdf_path,
            "chemfig": chemfig,
            "options": options,
            "timestamp": datetime.now().isoformat()
        }

        return jsonify({
            "success": True,
            "name": name,
            "smiles": smiles,
            "hash": content_hash,
            "svg_url": f"/images/{content_hash}.svg" if svg_path else None,
            "pdf_url": f"/images/{content_hash}.pdf" if pdf_path else None,
            "chemfig": chemfig,
            "cached": False,
            "source": "OPSIN",
            "timestamp": datetime.now().isoformat()
        })

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/', methods=['GET'])
def index():
    """Serve the test HTML file"""
    try:
        html_file = Path("test_m2cf_full.html")
        if html_file.exists():
            with open(html_file, 'r', encoding='utf-8') as f:
                return f.read(), 200, {'Content-Type': 'text/html'}
        else:
            return jsonify({
                "name": "Mol2ChemFig Server",
                "version": "1.0.0",
                "description": "Flask server wrapper for mol2chemfig Docker backend",
                "port": 5001,
                "endpoints": {
                    "health": "/health",
                    "generate": "/api/generate (POST)",
                    "layers": "/api/layers (POST)",
                    "search": "/api/search (GET/POST)",
                    "opsin": "/api/opsin (GET/POST)",
                    "generate_3d": "/api/generate-3d (POST)",
                    "images": "/images/<filename>",
                    "cache_stats": "/api/cache/stats",
                    "cache_clear": "/api/cache/clear (POST)",
                    "options_save": "/api/options/save (POST)",
                    "options_get": "/api/options/get (GET)",
                    "options_clear": "/api/options/clear (POST)"
                },
                "backend": MOL2CHEMFIG_BACKEND,
                "note": "test_m2cf_full.html not found in current directory"
            })
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    print("=" * 70)
    print("🚀 Mol2ChemFig Server Starting...")
    print("=" * 70)
    print(f"Port: 5001")
    print(f"Backend: {MOL2CHEMFIG_BACKEND}")
    print(f"Storage: {STORAGE_DIR.absolute()}")
    print("=" * 70)
    print("\nEndpoints:")
    print("  GET  /                  - Server info")
    print("  GET  /health            - Health check")
    print("  POST /api/generate      - Generate molecule image")
    print("  POST /api/layers        - Generate layered SVGs")
    print("  GET  /api/search        - Search by name")
    print("  GET  /api/opsin         - Convert name to 3D SMILES via OPSIN")
    print("  POST /api/generate-3d   - Generate 3D molecule with OPSIN")
    print("  GET  /images/<file>     - Serve generated images")
    print("  POST /api/options/save  - Save default chemfig options")
    print("  GET  /api/options/get   - Get saved options")
    print("  POST /api/options/clear - Clear saved options")
    print("=" * 70)
    
    app.run(host='0.0.0.0', port=5001, debug=True)
