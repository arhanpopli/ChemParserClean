"""
MoleculeViewer API - SMILES to SVG converter with nomenclature lookup.
Simple focused app for molecule visualization.
"""

from flask import Flask, jsonify, request, render_template, send_file, Response
from flask_cors import CORS
import json
import os
import hashlib
import time
from datetime import datetime, timedelta
from pathlib import Path
from app.chemistry import smiles_to_svg, nomenclature_to_smiles, get_molecule_info

# Load environment variables
try:
    from dotenv import load_dotenv
    load_dotenv()
except:
    pass

# Get the templates folder path relative to this file
template_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'templates'))
app = Flask(__name__, template_folder=template_dir)
CORS(app)

# ============================================
# IMAGE CACHE SYSTEM (24-hour expiration)
# ============================================
CACHE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'svg-cache'))
CACHE_EXPIRY_HOURS = 24

# Worldwide accessible URL (set via environment or request)
PUBLIC_BASE_URL = os.environ.get('PUBLIC_BASE_URL', 'http://localhost:5000')

# Create cache directory if it doesn't exist
os.makedirs(CACHE_DIR, exist_ok=True)

def get_cache_hash(content):
    """Generate hash for cache key"""
    return hashlib.md5(content.encode()).hexdigest()

def cleanup_old_cache():
    """Remove cache files older than 24 hours"""
    try:
        now = time.time()
        for filename in os.listdir(CACHE_DIR):
            filepath = os.path.join(CACHE_DIR, filename)
            if os.path.isfile(filepath):
                age_hours = (now - os.path.getmtime(filepath)) / 3600
                if age_hours > CACHE_EXPIRY_HOURS:
                    os.remove(filepath)
                    print(f"Cleaned old cache: {filename}")
    except Exception as e:
        print(f"Cache cleanup error: {e}")

def save_to_cache(svg_content, identifier):
    """Save SVG to cache and return URL and local path"""
    try:
        cache_hash = get_cache_hash(svg_content)
        filename = f"{identifier}_{cache_hash}.svg"
        filepath = os.path.join(CACHE_DIR, filename)
        
        # Only write if not already cached
        if not os.path.exists(filepath):
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(svg_content)
            print(f"Cached SVG: {filename}")
        
        # Return cache URL using public base URL (worldwide accessible)
        cache_url = f"{PUBLIC_BASE_URL}/cache/{filename}"
        return cache_url, filepath
    except Exception as e:
        print(f"Cache save error: {e}")
        return None, None


@app.route('/', methods=['GET'])
def home():
    """Serve the main HTML interface."""
    return render_template('index.html')


@app.route('/api/smiles-to-svg', methods=['POST'])
def smiles_to_svg_endpoint():
    """
    Convert SMILES to SVG with visualization options matching mol2chemfig Docker API.
    
    Request:
        {
            "smiles": "C1=CC=CC=C1",
            "width": 600,
            "height": 500,
            "options": {
                "show_carbons": false,        # Display carbon atom symbols
                "show_methyls": false,        # Show methyl group symbols (CH3)
                "aromatic_circles": true,     # Draw circles in aromatic rings
                "fancy_bonds": true,          # Fancy double/triple bond rendering
                "atom_numbers": false,        # Show atom indices
                "hydrogens": "keep",          # 'keep', 'add', or 'delete'
                "flip_horizontal": false,     # Flip X axis
                "flip_vertical": false,       # Flip Y axis
                "rotate": 0,                  # Rotation angle in degrees
                "recalculate_coordinates": false  # Recalculate 2D coordinates
            }
        }
    
    Response:
        {
            "error": null,
            "svg": "<svg>...</svg>",
            "smiles": "C1=CC=CC=C1",
            "info": {
                "molecular_weight": 78.11,
                "formula": "C6H6",
                ...
            }
        }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()
        width = data.get('width', 600)
        height = data.get('height', 500)
        options = data.get('options', {})
        
        # LOG INCOMING REQUEST
        print(f"\n{'='*70}")
        print(f"[MoleculeViewer API] Receiving SMILES request")
        print(f"{'='*70}")
        print(f"  SMILES: {smiles}")
        print(f"  Size: {width}x{height}px")
        print(f"  Options: {json.dumps(options, indent=2)}")
        print(f"{'='*70}\n")
        
        if not smiles:
            return jsonify({
                'error': 'SMILES required',
                'svg': None
            }), 400
        
        error, svg = smiles_to_svg(smiles, width, height, options)
        
        if error:
            print(f"[MoleculeViewer] Error converting SMILES: {error}\n")
            return jsonify({
                'error': error,
                'svg': None
            }), 400
        
        # Get molecule info
        error_info, info = get_molecule_info(smiles)
        
        print(f"[MoleculeViewer] Successfully rendered: {smiles}")
        print(f"   Formula: {info.get('formula', 'N/A')}")
        print(f"   Molecular Weight: {info.get('molecular_weight', 'N/A')}\n")
        
        return jsonify({
            'error': None,
            'svg': svg,
            'smiles': smiles,
            'info': info
        }), 200
        
    except Exception as e:
        print(f"\n[MoleculeViewer] SERVER ERROR: {str(e)}\n")
        return jsonify({
            'error': 'Server error: {}'.format(str(e)),
            'svg': None
        }), 500


@app.route('/api/nomenclature-to-smiles', methods=['POST'])
def nomenclature_to_smiles_endpoint():
    """
    Convert chemical name/nomenclature to SMILES with source information.
    
    Request:
        {
            "nomenclature": "aspirin"
        }
    
    Response:
        {
            "error": null,
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "nomenclature": "aspirin",
            "source": "ChemDoodle Database | OPSIN Parser | Fallback Dictionary | PubChem (CID: 2244)"
        }
    """
    try:
        data = request.get_json()
        nomenclature = data.get('nomenclature', '').strip()
        
        if not nomenclature:
            return jsonify({
                'error': 'Nomenclature required',
                'smiles': None,
                'source': None
            }), 400
        
        error, smiles, source = nomenclature_to_smiles(nomenclature)
        
        if smiles is None:
            return jsonify({
                'error': error,
                'smiles': None,
                'nomenclature': nomenclature,
                'source': None
            }), 404
        
        return jsonify({
            'error': None,
            'smiles': smiles,
            'nomenclature': nomenclature,
            'source': source
        }), 200
        
    except Exception as e:
        return jsonify({
            'error': 'Server error: {}'.format(str(e)),
            'smiles': None,
            'source': None
        }), 500


@app.route('/api/nomenclature-to-svg', methods=['POST'])
def nomenclature_to_svg_endpoint():
    """
    Convert chemical name directly to SVG with visualization options.
    
    Request:
        {
            "nomenclature": "benzene",
            "width": 600,
            "height": 500,
            "options": {
                "show_carbons": false,        # Display carbon atom symbols
                "show_methyls": false,        # Show methyl group symbols (CH3)
                "aromatic_circles": true,     # Draw circles in aromatic rings
                "fancy_bonds": true,          # Fancy double/triple bond rendering
                "atom_numbers": false,        # Show atom indices
                "hydrogens": "keep",          # 'keep', 'add', or 'delete'
                "flip_horizontal": false,     # Flip X axis
                "flip_vertical": false,       # Flip Y axis
                "rotate": 0,                  # Rotation angle in degrees
                "recalculate_coordinates": false  # Recalculate 2D coordinates
            }
        }
    
    Response:
        {
            "error": null,
            "svg": "<svg>...</svg>",
            "smiles": "C1=CC=CC=C1",
            "nomenclature": "benzene",
            "info": {...}
        }
    """
    try:
        data = request.get_json()
        nomenclature = data.get('nomenclature', '').strip()
        width = data.get('width', 600)
        height = data.get('height', 500)
        options = data.get('options', {})
        
        # LOG INCOMING NOMENCLATURE REQUEST
        print(f"\n{'='*70}")
        print(f"[MoleculeViewer API] Receiving nomenclature request")
        print(f"{'='*70}")
        print(f"  Nomenclature: {nomenclature}")
        print(f"  Size: {width}x{height}px")
        print(f"  Options: {json.dumps(options, indent=2)}")
        print(f"{'='*70}\n")
        
        if not nomenclature:
            return jsonify({
                'error': 'Nomenclature required',
                'svg': None
            }), 400
        
        # Step 1: Get SMILES from nomenclature
        error, smiles, source = nomenclature_to_smiles(nomenclature)
        
        if smiles is None:
            print(f"[MoleculeViewer] Nomenclature lookup failed: {error}\n")
            return jsonify({
                'error': error,
                'svg': None,
                'smiles': None,
                'nomenclature': nomenclature,
                'source': None
            }), 404
        
        print(f"Nomenclature -> SMILES: {nomenclature} -> {smiles}")
        print(f"  Source: {source}\n")
        
        # Step 2: Convert SMILES to SVG
        error, svg = smiles_to_svg(smiles, width, height, options)
        
        if error:
            print(f"[MoleculeViewer] Error converting SMILES: {error}\n")
            return jsonify({
                'error': error,
                'svg': None,
                'smiles': smiles,
                'nomenclature': nomenclature,
                'source': source
            }), 400
        
        # Get molecule info
        error_info, info = get_molecule_info(smiles)
        
        print(f"[MoleculeViewer] Successfully rendered: {nomenclature}")
        print(f"   SMILES: {smiles}")
        print(f"   Formula: {info.get('formula', 'N/A')}")
        print(f"   Molecular Weight: {info.get('molecular_weight', 'N/A')}\n")
        
        return jsonify({
            'error': None,
            'svg': svg,
            'smiles': smiles,
            'nomenclature': nomenclature,
            'source': source,
            'info': info
        }), 200
        
    except Exception as e:
        print(f"\n[MoleculeViewer] SERVER ERROR: {str(e)}\n")
        return jsonify({
            'error': 'Server error: {}'.format(str(e)),
            'svg': None
        }), 500


@app.route('/api/molecule-info', methods=['POST'])
def molecule_info_endpoint():
    """
    Get molecular information (weight, formula, etc.) from SMILES.
    
    Request:
        {
            "smiles": "C1=CC=CC=C1"
        }
    
    Response:
        {
            "error": null,
            "info": {
                "molecular_weight": 78.11,
                "formula": "C6H6",
                "num_atoms": 6,
                ...
            }
        }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()
        
        if not smiles:
            return jsonify({
                'error': 'SMILES required',
                'info': None
            }), 400
        
        error, info = get_molecule_info(smiles)
        
        if error:
            return jsonify({
                'error': error,
                'info': None
            }), 400
        
        return jsonify({
            'error': None,
            'info': info
        }), 200
        
    except Exception as e:
        return jsonify({
            'error': 'Server error: {}'.format(str(e)),
            'info': None
        }), 500


@app.route('/api/render-smiles', methods=['GET'])
def render_smiles_image():
    """
    Render SMILES as SVG image - optimized for Chrome extension.
    
    Query Parameters:
        smiles: SMILES string (required)
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        dark: Dark mode (true/false, default: false)
    
    Returns:
        SVG image with proper Content-Type header
    
    Example:
        GET /api/render-smiles?smiles=CCO
        GET /api/render-smiles?smiles=c1ccccc1&width=400&height=300
    """
    try:
        smiles = request.args.get('smiles', '').strip()
        width = int(request.args.get('width', 300))
        height = int(request.args.get('height', 200))
        dark_mode = request.args.get('dark', 'false').lower() == 'true'
        
        if not smiles:
            return jsonify({'error': 'SMILES required'}), 400
        
        # Call chemistry engine to render SVG
        error, svg = smiles_to_svg(
            smiles, 
            width, 
            height,
            options={
                'dark_mode': dark_mode,
                'aromatic_circles': True,
                'fancy_bonds': True
            }
        )
        
        if error:
            # Return error as SVG text for browser rendering
            error_svg = f'''<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
                <rect width="{width}" height="{height}" fill="{'#222' if dark_mode else '#fff'}"/>
                <text x="10" y="30" fill="{'#f00' if dark_mode else '#d00'}" font-family="monospace" font-size="12">Error: {error}</text>
                <text x="10" y="50" fill="{'#aaa' if dark_mode else '#999'}" font-family="monospace" font-size="10">SMILES: {smiles[:50]}</text>
            </svg>'''
            return error_svg, 200, {'Content-Type': 'image/svg+xml'}
        
        # Return SVG with proper headers for browser caching and CORS
        return svg, 200, {
            'Content-Type': 'image/svg+xml',
            'Cache-Control': 'public, max-age=86400',  # Cache for 1 day
            'Access-Control-Allow-Origin': '*',  # Allow any origin for Chrome extension
            'Access-Control-Allow-Methods': 'GET, OPTIONS',
        }
        
    except Exception as e:
        error_svg = f'''<svg width="300" height="200" xmlns="http://www.w3.org/2000/svg">
            <rect width="300" height="200" fill="#fff"/>
            <text x="10" y="30" fill="#d00" font-family="monospace" font-size="12">Server Error</text>
            <text x="10" y="50" fill="#999" font-family="monospace" font-size="10">{str(e)[:50]}</text>
        </svg>'''
        return error_svg, 500, {'Content-Type': 'image/svg+xml'}


@app.route('/api/render-nomenclature', methods=['GET'])
def render_nomenclature_image():
    """
    Render chemical nomenclature as SVG image - like CodeCogs service.
    
    Query Parameters:
        nomenclature: Chemical name (required)
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        dark: Dark mode (true/false, default: false)
    
    Returns:
        SVG image with proper Content-Type header
    
    Example:
        GET /api/render-nomenclature?nomenclature=acetone
        GET /api/render-nomenclature?nomenclature=benzene&width=400&height=300
    """
    try:
        nomenclature = request.args.get('nomenclature', '').strip()
        width = int(request.args.get('width', 300))
        height = int(request.args.get('height', 200))
        dark_mode = request.args.get('dark', 'false').lower() == 'true'
        
        if not nomenclature:
            return jsonify({'error': 'Nomenclature required'}), 400
        
        # Step 1: Convert nomenclature to SMILES
        error, smiles, source = nomenclature_to_smiles(nomenclature)
        
        if smiles is None:
            # Return error as SVG text for browser rendering
            error_svg = f'''<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
                <rect width="{width}" height="{height}" fill="{'#222' if dark_mode else '#fff'}"/>
                <text x="10" y="30" fill="{'#f00' if dark_mode else '#d00'}" font-family="monospace" font-size="12">Not Found</text>
                <text x="10" y="50" fill="{'#aaa' if dark_mode else '#999'}" font-family="monospace" font-size="10">{nomenclature}</text>
            </svg>'''
            return error_svg, 200, {'Content-Type': 'image/svg+xml'}
        
        # Step 2: Render SMILES to SVG
        error, svg = smiles_to_svg(
            smiles, 
            width, 
            height,
            options={
                'dark_mode': dark_mode,
                'aromatic_circles': True,
                'fancy_bonds': True
            }
        )
        
        if error:
            # Return error as SVG text for browser rendering
            error_svg = f'''<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
                <rect width="{width}" height="{height}" fill="{'#222' if dark_mode else '#fff'}"/>
                <text x="10" y="30" fill="{'#f00' if dark_mode else '#d00'}" font-family="monospace" font-size="12">Render Error</text>
                <text x="10" y="50" fill="{'#aaa' if dark_mode else '#999'}" font-family="monospace" font-size="10">SMILES: {smiles[:30]}</text>
            </svg>'''
            return error_svg, 200, {'Content-Type': 'image/svg+xml'}
        
        print(f"✅ [MoleculeViewer] Rendered via URL: {nomenclature} → {smiles}")
        
        # Return SVG with proper headers for browser caching and CORS
        return svg, 200, {
            'Content-Type': 'image/svg+xml',
            'Cache-Control': 'public, max-age=86400',  # Cache for 1 day
            'Access-Control-Allow-Origin': '*',  # Allow any origin for Chrome extension
            'Access-Control-Allow-Methods': 'GET, OPTIONS',
        }
        
    except Exception as e:
        error_svg = f'''<svg width="300" height="200" xmlns="http://www.w3.org/2000/svg">
            <rect width="300" height="200" fill="#fff"/>
            <text x="10" y="30" fill="#d00" font-family="monospace" font-size="12">Server Error</text>
            <text x="10" y="50" fill="#999" font-family="monospace" font-size="10">{str(e)[:50]}</text>
        </svg>'''
        return error_svg, 500, {'Content-Type': 'image/svg+xml'}


@app.route('/img/smiles', methods=['GET'])
def img_smiles():
    """Direct SMILES to SVG image endpoint with caching."""
    try:
        smiles = request.args.get('smiles', '').strip()
        width = request.args.get('width', '300')
        height = request.args.get('height', '200')
        return_json = request.args.get('json', 'true').lower() == 'true'
        
        if not smiles:
            return jsonify({'error': 'SMILES required'}) if return_json else ("SMILES required", 400)
        
        try:
            width = int(width)
            height = int(height)
        except:
            width, height = 300, 200
        
        error, svg = smiles_to_svg(smiles, width, height, {})
        if error:
            return jsonify({'error': error})
        
        # Cache the SVG
        cache_url, filepath = save_to_cache(svg, f"smiles_{smiles[:10]}")
        
        print(f"SMILES endpoint: {smiles}")
        print(f"Cache URL: {cache_url}")
        
        # Always return JSON with cache link and SVG
        return jsonify({
            'success': True,
            'smiles': smiles,
            'cache_url': cache_url,
            'image_url': cache_url,  # Direct download link
            'expires_in_hours': CACHE_EXPIRY_HOURS,
            'svg': svg  # Also include SVG for direct display
        }), 200
    except Exception as e:
        print(f"❌ Error in /img/smiles: {str(e)}")
        return jsonify({'error': str(e), 'success': False})


@app.route('/img/nomenclature', methods=['GET'])
def img_nomenclature():
    """Direct nomenclature to SVG image endpoint with caching."""
    try:
        nomenclature = request.args.get('nomenclature', '').strip()
        width = request.args.get('width', '300')
        height = request.args.get('height', '200')
        return_json = request.args.get('json', 'true').lower() == 'true'
        
        if not nomenclature:
            return jsonify({'error': 'Nomenclature required'}) if return_json else ("Nomenclature required", 400)
        
        try:
            width = int(width)
            height = int(height)
        except:
            width, height = 300, 200
        
        # Convert name to SMILES
        result = nomenclature_to_smiles(nomenclature)
        if isinstance(result, tuple):
            if len(result) == 3:
                error, smiles, source = result
            else:
                error, smiles = result
        else:
            return jsonify({'error': 'Invalid nomenclature response', 'success': False})
            
        if not smiles:
            error_msg = f"Cannot convert '{nomenclature}' to SMILES"
            return jsonify({'error': error_msg, 'success': False})
        
        # Render to SVG
        error, svg = smiles_to_svg(smiles, width, height, {})
        if error:
            return jsonify({'error': error, 'success': False})
        
        # Cache the SVG
        cache_url, filepath = save_to_cache(svg, f"nomenclature_{nomenclature[:10]}")
        
        print(f"Nomenclature endpoint: {nomenclature} -> {smiles}")
        print(f"Cache URL: {cache_url}")
        
        # Always return JSON with cache link and SVG
        return jsonify({
            'success': True,
            'nomenclature': nomenclature,
            'smiles': smiles,
            'cache_url': cache_url,
            'image_url': cache_url,  # Direct download link
            'expires_in_hours': CACHE_EXPIRY_HOURS,
            'svg': svg  # Also include SVG for direct display
        }), 200
    except Exception as e:
        print(f"Error in /img/nomenclature: {str(e)}")
        return jsonify({'error': str(e), 'success': False})


@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint."""
    return jsonify({'status': 'ok'}), 200


@app.route('/cache/<filename>', methods=['GET'])
def serve_cache(filename):
    """Serve cached SVG files"""
    try:
        # Security: prevent directory traversal
        if '..' in filename or '/' in filename:
            return "Invalid filename", 400
        
        filepath = os.path.join(CACHE_DIR, filename)
        if not os.path.exists(filepath):
            return "File not found", 404
        
        # Read and serve SVG directly with correct headers
        with open(filepath, 'r', encoding='utf-8') as f:
            svg_content = f.read()
        
        return Response(svg_content, mimetype='image/svg+xml', headers={
            'Content-Disposition': f'attachment; filename="{filename}"'
        })
    except Exception as e:
        return jsonify({'error': str(e), 'success': False})


@app.route('/cache/cleanup', methods=['POST'])
def cleanup_cache():
    """Manually trigger cache cleanup"""
    cleanup_old_cache()
    return jsonify({'status': 'cleaned'}), 200


if __name__ == '__main__':
    # Run cache cleanup on startup
    cleanup_old_cache()
    app.run(debug=False, host='0.0.0.0', port=5000, use_reloader=False)
