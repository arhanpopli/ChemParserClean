"""
MolView Python Server
A simplified Python Flask server that mimics some MolView functionality for local use
"""

from flask import Flask, request, jsonify, send_from_directory, render_template_string, redirect, url_for
import os
import requests
import json
from urllib.parse import quote_plus
import tempfile


app = Flask(__name__,
           static_folder='molview-2.4.6 - github repo, src',
           static_url_path='')

# Configuration
MOLVIEW_DIR = os.path.join(os.path.dirname(__file__), 'molview-2.4.6 - github repo, src')

@app.route('/')
def home():
    """Serve the main MolView page"""
    # Since we can't run PHP through Flask, we'll serve an embedded interface
    # that can load MolView from our new PHP-based server
    html_content = '''
    <!DOCTYPE html>
    <html>
    <head>
        <title>MolView Local Server</title>
        <meta charset="UTF-8">
        <style>
            body { margin: 0; padding: 0; background: #1a1a2e; }
            #container { width: 100vw; height: 100vh; }
            iframe { width: 100%; height: 100%; border: none; }
        </style>
    </head>
    <body>
        <div id="container">
            <iframe src="http://localhost:5003/Molview/molview-2.4.6 - github repo, src/"></iframe>
        </div>
    </body>
    </html>
    '''
    return html_content

@app.route('/api/molview-convert')
def molview_convert():
    """
    API endpoint to convert molecular identifiers to different formats
    This mimics the CIR (Chemical Identifier Resolver) functionality that MolView uses
    """
    query = request.args.get('query', '')
    if not query:
        return jsonify({'error': 'Query parameter required'}), 400
    
    # Try to get SMILES from various sources
    try:
        # Try PubChem first
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote_plus(query)}/property/CanonicalSMILES/JSON"
        response = requests.get(pubchem_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                return jsonify({
                    'smiles': smiles,
                    'source': 'pubchem',
                    'name': query
                })
        
        # If PubChem fails, try other methods
        # For now, return error - in a full implementation, we'd try more sources
        return jsonify({'error': 'Could not resolve identifier', 'query': query}), 404
        
    except Exception as e:
        return jsonify({'error': str(e), 'query': query}), 500

@app.route('/api/get-structure/<identifier>')
def get_structure(identifier):
    """Get molecular structure data by identifier (CID, PDB ID, etc.)"""
    identifier_lower = identifier.lower()
    
    try:
        if identifier.startswith('cid:'):
            # PubChem CID
            cid = identifier[4:]  # Remove 'cid:' prefix
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
            response = requests.get(url)
            if response.status_code == 200:
                return response.text, 200, {'Content-Type': 'chemical/x-mdl-sdfile'}
        
        elif len(identifier) == 4 and identifier[0].isdigit():  # PDB ID format
            # PDB ID
            url = f"https://files.rcsb.org/download/{identifier}.cif"
            response = requests.get(url)
            if response.status_code == 200:
                return response.text, 200, {'Content-Type': 'chemical/x-cif'}
        
        else:
            # Treat as compound name
            # First get CID from name
            cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote_plus(identifier)}/cids/JSON"
            response = requests.get(cid_url)
            
            if response.status_code == 200:
                data = response.json()
                if data.get('IdentifierList', {}).get('CID'):
                    cid = str(data['IdentifierList']['CID'][0])
                    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
                    sdf_response = requests.get(sdf_url)
                    
                    if sdf_response.status_code == 200:
                        return sdf_response.text, 200, {'Content-Type': 'chemical/x-mdl-sdfile'}
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
    return jsonify({'error': 'Could not retrieve structure'}), 404

@app.route('/molview-proxy')
def molview_proxy():
    """Proxy to get around CORS issues when loading from external sources"""
    url = request.args.get('url')
    if not url:
        return jsonify({'error': 'URL parameter required'}), 400
    
    try:
        response = requests.get(url, timeout=15)
        return response.content, response.status_code, {'Content-Type': response.headers.get('Content-Type', 'application/octet-stream')}
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/molview-embedded')
def molview_embedded():
    """Embedded version of MolView that can be used in iframes"""
    # This would render a simplified version of MolView optimized for embedding
    query = request.args.get('q', '')
    smiles = request.args.get('smiles', '')
    cid = request.args.get('cid', '')
    
    # Create a minimal embeddable interface
    embed_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>MolView Embedded</title>
        <meta charset="UTF-8">
        <style>
            body {{ margin: 0; padding: 0; background: #1a1a2e; }}
            #container {{ width: 100vw; height: 100vh; }}
        </style>
    </head>
    <body>
        <div id="container">
            <h3 style="color: white; text-align: center; padding: 20px;">MolView Embedded (Placeholder)</h3>
            <p style="color: #aaa; text-align: center;">Query: {query or smiles or cid or 'None'}</p>
        </div>
        
        <!-- In a full implementation, we would load the actual MolView interface with these parameters -->
        <script>
            // This would initialize MolView with the provided parameters
            console.log('MolView embedded with parameters:', {{
                query: '{query}',
                smiles: '{smiles}',
                cid: '{cid}'
            }});
        </script>
    </body>
    </html>
    """
    return embed_html

if __name__ == '__main__':
    print("Starting MolView-compatible local server...")
    print("MolView will be available at: http://localhost:5004/")
    print("API endpoints available:")
    print("  - /api/molview-convert?query=[compound]")
    print("  - /api/get-structure/[identifier]")
    print("  - /molview-proxy?url=[url]")
    print("  - /molview-embedded?q=[query]&smiles=[smiles]&cid=[cid]")
    print()
    print("Press Ctrl+C to stop the server")
    app.run(host='localhost', port=5004, debug=False)