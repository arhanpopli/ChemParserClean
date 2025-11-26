#!/usr/bin/env python3
"""
MolView Local Server
This server provides a local version of MolView to avoid CORS issues
when using it in the chem-extension.
"""

import os
import sys
import json
from flask import Flask, request, jsonify, send_from_directory, render_template_string, abort, redirect
import requests
from urllib.parse import urlparse, parse_qs
import re
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Get the directory containing this script
script_dir = os.path.dirname(os.path.abspath(__file__))
molview_static_dir = os.path.join(script_dir, "downloaded website UI")

app = Flask(__name__, static_folder=molview_static_dir, static_url_path='')

# Create a proper embed endpoint similar to embed.molview.org/v1
@app.route('/embed/')
def embed_handler():
    """
    Serve the embed version of MolView for clean integration
    Accepts query parameters like ?smiles=, ?cid=, ?pdbid=, etc.
    Similar to embed.molview.org/v1/
    """
    # Get query parameters
    smiles = request.args.get('smiles')
    cid = request.args.get('cid')
    pdbid = request.args.get('pdbid')
    codid = request.args.get('codid')
    query = request.args.get('q', request.args.get('query'))

    # Create embed page that loads the specified molecule
    embed_content = f'''
<!DOCTYPE html>
<html>

<head>
    <meta charset="UTF-8" />
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <meta name="viewport" content="width=device-width, user-scalable=no" />
    <title>MolView Embed</title>

    <!-- CSS -->
    <link type="text/css" rel="stylesheet" href="https://fonts.googleapis.com/css?family=Open+Sans:400,700" />
    <link type="text/css" rel="stylesheet" href="/css/molview-app.min.css" media="screen" />
    <link type="text/css" rel="stylesheet" href="/css/molview-desktop.min.css" media="screen" />

    <!-- JS -->
    <script type="text/javascript" src="/js/molview-base.min.js"></script>
    <script type="text/javascript" src="/js/molview-applib.min.js"></script>
    <script type="text/javascript" src="/js/molview-datasets.min.js"></script>
    <script type="text/javascript" src="/js/molview-core.min.js"></script>
    <script type="text/javascript" src="/js/molview-molpad.min.js"></script>
    <script type="text/javascript" src="/js/molview-app.min.js"></script>

    <!-- Data injection -->
    <script type="text/javascript">
        Model.JSmol.hq = true;
        MolView.touch = false;
        MolView.mobile = false;
        MolView.layout = "layout-model";
        MolView.embed = true; // This is an embed instance

        Request.CIR.available = true;

        // Load the molecule after initialization
        $(document).ready(function() {{
            setTimeout(function() {{
                // Load based on query parameters
                if("{smiles}") {{
                    Actions.load_smiles("{smiles}");
                }} else if("{cid}") {{
                    Actions.load_cid("{cid}");
                }} else if("{pdbid}") {{
                    Actions.load_pdbid("{pdbid}");
                }} else if("{codid}") {{
                    Actions.load_codid("{codid}");
                }} else if("{query}") {{
                    // Search for the query term
                    Search.execute("{query}", "all");
                }}
            }}, 1000);
        }});
    </script>
</head>
<body id="model" style="background:#000000; margin: 0; padding: 0; overflow: hidden;">
    <div id="progress" style="display: block; position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">
        <canvas id="progress-canvas"></canvas>
    </div>
    <div id="chemdoodle" class="render-engine full-cover" style="display: none;"><canvas id="chemdoodle-canvas"></canvas></div>
    <div id="jsmol" class="render-engine full-cover" style="display: none;"></div>
    <div id="glmol" class="render-engine full-cover" style="display: none;"></div>
    <div id="messages"></div>

    <!-- Hide all UI elements for clean embed -->
    <style>
        #menu, #search, #main-layer, #progress, .toolbar, .btn, .dropdown, #autocomplete-dropdown-wrapper, #sketcher {{
            display: none !important;
        }}
        #model {{
            margin: 0;
            padding: 0;
            width: 100vw;
            height: 100vh;
        }}
        body {{
            overflow: hidden;
        }}
    </style>
</body>
</html>
    '''
    return embed_content

# Serve the main index page
@app.route('/')
def index():
    return send_from_directory(molview_static_dir, 'index.html')

# API endpoints that MolView uses - these will proxy to external services
@app.route('/api/cod/cif/<codid>')
def api_cod_cif(codid):
    """Forward COD CIF requests"""
    try:
        url = f"https://www.crystallography.net/cod/cif/{codid}.cif"
        response = requests.get(url, timeout=10)
        return response.content, response.status_code, {'Content-Type': 'text/plain'}
    except Exception as e:
        logger.error(f"Error fetching COD CIF {codid}: {e}")
        return "Not found", 404

@app.route('/api/cod/<action>/<query>')
def api_cod(action, query):
    """Forward COD requests"""
    try:
        url = f"https://www.crystallography.net/cod/{action}/{query}"
        response = requests.get(url, timeout=10)
        return response.json() if response.headers.get('content-type', '').startswith('application/json') else response.text
    except Exception as e:
        logger.error(f"Error fetching COD {action}/{query}: {e}")
        return {"error": "Unable to fetch data"}, 500

@app.route('/api/nist/lookup/<cas>')
def api_nist_lookup(cas):
    """Forward NIST lookup requests"""
    try:
        # This would need to be implemented based on NIST API
        return {"message": "NIST lookup temporarily unavailable"}, 501
    except Exception as e:
        logger.error(f"Error fetching NIST lookup {cas}: {e}")
        return {"error": "Unable to fetch data"}, 500

@app.route('/api/image/<dbid>/<identifier>')
def api_image(dbid, identifier):
    """Forward image requests"""
    try:
        if dbid.lower() == "pubchem":
            # Redirect to PubChem image API
            redirect_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{identifier}/png"
            return redirect(redirect_url)
        return {"error": "Unknown database"}, 400
    except Exception as e:
        logger.error(f"Error fetching image {dbid}/{identifier}: {e}")
        return {"error": "Unable to fetch image"}, 500

# Page endpoints
@app.route('/<page>')
def page_handler(page):
    if page in ['readme', 'changelog', 'license', 'legal', 'internetExplorer', 'htmlCanvas', 'tracking']:
        # Return a simple page or template for these special pages
        return f"<h1>{page.title()} Page</h1><p>This is the {page} page.</p>"
    else:
        # For other requests, serve the main index
        return send_from_directory(molview_static_dir, 'index.html')

# Fallback for any other API requests
@app.route('/api/<path:subpath>')
def api_fallback(subpath):
    """Handle other API requests"""
    return {"status": "API endpoint not implemented", "path": subpath}, 501

# Generic page handler
@app.route('/page.php')
def page_php():
    error_id = request.args.get('id', 'unknown')
    return f"<h1>Error Page: {error_id}</h1><p>Error occurred at: {request.args}</p>"

# Serve static files (CSS, JS, images, etc.) - THIS MUST BE LAST
# This catch-all route should be defined last to avoid interfering with other routes
@app.route('/<path:filename>')
def serve_static(filename):
    return send_from_directory(molview_static_dir, filename)

if __name__ == '__main__':
    # Set up the server
    port = int(os.environ.get('MOLVIEW_PORT', 5003))
    logger.info(f"Starting MolView server on port {port}")
    logger.info(f"Serving from directory: {molview_static_dir}")
    
    # Run the server
    app.run(host='localhost', port=port, debug=False)