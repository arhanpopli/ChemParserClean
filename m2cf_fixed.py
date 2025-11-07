import tempfile
import os
import shutil
import json
from pprint import pprint
try:
    import urllib2 as urllib_request
    import urllib as urllib_parse
except ImportError:
    import urllib.request as urllib_request
    import urllib.parse as urllib_parse

from flask import Blueprint, jsonify, request
import pubchempy as pcp

from rdkit import Chem
from rdkit.Chem import AllChem

from chemistry.utils import (
    smiles_mol_to_chemfig,
    convert_mol_format,
    update_chemfig,
    get_smiles,
    combine_args,
    chemfig_to_pdf
)

from chemistry.kekule_parser import generate_latex, update_reaction_chemfig

# Import SVG generation functions
try:
    from rdkit.Chem import Draw
    HAS_SVG_SUPPORT = True
except ImportError:
    HAS_SVG_SUPPORT = False

m2cf = Blueprint('m2cf', __name__, url_prefix="/m2cf")

# Extract SVG differences - returns only the NEW elements added by an option
def extract_svg_diff(base_svg, option_svg):
    """
    Extract only the visual differences between base and option SVG.
    Returns an SVG with only the NEW elements (labels, circles, etc.) that the option adds.
    This creates true "layers" that can be stacked on top of the base.
    
    CRITICAL: Compares FULL element tags (including coordinates) to detect differences.
    """
    import re
    
    try:
        # Extract viewBox and dimensions from base SVG
        viewbox_match = re.search(r'viewBox=["\']([^"\']+)["\']', base_svg)
        width_match = re.search(r'width=["\']([^"\']+)["\']', base_svg)
        height_match = re.search(r'height=["\']([^"\']+)["\']', base_svg)
        
        viewbox = viewbox_match.group(1) if viewbox_match else "0 0 100 100"
        width = width_match.group(1) if width_match else "100pt"
        height = height_match.group(1) if height_match else "100pt"
        
        # Extract FULL element tags (including all attributes like x, y, transform, etc.)
        def extract_full_elements(svg):
            # Match complete path tags
            paths = set(re.findall(r'<path[^>]+/>', svg))
            # Match complete text tags (opening tag + content + closing tag)
            texts = set(re.findall(r'<text[^>]*>[^<]*</text>', svg))
            # Match complete circle tags
            circles = set(re.findall(r'<circle[^>]+/>', svg))
            return paths, texts, circles
        
        base_paths, base_texts, base_circles = extract_full_elements(base_svg)
        opt_paths, opt_texts, opt_circles = extract_full_elements(option_svg)
        
        # Find NEW elements (in option but not in base)
        # These are elements with different positions/attributes OR completely new content
        new_paths = opt_paths - base_paths
        new_texts = opt_texts - base_texts
        new_circles = opt_circles - base_circles
        
        print("DEBUG: Base has %d paths, %d texts, %d circles" % (len(base_paths), len(base_texts), len(base_circles)))
        print("DEBUG: Option has %d paths, %d texts, %d circles" % (len(opt_paths), len(opt_texts), len(opt_circles)))
        print("DEBUG: Diff has %d paths, %d texts, %d circles" % (len(new_paths), len(new_texts), len(new_circles)))
        
        # Build diff SVG with only new elements
        diff_svg_parts = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"',
            ' width="{}" height="{}" viewBox="{}">'.format(width, height, viewbox),
            '<defs>',
            '<style type="text/css"><![CDATA[',
            'text { font-family: "Computer Modern", serif; }',
            ']]></style>',
            '</defs>',
            '<g>'
        ]
        
        # Add new text elements (CH3, C labels, atom numbers, etc.)
        for text_element in new_texts:
            diff_svg_parts.append(text_element)
        
        # Add new circles (aromatic rings)
        for circle in new_circles:
            diff_svg_parts.append(circle)
        
        # Add new paths (fancy bonds, etc.)
        for path in new_paths:
            diff_svg_parts.append(path)
        
        diff_svg_parts.append('</g>')
        diff_svg_parts.append('</svg>')
        
        diff_result = '\n'.join(diff_svg_parts)
        print("DEBUG: Diff SVG size: %d chars" % len(diff_result))
        
        return diff_result
        
    except Exception as e:
        print("SVG diff extraction error:", e)
        import traceback
        traceback.print_exc()
        return option_svg  # Fallback to full option SVG

# SVG generation using pdflatex + dvisvgm --pdf (EXACT same rendering as PDF!)
def chemfig_to_svg(chemfig_code):
    """
    Convert chemfig LaTeX code to SVG using latex + dvisvgm.
    Uses the EXACT SAME LaTeX template as pdfgen.py with proper bond spacing settings.
    Key settings: \\setdoublesep{3pt}, \\setcrambond{2.5pt}{0.4pt}{1.0pt}
    """
    tempdir = tempfile.mkdtemp()
    tex_file = os.path.join(tempdir, 'molecule.tex')
    dvi_file = os.path.join(tempdir, 'molecule.dvi')
    svg_file = os.path.join(tempdir, 'molecule.svg')
    
    # Path to mol2chemfig.sty (same as pdfgen.py)
    m2pkg_path = '/usr/src/app/src/mol2chemfig'
    pkg = '/mol2chemfig.sty'
    
    # LaTeX template - EXACT SAME as original pdfgen.py (lines 128-145)
    # This includes critical bond spacing parameters!
    latex_content = r'''\documentclass{minimal}
\usepackage{xcolor, mol2chemfig}
\usepackage[helvet]{sfmath}
\setcrambond{2.5pt}{0.4pt}{1.0pt}
\setbondoffset{1pt}
\setdoublesep{3pt}
\setatomsep{16pt}
\renewcommand{\printatom}[1]{\fontsize{8pt}{10pt}\selectfont{\ensuremath{\mathsf{#1}}}}
\setlength{\parindent}{0pt}
\begin{document}
%s
\end{document}''' % chemfig_code
    
    try:
        # Create symlink to mol2chemfig.sty (same as pdfgen.py does)
        try:
            os.symlink(m2pkg_path + pkg, tempdir + pkg)
        except (OSError, AttributeError):
            # Fallback if symlink fails
            shutil.copy(m2pkg_path + pkg, tempdir + pkg)
        
        # Write LaTeX file
        with open(tex_file, 'w') as f:
            f.write(latex_content)
        
        # Compile using latex (not pdflatex - dvisvgm needs DVI)
        latex_cmd = 'latex -interaction=nonstopmode %s > /dev/null 2>&1' % tex_file
        curdir = os.getcwd()
        os.chdir(tempdir)
        os.system(latex_cmd)
        
        if not os.path.exists(dvi_file):
            os.chdir(curdir)
            return None, "latex compilation failed"
        
        # Convert DVI to SVG - extract vector with proper fonts
        dvisvgm_cmd = 'dvisvgm --font-format=woff --exact %s -o %s 2>&1' % (dvi_file, svg_file)
        result = os.system(dvisvgm_cmd)
        
        if not os.path.exists(svg_file):
            os.chdir(curdir)
            return None, "PDF to SVG conversion failed (exit code: %d)" % result
        
        # Read SVG file
        with open(svg_file, 'r') as f:
            svg_content = f.read()
        
        os.chdir(curdir)
        
        # Minimal post-processing - just ensure proper display
        if '<svg' in svg_content and '<defs>' in svg_content:
            # Insert minimal style for consistent browser rendering
            style_block = '''<style type="text/css">
    <![CDATA[
    text { font-family: "Computer Modern", serif; }
    ]]>
</style>'''
            svg_content = svg_content.replace('<defs>', '<defs>' + style_block)
        
        return svg_content, None
        
    except Exception as e:
        return None, str(e)
    finally:
        # Clean up temp directory
        shutil.rmtree(tempdir, ignore_errors=True)

# Common chemical name to SMILES lookup (for names PubChem doesn't recognize)
COMMON_NAMES = {
    'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'acetylsalicylic acid': 'CC(=O)Oc1ccccc1C(=O)O',
    'benzene': 'c1ccccc1',
    'toluene': 'Cc1ccccc1',
    'ethanol': 'CCO',
    'methane': 'C',
    'ethane': 'CC',
    'propane': 'CCC',
    'butane': 'CCCC',
    'acetone': 'CC(=O)C',
    'acetic acid': 'CC(=O)O',
    'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'glucose': 'C([C@@H]1[C@H]([C@@H]([C@H](C(=O)O1)O)O)O)O',
    'water': 'O',
}


def get_smiles_from_opsin(compound_name):
    """
    Use OPSIN JAR to convert IUPAC names to SMILES (fallback when PubChem fails)
    Returns (error, smiles)
    """
    try:
        import subprocess
        import os
        
        # Check if OPSIN JAR exists
        opsin_jar = "/usr/src/app/opsin-cli.jar"
        if not os.path.exists(opsin_jar):
            return "OPSIN not available", None
        
        # Call OPSIN to convert name to SMILES
        result = subprocess.check_output([
            "java", "-jar", opsin_jar,
            "-o", "smi",  # Output SMILES format
            compound_name
        ], stderr=subprocess.STDOUT, timeout=10)
        
        smiles = result.decode('utf-8').strip()
        if smiles and not smiles.startswith("Error"):
            return None, smiles
        else:
            return "OPSIN could not parse '{}'".format(compound_name), None
    except Exception as e:
        return "OPSIN error: {}".format(str(e)), None


def get_smiles_from_pubchem_api(compound_name):
    """
    Fetch SMILES from PubChem REST API using compound name
    Falls back to local lookup table, then OPSIN for names not in PubChem
    Returns (error, smiles)
    """
    try:
        # First check local common names table (case-insensitive)
        lower_name = compound_name.lower().strip()
        if lower_name in COMMON_NAMES:
            return None, COMMON_NAMES[lower_name]
        
        # URL encode the compound name to handle spaces and special characters
        try:
            from urllib import quote
        except ImportError:
            from urllib.parse import quote
        
        encoded_name = quote(compound_name.encode('utf-8'))
        
        # Try compound name search first
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES/JSON".format(encoded_name)
        try:
            response = urllib_request.urlopen(url, timeout=10)
            data = json.loads(response.read().decode('utf-8'))
            
            if 'properties' in data and len(data['properties']) > 0:
                smiles = data['properties'][0].get('CanonicalSMILES')
                if smiles:
                    return None, smiles
        except:
            pass
        
        # Fallback: Try searching by substance name (more flexible)
        url_substance = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/{}/cids/JSON?cids_type=same".format(encoded_name)
        try:
            response = urllib_request.urlopen(url_substance, timeout=10)
            data = json.loads(response.read().decode('utf-8'))
            
            if 'IdentifierList' in data and 'CID' in data['IdentifierList'] and len(data['IdentifierList']['CID']) > 0:
                cid = data['IdentifierList']['CID'][0]
                # Now get the SMILES for this CID
                url_cid = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/CanonicalSMILES/JSON".format(cid)
                response = urllib_request.urlopen(url_cid, timeout=10)
                cid_data = json.loads(response.read().decode('utf-8'))
                
                if 'properties' in cid_data and len(cid_data['properties']) > 0:
                    smiles = cid_data['properties'][0].get('CanonicalSMILES')
                    if smiles:
                        return None, smiles
        except:
            pass
        
        # Final fallback: Try OPSIN for IUPAC names
        error_opsin, smiles_opsin = get_smiles_from_opsin(compound_name)
        if smiles_opsin:
            return None, smiles_opsin
        
        return "Compound '{}' not found in PubChem, common names, or OPSIN databases".format(compound_name), None
    except Exception as e:
        return "Error searching PubChem: {}".format(str(e)), None


@m2cf.route('/convert', methods=["POST"])
def convert():
    data = request.get_json()
    mol_block = data['data']
    args = None
    if "selections" in data:
        options = data['selections']
        angle = str(data['angle'])
        indentation = str(data['indentation'])
        h2 = data['h2']
        args = combine_args(options, angle, indentation, h2)
    # Validate a drawn structure by RDKit
    mol = Chem.MolFromMolBlock(mol_block)
    if mol is None:
        return jsonify(
            {
                "error": "Chemfig cannot be generated. Please check structure.",
            }
        )
    chemfig, pdflink, error = convert_mol_format(mol_block, args=args)
    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink
        }
    )


@m2cf.route('/search', methods=["POST", "GET"])
def search():
    # Support both POST (JSON body) and GET (query param)
    if request.method == "GET":
        search_term = request.args.get('searchTerm', '').strip()
    else:
        data = request.get_json()
        search_term = data['searchTerm'].strip()

    # Use PubChem REST API to convert nomenclature to SMILES
    error, smiles = get_smiles_from_pubchem_api(search_term)
    
    if smiles is None:
        return jsonify(
            {
                "error": error,
                "smiles": None,
                "chemfig": None,
                "pdflink": None
            }
        )

    # Successfully got SMILES - now convert to chemfig
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify(
                {
                    "error": "Invalid SMILES from PubChem: {}".format(smiles),
                    "smiles": smiles,
                    "chemfig": None,
                    "pdflink": None
                }
            )

        Chem.Kekulize(mol)
        smiles2 = Chem.MolToSmiles(mol, kekuleSmiles=True)

        mol2 = Chem.MolFromSmiles(smiles2)
        AllChem.Compute2DCoords(mol)
        mol_block = Chem.MolToMolBlock(mol2)

        chemfig, pdflink, error = smiles_mol_to_chemfig("-w",
                                                        '-i direct {}'
                                                        .format(smiles))
        return jsonify(
            {
                "error": error,
                "smiles": smiles,
                "chemfig": chemfig,
                "pdflink": pdflink,
                "molblock": mol_block
            }
        )
    except Exception as e:
        return jsonify(
            {
                "error": "Failed to generate chemfig: {}".format(str(e)),
                "smiles": smiles,
                "chemfig": None,
                "pdflink": None
            }
        )


@m2cf.route('/submit', methods=["POST"])
def submit():
    data = request.get_json()
    text_area_data = data['textAreaData'].strip()

    # for molfiles
    if "END" in text_area_data:
        try:
            chemfig, pdflink, error = convert_mol_format(text_area_data)
            chem_data = text_area_data
            chem_format = "mol"
        except Exception as e:
            return jsonify(
                {
                    "error": "Sorry, Chemfig cannot be generated. Check your\
                    MOL format",
                    "chemfig": None,
                    "pdflink": None
                }
            )

    # for smiles
    else:
        try:
            mol = Chem.MolFromSmiles(text_area_data)
            Chem.Kekulize(mol)
        # if a user inputes other than MOL, smiles or chemfig
        except Exception as e:
            print(e)

            return jsonify(
                {
                    "error": "Sorry, Chemfig cannot be generated",
                    "chemfig": None,
                    "pdflink": None
                }
            )

        smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
        chemfig, pdflink, error = smiles_mol_to_chemfig("-w",
                                                        '-i direct {}'
                                                        .format(smiles))
        chem_data = smiles
        chem_format = 'smiles'

    # Generate SVG from chemfig code
    svglink = None
    if chemfig:
        svg_content, svg_error = chemfig_to_svg(chemfig)
        if svg_content and not svg_error:
            svglink = svg_content
        else:
            print("SVG generation error:", svg_error)

    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink,
            "svglink": svglink,
            "chem_format": chem_format,
            "chem_data": chem_data
        }
    )

# @m2cf.route('/submit', methods=["POST"])
# def submit():
    # data = request.get_json()
    # chem_data = data["data"]
    # data_type = data["dataType"]
    # if data_type == "mol":
        # chemfig, pdflink, error = convert_mol_format(chem_data)
        # return jsonify(
            # {
                # "error": error,
                # "chemfig": chemfig,
                # "pdflink": pdflink
            # }
        # )


@m2cf.route('/update_chemfig', methods=["POST"])
def modify_chemfig():

    error = None

    data = request.get_json()
    chemfig = data["textAreaData"]

    try:
        pdf_link = update_chemfig(chemfig.strip())
        return jsonify(
            {
                "pdflink": pdf_link,
                "error": error
            }
        )
    except Exception as e:
        error = e
        return jsonify(
            {
                "pdflink": "",
                "error": error
            }
        )


@m2cf.route('/apply', methods=["POST"])
def apply():
    error = None
    data = request.get_json()
    chem_data = data['chem_data']
    chem_format = data['chem_format']
    options = data['selections']
    angle = str(data['angle'])
    indentation = str(data['indentation'])
    h2 = data['h2']

    args = combine_args(options, angle, indentation, h2)

    if chem_format == "mol":
        chemfig, pdflink, error = convert_mol_format(chem_data, args=args)
    else:
        chemfig, pdflink, error = smiles_mol_to_chemfig("-w "
                                                        + args
                                                        + " -i direct {}"
                                                        .format(chem_data))

    # Generate SVG from chemfig code
    svglink = None
    if chemfig:
        svg_content, svg_error = chemfig_to_svg(chemfig)
        if svg_content and not svg_error:
            svglink = svg_content

    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink,
            "svglink": svglink
        }
    )


@m2cf.route('/reset', methods=["POST"])
def reset():
    data = request.get_json()
    chem_data = data['chem_data']
    chem_format = data['chem_format']

    if chem_format == "mol":
        chemfig, pdflink, error = convert_mol_format(chem_data, args=None)
    else:
        chemfig, pdflink, error = smiles_mol_to_chemfig("-w "
                                                        + " -i direct {}"
                                                        .format(chem_data))

    # Generate SVG from chemfig code
    svglink = None
    if chemfig:
        svg_content, svg_error = chemfig_to_svg(chemfig)
        if svg_content and not svg_error:
            svglink = svg_content

    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink,
            "svglink": svglink
        }
    )


@m2cf.route('/layers', methods=["POST"])
def layers():
    """
    Generate separate SVG layers for layered composition.
    Returns base molecule + individual layers for each option.
    This enables efficient client-side composition instead of 2^8 combinations.
    """
    error = None
    data = request.get_json()
    chem_data = data['chem_data']
    chem_format = data['chem_format']
    
    # Options to generate layers for (excluding transformations like angle/flip)
    layer_options = [
        {'name': 'aromatic', 'flag': '-o', 'description': 'Aromatic circles'},
        {'name': 'carbon', 'flag': '-c', 'description': 'Show carbon'},
        {'name': 'methyl', 'flag': '-m', 'description': 'Show methyl'},
        {'name': 'numbers', 'flag': '-n', 'description': 'Atom numbers'},
        {'name': 'fancy', 'flag': '-f', 'description': 'Fancy bonds'},
        {'name': 'compact', 'flag': '-z', 'description': 'Compact view'}
    ]
    
    layers_result = {}
    base_svg_for_diff = None  # Store base SVG for diff extraction
    
    # Generate BASE layer (no options)
    try:
        if chem_format == "mol":
            base_chemfig, _, base_error = convert_mol_format(chem_data, args=None)
        else:
            base_chemfig, _, base_error = smiles_mol_to_chemfig("-w -i direct {}".format(chem_data))
        
        if base_chemfig and not base_error:
            base_svg, svg_error = chemfig_to_svg(base_chemfig)
            if base_svg and not svg_error:
                layers_result['base'] = base_svg
                layers_result['base_chemfig'] = base_chemfig
                base_svg_for_diff = base_svg  # Store for diffing
            else:
                error = svg_error
        else:
            error = base_error
    except Exception as e:
        error = str(e)
    
    # Generate individual OPTION layers (FULL molecule with each option)
    # Strategy: Return complete molecule renderings, use CSS to switch between them
    for option in layer_options:
        try:
            args = combine_args([option['flag']], '0', '4', 'keep')
            
            if chem_format == "mol":
                layer_chemfig, _, layer_error = convert_mol_format(chem_data, args=args)
            else:
                layer_chemfig, _, layer_error = smiles_mol_to_chemfig("-w " + args + " -i direct {}".format(chem_data))
            
            if layer_chemfig and not layer_error:
                layer_svg_full, svg_error = chemfig_to_svg(layer_chemfig)
                if layer_svg_full and not svg_error:
                    # Store FULL SVG (not diff) - simpler and more reliable
                    layers_result[option['name']] = layer_svg_full
                    layers_result[option['name'] + '_chemfig'] = layer_chemfig
                    print("Generated FULL layer for {}: {} chars".format(option['name'], len(layer_svg_full)))
        except Exception as e:
            print("Error generating {} layer: {}".format(option['name'], e))
    
    return jsonify({
        "error": error,
        "layers": layers_result,
        "chem_data": chem_data,
        "chem_format": chem_format
    })


@m2cf.route('/svg', methods=["POST"])
def svg():
    """
    Convert SMILES or MOL to SVG with transparent background.
    Direct SVG generation using RDKit.
    """
    if not HAS_SVG_SUPPORT:
        return jsonify(
            {
                "error": "SVG generation not available (RDKit Draw module not found)",
                "svg": None
            }
        ), 500
    
    data = request.get_json()
    text_area_data = data['textAreaData'].strip()
    
    try:
        # Parse SMILES or MOL
        if "END" in text_area_data:
            # MOL block format
            mol = Chem.MolFromMolBlock(text_area_data)
        else:
            # SMILES format
            mol = Chem.MolFromSmiles(text_area_data)
        
        if mol is None:
            return jsonify(
                {
                    "error": "Invalid chemical structure format",
                    "svg": None
                }
            ), 400
        
        # Generate 2D coordinates if needed
        AllChem.Compute2DCoords(mol)
        
        # Draw to SVG with transparent background
        svg_drawer = Draw.MolDraw2DSVG(400, 300)
        svg_drawer.DrawMolecule(mol)
        svg_drawer.FinishDrawing()
        
        svg_text = svg_drawer.GetDrawingText()
        
        # Remove or modify background to be transparent
        svg_text = svg_text.replace('fill="#FFFFFF"', 'fill="none"')
        svg_text = svg_text.replace('fill="#ffffff"', 'fill="none"')
        svg_text = svg_text.replace('fill="white"', 'fill="none"')
        # Also remove the rect element that creates the background
        import re
        svg_text = re.sub(r'<rect[^>]*style=\'opacity:1\.0;fill:#FFFFFF[^>]*>\s*</rect>', '', svg_text)
        svg_text = re.sub(r'<rect[^>]*style=\'opacity:1\.0;fill:#ffffff[^>]*>\s*</rect>', '', svg_text)
        
        return jsonify(
            {
                "error": None,
                "svg": svg_text,
                "format": "svg",
                "width": 400,
                "height": 300
            }
        ), 200
        
    except Exception as e:
        return jsonify(
            {
                "error": "SVG generation failed: {}".format(str(e)),
                "svg": None
            }
        ), 500


@m2cf.route('/nomenclature-to-svg', methods=["POST"])
def nomenclature_to_svg():
    """
    Unified endpoint: nomenclature/name -> SMILES -> SVG
    Combines nomenclature lookup with SVG generation.
    """
    if not HAS_SVG_SUPPORT:
        return jsonify(
            {
                "error": "SVG generation not available",
                "svg": None,
                "smiles": None,
                "nomenclature": None
            }
        ), 500
    
    data = request.get_json()
    nomenclature = data.get('nomenclature', '').strip()
    
    if not nomenclature:
        return jsonify(
            {
                "error": "nomenclature parameter required",
                "svg": None,
                "smiles": None
            }
        ), 400
    
    try:
        # Step 1: Convert nomenclature to SMILES using PubChem
        error, smiles = get_smiles_from_pubchem_api(nomenclature)
        
        if smiles is None:
            return jsonify(
                {
                    "error": error or "Nomenclature not found in database",
                    "svg": None,
                    "smiles": None,
                    "nomenclature": nomenclature
                }
            ), 404
        
        # Step 2: Convert SMILES to SVG
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify(
                {
                    "error": "Invalid SMILES from nomenclature lookup: {}".format(smiles),
                    "svg": None,
                    "smiles": smiles,
                    "nomenclature": nomenclature
                }
            ), 400
        
        AllChem.Compute2DCoords(mol)
        
        # Draw to SVG
        svg_drawer = Draw.MolDraw2DSVG(400, 300)
        svg_drawer.DrawMolecule(mol)
        svg_drawer.FinishDrawing()
        
        svg_text = svg_drawer.GetDrawingText()
        
        # Make background transparent
        svg_text = svg_text.replace('fill="#FFFFFF"', 'fill="none"')
        svg_text = svg_text.replace('fill="#ffffff"', 'fill="none"')
        svg_text = svg_text.replace('fill="white"', 'fill="none"')
        # Also remove the rect element that creates the background
        import re
        svg_text = re.sub(r'<rect[^>]*style=\'opacity:1\.0;fill:#FFFFFF[^>]*>\s*</rect>', '', svg_text)
        svg_text = re.sub(r'<rect[^>]*style=\'opacity:1\.0;fill:#ffffff[^>]*>\s*</rect>', '', svg_text)
        
        return jsonify(
            {
                "error": None,
                "svg": svg_text,
                "smiles": smiles,
                "nomenclature": nomenclature,
                "format": "svg",
                "width": 400,
                "height": 300
            }
        ), 200
        
    except Exception as e:
        return jsonify(
            {
                "error": "Processing failed: {}".format(str(e)),
                "svg": None,
                "smiles": None,
                "nomenclature": nomenclature
            }
        ), 500


# REACTION ###

@m2cf.route("/reaction/convert", methods=["POST"])
def convert_reaction():
    data = request.get_json()
    json_doc = json.loads(data['docJSON'])
    mol_files = data["mol_files"]

    reaction_chemfig, txt_chemfig = generate_latex(json_doc, mol_files, None)
    pdf_link = chemfig_to_pdf(reaction_chemfig)

    return jsonify(
        {
            "OK": "ALL IS GOOD",
            "chemfig": txt_chemfig,
            "pdflink": pdf_link

        }
    )


@m2cf.route("/reaction/apply", methods=["POST"])
def apply_reaction():
    data = request.get_json()
    json_doc = json.loads(data['docJSON'])
    mol_files = data["mol_files"]
    options = data['selections']
    angle = str(data['angle'])
    indentation = str(data['indentation'])
    h2 = data['h2']

    args = combine_args(options, angle, indentation, h2)

    reaction_chemfig, txt_chemfig = generate_latex(json_doc, mol_files, args)
    pdf_link = chemfig_to_pdf(reaction_chemfig)

    return jsonify(
        {
            "OK": "APPLY",
            "chemfig": txt_chemfig,
            "pdflink": pdf_link
        }
    )


@m2cf.route('/reaction/reset', methods=["POST"])
def reaction_reset():
    data = request.get_json()
    data = request.get_json()
    json_doc = json.loads(data['docJSON'])
    mol_files = data["mol_files"]

    reaction_chemfig, txt_chemfig = generate_latex(json_doc, mol_files, None)
    pdf_link = chemfig_to_pdf(reaction_chemfig)

    return jsonify(
        {
            "OK": "APPLY",
            "chemfig": txt_chemfig,
            "pdflink": pdf_link
        }
    )


@m2cf.route('/reaction/update_chemfig', methods=["POST"])
def reaction_chemfig():

    error = None

    data = request.get_json()
    chemfig = data["textAreaData"]
    chemfig = update_reaction_chemfig(chemfig.strip())

    try:
        pdf_link = update_chemfig(chemfig)
        return jsonify(
            {
                "pdflink": pdf_link,
                "error": error
            }
        )
    except Exception as e:
        error = e
        return jsonify(
            {
                "pdflink": "",
                "error": error
            }
        )
