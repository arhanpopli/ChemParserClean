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


def get_smiles_from_pubchem_api(compound_name):
    """
    Fetch SMILES from PubChem REST API using compound name
    Returns (error, smiles)
    """
    try:
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
        
        return "Compound '{}' not found in PubChem database. Try pasting SMILES directly.".format(compound_name), None
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


@m2cf.route('/search', methods=["POST"])
def search():
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

    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink,
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

    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink
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

    return jsonify(
        {
            "error": error,
            "chemfig": chemfig,
            "pdflink": pdflink
        }
    )


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
