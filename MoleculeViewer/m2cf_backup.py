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

# MoleculeViewer integration (OPSIN + SVG rendering)
MV_IMPORT_ERROR = None
MV_INTEGRATION_AVAILABLE = False
try:
    from MoleculeViewer.app import config as mv_config
    from MoleculeViewer.app.chemistry import (
        nomenclature_to_smiles as mv_nomenclature_to_smiles,
        smiles_to_svg as mv_smiles_to_svg,
        get_molecule_info as mv_get_molecule_info
    )
    MV_INTEGRATION_AVAILABLE = True
    mv_config.ENABLE_CHEMDOODLE = False  # Skip ChemDoodle DB per integration requirements
except Exception as exc:  # pragma: no cover - best effort optional dependency
    MV_IMPORT_ERROR = exc
    mv_config = None
    mv_nomenclature_to_smiles = None
    mv_smiles_to_svg = None
    mv_get_molecule_info = None


DEFAULT_SVG_WIDTH = 720
DEFAULT_SVG_HEIGHT = 540
DEFAULT_SVG_OPTIONS = {
    "fancy_bonds": True,
    "aromatic_circles": True,
    "hydrogens": "keep"
}

# Import SVG generation functions
try:
    from rdkit.Chem import Draw
    HAS_SVG_SUPPORT = True
except ImportError:
    HAS_SVG_SUPPORT = False

m2cf = Blueprint('m2cf', __name__, url_prefix="/m2cf")

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
    """Unified nomenclature search combining MoleculeViewer (OPSIN + PubChem) and mol2chemfig."""

    def _to_int(value, default):
        try:
            if value in (None, ''):
                return default
            return int(value)
        except (TypeError, ValueError):
            return default

    svg_options = {}
    width = DEFAULT_SVG_WIDTH
    height = DEFAULT_SVG_HEIGHT

    if request.method == "GET":
        search_term = request.args.get('searchTerm', '').strip()
        width = _to_int(request.args.get('width') or request.args.get('svgWidth'), DEFAULT_SVG_WIDTH)
        height = _to_int(request.args.get('height') or request.args.get('svgHeight'), DEFAULT_SVG_HEIGHT)
        raw_options = request.args.get('svgOptions') or request.args.get('svg_options')
        if raw_options:
            try:
                svg_options = json.loads(raw_options)
            except Exception:
                svg_options = {}
    else:
        data = request.get_json() or {}
        search_term = (data.get('searchTerm') or '').strip()
        width = _to_int(data.get('width'), DEFAULT_SVG_WIDTH)
        height = _to_int(data.get('height'), DEFAULT_SVG_HEIGHT)
        svg_options = data.get('svgOptions') or data.get('svg_options') or {}

    if not search_term:
        return jsonify(
            {
                "error": "Please provide a searchTerm",
                "smiles": None,
                "chemfig": None,
                "pdflink": None,
                "svg": None,
                "molblock": None,
                "source": None,
                "info": None,
                "warnings": ["Missing search term"]
            }
        )

    merged_svg_options = DEFAULT_SVG_OPTIONS.copy()
    if isinstance(svg_options, dict):
        merged_svg_options.update(svg_options)

    smiles = None
    source = None
    svg_markup = None
    mol_info = None
    warnings = []

    if MV_INTEGRATION_AVAILABLE and mv_nomenclature_to_smiles:
        try:
            mv_error, mv_smiles, mv_source = mv_nomenclature_to_smiles(search_term)
        except Exception as exc:  # pragma: no cover - defensive guard
            mv_error = f"MoleculeViewer lookup failed: {exc}"
            mv_smiles = None
            mv_source = None

        if mv_smiles:
            smiles = mv_smiles
            source = mv_source
        if mv_error and not mv_smiles:
            warnings.append(mv_error)

    if smiles is None:
        fallback_error, fallback_smiles = get_smiles_from_pubchem_api(search_term)
        if fallback_smiles:
            smiles = fallback_smiles
            source = source or "PubChem API"
        if fallback_error and fallback_smiles is None:
            warnings.append(fallback_error)

    if smiles is None:
        return jsonify(
            {
                "error": warnings[-1] if warnings else f"No SMILES found for '{search_term}'",
                "smiles": None,
                "chemfig": None,
                "pdflink": None,
                "svg": None,
                "molblock": None,
                "source": None,
                "info": None,
                "warnings": warnings
            }
        )

    # Skip MoleculeViewer SVG rendering - use mol2chemfig's native chemfig + PDF output
    # SVG generation can be added later if needed via mol2chemfig's own pipeline

    chemfig = None
    pdflink = None
    mol_block = None
    chemfig_error = None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            chemfig_error = f"Invalid SMILES returned for '{search_term}': {smiles}"
        else:
            Chem.Kekulize(mol)
            smiles_kekule = Chem.MolToSmiles(mol, kekuleSmiles=True)
            mol2 = Chem.MolFromSmiles(smiles_kekule)
            AllChem.Compute2DCoords(mol)
            mol_block = Chem.MolToMolBlock(mol2)
            chemfig, pdflink, chemfig_error = smiles_mol_to_chemfig(
                "-w",
                '-i direct {}'.format(smiles_kekule)
            )
            smiles = smiles_kekule or smiles
    except Exception as exc:
        chemfig_error = f"Failed to generate chemfig: {exc}"

    if chemfig_error and chemfig is None:
        warnings.append(chemfig_error)

    return jsonify(
        {
            "error": chemfig_error if chemfig_error and chemfig is None else None,
            "smiles": smiles,
            "chemfig": chemfig,
            "pdflink": pdflink,
            "svg": svg_markup,
            "molblock": mol_block,
            "source": source,
            "info": mol_info,
            "warnings": warnings,
            "svg_options": merged_svg_options
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
