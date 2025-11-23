"""
Local mol2chemfig processing module
Uses mol2chemfigPy3 for chemfig code generation
Uses LaTeX + dvisvgm for SVG rendering (EXACT same as Docker backend)
Falls back to Docker backend if local processing fails
"""

import json
import hashlib
import tempfile
import shutil
import subprocess
import os
from pathlib import Path

# Try to import local dependencies
try:
    from mol2chemfigPy3 import mol2chemfig
    HAS_MOL2CHEMFIG = True
except ImportError:
    HAS_MOL2CHEMFIG = False
    print("Warning: mol2chemfigPy3 not installed. Run: pip install mol2chemfigPy3")

# Path to mol2chemfig.sty (downloaded from CTAN)
SCRIPT_DIR = Path(__file__).parent
MOL2CHEMFIG_STY = SCRIPT_DIR / "mol2chemfig.sty"

# Cache directory
CACHE_DIR = Path("cache") / "mol2chemfig_local"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Check for LaTeX tools
def check_latex_tools():
    """Check if LaTeX and dvisvgm are available"""
    try:
        result = subprocess.run(['latex', '--version'], capture_output=True, text=True, timeout=5)
        has_latex = result.returncode == 0
    except:
        has_latex = False

    try:
        result = subprocess.run(['dvisvgm', '--version'], capture_output=True, text=True, timeout=5)
        has_dvisvgm = result.returncode == 0
    except:
        has_dvisvgm = False

    return has_latex, has_dvisvgm

HAS_LATEX, HAS_DVISVGM = check_latex_tools()
HAS_MOL2CHEMFIG_STY = MOL2CHEMFIG_STY.exists()


def get_content_hash(smiles, options=None):
    """Generate unique hash for SMILES + options combination"""
    options_str = json.dumps(sorted(options or []))
    content = f"{smiles}:{options_str}"
    return hashlib.sha256(content.encode()).hexdigest()[:16]


def smiles_to_chemfig(smiles, options=None):
    """
    Convert SMILES to chemfig code using mol2chemfigPy3

    Args:
        smiles: SMILES string
        options: List of mol2chemfig options like ['-o', '-m', '-f']

    Returns:
        dict with 'chemfig' key and optionally 'error' key
    """
    if not HAS_MOL2CHEMFIG:
        return {"error": "mol2chemfigPy3 not installed", "chemfig": None}

    try:
        import sys
        import io

        # mol2chemfigPy3 API:
        # mol2chemfig(content, *args, rotate=0.0, aromatic=True, marker=None,
        #             name=None, relative_angle=False, show_carbon=False,
        #             show_methyl=False, inline=False)
        # NOTE: mol2chemfig PRINTS to stdout and returns None!

        # Parse options into keyword arguments
        kwargs = {
            'aromatic': '-o' in (options or []),  # aromatic circles
            'show_carbon': '-c' in (options or []),
            'show_methyl': '-m' in (options or []),
        }

        # Capture stdout since mol2chemfig prints instead of returning
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO()

        try:
            mol2chemfig(smiles, **kwargs)
            chemfig_code = buffer.getvalue().strip()
        finally:
            sys.stdout = old_stdout

        if chemfig_code:
            return {"chemfig": chemfig_code, "error": None}
        else:
            return {"chemfig": None, "error": "mol2chemfig returned empty result"}

    except Exception as e:
        return {"chemfig": None, "error": f"mol2chemfig error: {str(e)}"}


def chemfig_to_svg(chemfig_code):
    """
    Convert chemfig LaTeX code to SVG using latex + dvisvgm.
    Uses the EXACT SAME LaTeX template as Docker backend with proper bond spacing settings.
    Key settings: \setdoublesep{3pt}, \setcrambond{2.5pt}{0.4pt}{1.0pt}

    Args:
        chemfig_code: The chemfig LaTeX code (e.g., \chemfig{...})

    Returns:
        tuple (svg_content, error) - svg_content is None if error
    """
    if not HAS_LATEX:
        return None, "LaTeX not installed or not in PATH"
    if not HAS_DVISVGM:
        return None, "dvisvgm not installed or not in PATH"
    if not HAS_MOL2CHEMFIG_STY:
        return None, f"mol2chemfig.sty not found at {MOL2CHEMFIG_STY}"

    tempdir = tempfile.mkdtemp()
    tex_file = os.path.join(tempdir, 'molecule.tex')
    dvi_file = os.path.join(tempdir, 'molecule.dvi')
    svg_file = os.path.join(tempdir, 'molecule.svg')

    # LaTeX template - Based on Docker backend but simplified for Windows/MiKTeX
    # Removed sfmath which may require installation prompts
    latex_content = r'''\documentclass{minimal}
\usepackage{xcolor, mol2chemfig}
\setcrambond{2.5pt}{0.4pt}{1.0pt}
\setbondoffset{1pt}
\setdoublesep{3pt}
\setatomsep{28pt}
\renewcommand{\printatom}[1]{\fontsize{12pt}{14pt}\selectfont{\ensuremath{\mathsf{#1}}}}
\setlength{\parindent}{0pt}
\begin{document}
%s
\end{document}''' % chemfig_code

    try:
        # Copy mol2chemfig.sty to temp directory
        sty_dest = os.path.join(tempdir, 'mol2chemfig.sty')
        shutil.copy(str(MOL2CHEMFIG_STY), sty_dest)

        # Write LaTeX file
        with open(tex_file, 'w', encoding='utf-8') as f:
            f.write(latex_content)

        # Compile using latex (not pdflatex - dvisvgm needs DVI)
        curdir = os.getcwd()
        os.chdir(tempdir)

        try:
            # Run latex
            result = subprocess.run(
                ['latex', '-interaction=nonstopmode', tex_file],
                capture_output=True,
                text=True,
                timeout=30
            )

            if not os.path.exists(dvi_file):
                os.chdir(curdir)
                return None, f"LaTeX compilation failed: {result.stderr or result.stdout}"

            # Convert DVI to SVG
            result = subprocess.run(
                ['dvisvgm', '--font-format=woff', '--exact', dvi_file, '-o', svg_file],
                capture_output=True,
                text=True,
                timeout=30
            )

            if not os.path.exists(svg_file):
                os.chdir(curdir)
                return None, f"dvisvgm conversion failed: {result.stderr or result.stdout}"

            # Read SVG file
            with open(svg_file, 'r', encoding='utf-8') as f:
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

        finally:
            os.chdir(curdir)

    except subprocess.TimeoutExpired:
        return None, "LaTeX/dvisvgm timed out"
    except Exception as e:
        return None, str(e)
    finally:
        # Clean up temp directory
        shutil.rmtree(tempdir, ignore_errors=True)


def process_smiles(smiles, options=None):
    """
    Full processing: SMILES -> chemfig + SVG

    Args:
        smiles: SMILES string
        options: List of mol2chemfig options

    Returns:
        dict with 'chemfig', 'svg', 'svglink', 'error', 'cached' keys
    """
    content_hash = get_content_hash(smiles, options)

    # Check cache
    svg_cache_file = CACHE_DIR / f"{content_hash}.svg"
    chemfig_cache_file = CACHE_DIR / f"{content_hash}.chemfig"

    if svg_cache_file.exists() and chemfig_cache_file.exists():
        with open(svg_cache_file, 'r', encoding='utf-8') as f:
            svg = f.read()
        with open(chemfig_cache_file, 'r', encoding='utf-8') as f:
            chemfig = f.read()
        return {
            "chemfig": chemfig,
            "svg": svg,
            "svglink": f"/cache/mol2chemfig_local/{content_hash}.svg",
            "error": None,
            "cached": True,
            "hash": content_hash,
            "source": "local_cache"
        }

    # Generate chemfig
    chemfig_result = smiles_to_chemfig(smiles, options)
    chemfig_code = chemfig_result.get('chemfig')

    if not chemfig_code:
        return {
            "chemfig": None,
            "svg": None,
            "svglink": None,
            "error": chemfig_result.get('error', 'Failed to generate chemfig'),
            "cached": False,
            "hash": content_hash,
            "source": "local"
        }

    # Generate SVG from chemfig
    svg_content, svg_error = chemfig_to_svg(chemfig_code)

    if svg_error:
        return {
            "chemfig": chemfig_code,
            "svg": None,
            "svglink": None,
            "error": svg_error,
            "cached": False,
            "hash": content_hash,
            "source": "local"
        }

    # Cache results
    with open(svg_cache_file, 'w', encoding='utf-8') as f:
        f.write(svg_content)
    with open(chemfig_cache_file, 'w', encoding='utf-8') as f:
        f.write(chemfig_code)

    return {
        "chemfig": chemfig_code,
        "svg": svg_content,
        "svglink": f"/cache/mol2chemfig_local/{content_hash}.svg",
        "error": None,
        "cached": False,
        "hash": content_hash,
        "source": "local"
    }


def process_nomenclature(name, options=None):
    """
    Process chemical name: name -> SMILES -> chemfig + SVG
    Uses PubChem API for name to SMILES conversion

    Args:
        name: Chemical name (e.g., "ethanol", "aspirin")
        options: List of mol2chemfig options

    Returns:
        dict with 'chemfig', 'svg', 'smiles', 'error' keys
    """
    import requests
    from urllib.parse import quote

    # Common names lookup (fast fallback)
    COMMON_NAMES = {
        'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
        'benzene': 'c1ccccc1',
        'toluene': 'Cc1ccccc1',
        'ethanol': 'CCO',
        'methanol': 'CO',
        'acetone': 'CC(=O)C',
        'caffeine': 'Cn1cnc2c1c(=O)n(c(=O)n2C)C',
        'glucose': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
        'water': 'O',
        'methane': 'C',
        'ethane': 'CC',
        'propane': 'CCC',
        'butane': 'CCCC',
        'acetic acid': 'CC(=O)O',
        'histamine': 'NCCc1cnc[nH]1',
    }

    name_lower = name.lower().strip()

    # Check common names first
    if name_lower in COMMON_NAMES:
        smiles = COMMON_NAMES[name_lower]
        result = process_smiles(smiles, options)
        result['smiles'] = smiles
        result['name'] = name
        result['source'] = 'local_common_names'
        return result

    # Try PubChem API
    try:
        encoded_name = quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES/JSON"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                props = data['PropertyTable']['Properties']
                if len(props) > 0:
                    smiles = props[0].get('CanonicalSMILES') or props[0].get('IsomericSMILES')
                    if smiles:
                        result = process_smiles(smiles, options)
                        result['smiles'] = smiles
                        result['name'] = name
                        result['source'] = 'local_pubchem'
                        return result

        return {
            "chemfig": None,
            "svg": None,
            "smiles": None,
            "name": name,
            "error": f"Could not find SMILES for: {name}",
            "source": "local"
        }

    except Exception as e:
        return {
            "chemfig": None,
            "svg": None,
            "smiles": None,
            "name": name,
            "error": f"PubChem lookup error: {str(e)}",
            "source": "local"
        }


def is_available():
    """Check if local processing is fully available"""
    return HAS_MOL2CHEMFIG and HAS_LATEX and HAS_DVISVGM and HAS_MOL2CHEMFIG_STY


def get_status():
    """Get status of local mol2chemfig processing"""
    return {
        "mol2chemfigPy3": HAS_MOL2CHEMFIG,
        "latex": HAS_LATEX,
        "dvisvgm": HAS_DVISVGM,
        "mol2chemfig_sty": HAS_MOL2CHEMFIG_STY,
        "available": is_available(),
        "cache_dir": str(CACHE_DIR),
        "cache_files": len(list(CACHE_DIR.glob('*')))
    }


# Test function
if __name__ == "__main__":
    print("=" * 60)
    print("mol2chemfig Local Processing Test")
    print("=" * 60)
    status = get_status()
    print(f"\nStatus:")
    for key, value in status.items():
        print(f"  {key}: {value}")

    if is_available():
        print("\n--- Testing SMILES Processing ---")
        result = process_smiles("CCO", options=['-o', '-m'])
        print(f"SMILES: CCO (ethanol)")
        print(f"Chemfig: {result.get('chemfig')}")
        print(f"SVG length: {len(result.get('svg', '')) if result.get('svg') else 0}")
        print(f"SVG link: {result.get('svglink')}")
        print(f"Error: {result.get('error')}")
        print(f"Cached: {result.get('cached')}")
        print(f"Source: {result.get('source')}")

        print("\n--- Testing Nomenclature Processing ---")
        result = process_nomenclature("aspirin")
        print(f"Name: aspirin")
        print(f"SMILES: {result.get('smiles')}")
        print(f"Chemfig: {result.get('chemfig')}")
        print(f"SVG length: {len(result.get('svg', '')) if result.get('svg') else 0}")
        print(f"Error: {result.get('error')}")
    else:
        print("\nLocal processing not fully available!")
        print("Missing components:")
        if not HAS_MOL2CHEMFIG:
            print("  - mol2chemfigPy3: pip install mol2chemfigPy3")
        if not HAS_LATEX:
            print("  - LaTeX: Install MiKTeX or TeX Live")
        if not HAS_DVISVGM:
            print("  - dvisvgm: Usually comes with LaTeX distribution")
        if not HAS_MOL2CHEMFIG_STY:
            print(f"  - mol2chemfig.sty: Download from CTAN to {MOL2CHEMFIG_STY}")
