'''
Generate SVG from RDKit molecules with transparent background.
Alternative to PDF/image generation for web display.
'''
import os
import shutil
import tempfile
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import StringIO  # Python 2 compatible

def smiles_to_svg(smiles, width=400, height=300):
    """
    Convert SMILES directly to SVG using RDKit.
    Returns SVG string with transparent background.
    
    Args:
        smiles: SMILES string
        width: SVG width in pixels
        height: SVG height in pixels
    
    Returns:
        (success, svg_string, error_msg)
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, None, "Invalid SMILES: {}".format(smiles)
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Draw to SVG with transparent background
        drawer = Draw.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText()
        
        # Remove background rectangle for transparency
        svg = svg.replace('style="fill:#FFFFFF;stroke:none" width=', 
                         'style="fill:none;stroke:none" width=')
        svg = svg.replace('style="fill:#ffffff;stroke:none" width=', 
                         'style="fill:none;stroke:none" width=')
        
        return True, svg, None
        
    except Exception as e:
        return False, None, str(e)


def chemfig_to_svg_latex(chemfig_code, width=400, height=300):
    """
    Convert ChemFig LaTeX code to SVG using pdflatex + pdf2svg.
    
    Args:
        chemfig_code: LaTeX chemfig code string
        width: Desired SVG width
        height: Desired SVG height
    
    Returns:
        (success, svg_string, error_msg)
    """
    try:
        tempdir = tempfile.mkdtemp()
        m2pkg_path = '/usr/src/app/src/mol2chemfig'
        pkg = '/mol2chemfig.sty'
        
        # Create symlink to chemfig style
        try:
            os.symlink(m2pkg_path + pkg, tempdir + pkg)
        except:
            pass  # Symlink might already exist
        
        # Create LaTeX document
        latex_template = r'''
\documentclass[border=0pt,tikz]{standalone}
\usepackage{tikz}
\usepackage{mol2chemfig}
\begin{document}
\begin{tikzpicture}
%s
\end{tikzpicture}
\end{document}
''' % chemfig_code
        
        latexfn = os.path.join(tempdir, 'molecule.tex')
        pdfname = os.path.join(tempdir, 'molecule.pdf')
        svgname = os.path.join(tempdir, 'molecule.svg')
        
        # Write LaTeX file
        with open(latexfn, 'w') as f:
            f.write(latex_template)
        
        curdir = os.getcwd()
        os.chdir(tempdir)
        
        # Compile LaTeX to PDF
        latexcmd = 'pdflatex -interaction=nonstopmode molecule.tex > /dev/null 2>&1'
        ret = os.system(latexcmd)
        
        if ret != 0:
            os.chdir(curdir)
            shutil.rmtree(tempdir)
            return False, None, "LaTeX compilation failed"
        
        # Convert PDF to SVG using pdf2svg
        pdf2svgcmd = 'pdf2svg {} {} -'.format(pdfname, svgname)
        ret = os.system(pdf2svgcmd)
        
        svg = None
        if ret == 0 or os.path.exists(svgname):
            try:
                with open(svgname, 'r') as f:
                    svg = f.read()
                    
                # Make background transparent
                svg = svg.replace('fill="white"', 'fill="none"')
                svg = svg.replace('fill="#ffffff"', 'fill="none"')
                svg = svg.replace('fill="#FFFFFF"', 'fill="none"')
            except:
                pass
        
        os.chdir(curdir)
        shutil.rmtree(tempdir)
        
        if svg:
            return True, svg, None
        else:
            return False, None, "SVG conversion failed"
            
    except Exception as e:
        return False, None, "SVG generation error: {}".format(str(e))
