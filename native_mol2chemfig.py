"""
Native mol2chemfig renderer - using subprocess to call mol2chemfig CLI wrapper.
"""
import subprocess
import tempfile
import os
import shutil
import traceback

def run_mol2chemfig(smiles, options=None):
    """
    Convert SMILES to SVG using mol2chemfig command-line tool.
    
    Args:
        smiles (str): SMILES string
        options (dict): Options for mol2chemfig
        
    Returns:
        str: SVG content or None if failed
    """
    if options is None:
        options = {}
        
    print(f"[Native] Rendering SMILES: {smiles}")
    print(f"[Native] Options: {options}")
    
    # Build command line - use our CLI wrapper
    wrapper_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig_cli.py')
    cmd = ['python', wrapper_path]
    
    if options.get('m2cfShowCarbons'): cmd.append('-c')
    if options.get('m2cfAromaticCircles'): cmd.append('--aromatic-circles')
    if options.get('m2cfShowMethyls'): cmd.append('--show-methyls')
    if options.get('m2cfAtomNumbers'): cmd.append('-n')
    if options.get('m2cfFancyBonds'): cmd.append('--fancy-bonds')
    if options.get('m2cfFlipHorizontal'): cmd.append('--flip-horizontal')
    if options.get('m2cfFlipVertical'): cmd.append('--flip-vertical')
    if options.get('m2cfRotate', 0) != 0: 
        cmd.extend(['--rotate', str(options['m2cfRotate'])])
    
    hydro_mode = options.get('m2cfHydrogensMode', 'keep')
    if hydro_mode == 'add': cmd.append('-h')
    elif hydro_mode == 'delete': cmd.append('-d')
    
    print(f"[Native] Command: {' '.join(cmd)}")
    print(f"[Native] SMILES (stdin): {smiles}")
    
    try:
        # Run mol2chemfig and get LaTeX output
        proc = subprocess.run(
            cmd,
            input=smiles,
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if proc.returncode != 0:
            print(f"[Native] mol2chemfig failed:")
            print(proc.stderr)
            return None
            
        chemfig_code = proc.stdout.strip()
        print(f"[Native] LaTeX Code:\n{chemfig_code}")
        
        # Wrap in standalone document
        latex_doc = r"""\documentclass[border=2pt]{standalone}
\usepackage[utf8]{inputenc}
\usepackage{mol2chemfig}
\usepackage{chemfig}
\setatomsep{2em}
\renewcommand*\printatom[1]{\ensuremath{\mathsf{#1}}}
\begin{document}
%s
\end{document}
""" % chemfig_code
        
        # Compile to PDF then SVG
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy mol2chemfig.sty to tmpdir
            sty_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig.sty')
            if os.path.exists(sty_path):
                shutil.copy(sty_path, tmpdir)
            else:
                print(f"[Native] Warning: mol2chemfig.sty not found at {sty_path}")

            tex_file = os.path.join(tmpdir, "mol.tex")
            with open(tex_file, 'w', encoding='utf-8') as f:
                f.write(latex_doc)
            
            # pdflatex
            print(f"[Native] Compiling LaTeX...")
            proc_latex = subprocess.run(
                ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'mol.tex'],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if proc_latex.returncode != 0:
                print("[Native] pdflatex failed:")
                print("--- STDOUT ---")
                print(proc_latex.stdout)
                
                # Save failed tex file
                failed_tex = os.path.join(os.path.dirname(__file__), "failed_mol.tex")
                shutil.copy(tex_file, failed_tex)
                print(f"[Native] Saved failed LaTeX to {failed_tex}")
                
                return None
            
            # dvisvgm
            print(f"[Native] Converting to SVG...")
            proc_svg = subprocess.run(
                ['dvisvgm', '--pdf', 'mol.pdf', '-n', '-o', 'mol.svg'],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if proc_svg.returncode != 0:
                print("[Native] dvisvgm failed")
                print(proc_svg.stderr)
                return None
            
            svg_file = os.path.join(tmpdir, "mol.svg")
            if os.path.exists(svg_file):
                with open(svg_file, 'r', encoding='utf-8') as f:
                    svg_content = f.read()
                print(f"[Native] ✅ Success! SVG size: {len(svg_content)} bytes")
                return svg_content
                
        return None
        
    except Exception as e:
        print(f"[Native] Error: {e}")
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Test
    svg = run_mol2chemfig("c1ccccc1", {'m2cfAromaticCircles': True, 'm2cfShowCarbons': True})
    if svg:
        print("\n✅ Test passed!")
        with open("test_output.svg", "w") as f:
            f.write(svg)
        print("Saved to test_output.svg")
    else:
        print("\n❌ Test failed")
