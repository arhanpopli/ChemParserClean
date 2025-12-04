"""
Test mol2chemfig LaTeX compilation with a minimal example.
"""
import subprocess
import tempfile
import os

# Simple test LaTeX document
latex_doc = r"""\documentclass[border=2pt]{standalone}
\usepackage[utf8]{inputenc}
\usepackage{chemfig}
\begin{document}
\chemfig{C-C-O}
\end{document}
"""

print("Testing LaTeX compilation...")

with tempfile.TemporaryDirectory() as tmpdir:
    tex_file = os.path.join(tmpdir, "test.tex")
    with open(tex_file, 'w', encoding='utf-8') as f:
        f.write(latex_doc)
    
    print(f"Compiling {tex_file}...")
    proc = subprocess.run(
        ['pdflatex', '-interaction=nonstopmode', 'test.tex'],
        cwd=tmpdir,
        capture_output=True,
        text=True,
        timeout=10
    )
    
    if proc.returncode == 0:
        print("✅ pdflatex succeeded!")
        
        # Try dvisvgm
        print("Converting to SVG...")
        proc_svg = subprocess.run(
            ['dvisvgm', '--pdf', 'test.pdf', '-n', '-o', 'test.svg'],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if proc_svg.returncode == 0:
            print("✅ dvisvgm succeeded!")
            svg_file = os.path.join(tmpdir, "test.svg")
            if os.path.exists(svg_file):
                with open(svg_file, 'r') as f:
                    svg = f.read()
                print(f"✅ SVG generated ({len(svg)} bytes)")
                
                # Save it
                with open("test_latex_output.svg", "w") as f:
                    f.write(svg)
                print("Saved to test_latex_output.svg")
        else:
            print("❌ dvisvgm failed:")
            print(proc_svg.stderr)
    else:
        print("❌ pdflatex failed:")
        print(proc.stdout[-1000:])
