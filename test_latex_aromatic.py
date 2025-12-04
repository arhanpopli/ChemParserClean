"""
Test mol2chemfig aromatic circles specifically.
"""
import subprocess
import tempfile
import os
import shutil

# LaTeX with aromatic circle macro usage
# Based on what mol2chemfig generates for benzene
latex_doc = r"""\documentclass[border=2pt]{standalone}
\usepackage[utf8]{inputenc}
\usepackage{mol2chemfig}
\usepackage{chemfig}
\begin{document}
% Benzene with aromatic circle
\chemfig{
           -[:30,,2]
     -[:90,,,2]
    -[:150]
    -[:210]
    -[:270]
    -[:330]
}
\mcfcringle{6}
\end{document}
"""

print("Testing aromatic circles compilation...")

with tempfile.TemporaryDirectory() as tmpdir:
    # Copy mol2chemfig.sty
    sty_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig.sty')
    shutil.copy(sty_path, tmpdir)

    tex_file = os.path.join(tmpdir, "test_aromatic.tex")
    with open(tex_file, 'w', encoding='utf-8') as f:
        f.write(latex_doc)
    
    print(f"Compiling {tex_file}...")
    proc = subprocess.run(
        ['pdflatex', '-interaction=nonstopmode', 'test_aromatic.tex'],
        cwd=tmpdir,
        capture_output=True,
        text=True,
        timeout=10
    )
    
    if proc.returncode == 0:
        print("✅ pdflatex succeeded!")
    else:
        print("❌ pdflatex failed:")
        print(proc.stdout[-2000:])
