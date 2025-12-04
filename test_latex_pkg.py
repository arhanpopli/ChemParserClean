"""
Test mol2chemfig package loading.
"""
import subprocess
import tempfile
import os
import shutil

# Minimal LaTeX using mol2chemfig package
latex_doc = r"""\documentclass[border=2pt]{standalone}
\usepackage[utf8]{inputenc}
\usepackage{mol2chemfig}
\begin{document}
\chemfig{C-C}
\end{document}
"""

print("Testing mol2chemfig package loading...")

with tempfile.TemporaryDirectory() as tmpdir:
    # Copy mol2chemfig.sty
    sty_path = os.path.join(os.path.dirname(__file__), 'mol2chemfig.sty')
    shutil.copy(sty_path, tmpdir)
    print(f"Copied mol2chemfig.sty to {tmpdir}")

    tex_file = os.path.join(tmpdir, "test_pkg.tex")
    with open(tex_file, 'w', encoding='utf-8') as f:
        f.write(latex_doc)
    
    print(f"Compiling {tex_file}...")
    # Remove -halt-on-error to see if it asks for something
    # But use a timeout so it doesn't hang forever
    try:
        proc = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', 'test_pkg.tex'],
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
            
            log_file = os.path.join(tmpdir, "test_pkg.log")
            if os.path.exists(log_file):
                print("\n--- LOG FILE ---")
                with open(log_file, 'r', errors='ignore') as f:
                    print(f.read()[-2000:])
                    
    except subprocess.TimeoutExpired:
        print("❌ Timeout! pdflatex hung.")
