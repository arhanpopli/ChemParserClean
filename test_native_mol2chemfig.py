"""Test native mol2chemfig (no Docker)"""
import subprocess
import os
import tempfile

def test_mol2chemfig_native():
    """Test if mol2chemfig works natively"""
    
    # Create a simple LaTeX file with mol2chemfig
    latex_content = r"""
\documentclass{article}
\usepackage{mol2chemfig}
\usepackage{chemfig}
\begin{document}
\chemfig{CCO}
\end{document}
"""
    
    # Create temp directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tex_file = os.path.join(tmpdir, "test.tex")
        
        # Write LaTeX file
        with open(tex_file, 'w') as f:
            f.write(latex_content)
        
        print(f"📄 Created test file: {tex_file}")
        
        # Run pdflatex
        print("🔨 Running pdflatex...")
        result = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', 'test.tex'],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            pdf_file = os.path.join(tmpdir, "test.pdf")
            if os.path.exists(pdf_file):
                print(f"✅ SUCCESS! PDF created: {pdf_file}")
                print(f"   Size: {os.path.getsize(pdf_file)} bytes")
                return True
            else:
                print("❌ FAILED: PDF not created")
                return False
        else:
            print("❌ FAILED: pdflatex error")
            print(result.stdout)
            print(result.stderr)
            return False

if __name__ == "__main__":
    print("Testing native mol2chemfig (no Docker)...\n")
    success = test_mol2chemfig_native()
    
    if success:
        print("\n🎉 mol2chemfig works natively!")
        print("✅ You can eliminate Docker completely")
    else:
        print("\n⚠️ mol2chemfig test failed")
        print("   Check if chemfig package is installed")
