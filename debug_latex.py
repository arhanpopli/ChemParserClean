import subprocess
import os

print("Compiling failed_mol.tex...")
proc = subprocess.run(
    ['pdflatex', '-interaction=nonstopmode', 'failed_mol.tex'],
    cwd=os.path.dirname(__file__),
    capture_output=True,
    text=True
)

print("Return code:", proc.returncode)
print("STDOUT (last 1000 chars):")
print(proc.stdout[-1000:])
