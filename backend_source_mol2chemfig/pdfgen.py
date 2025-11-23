'''
generate a pdf from a parsed mol2chemfig molecule.
return the result in a string.
'''
import os
import shutil
import subprocess
from tempfile import mkdtemp

latexfn = 'molecule.tex'
pdfname = 'molecule.pdf'
image_name = "molecule.jpg"
# imagecmd = "pdf2svg {} {}".format(pdfname, image_name)
# imagecmd = "convert -density 300 {} -quality 100 {}".format(pdfname, image_name)
# latexcmd = 'pdflatex -interaction=nonstopmode %s > /dev/null' % latexfn

m2pkg_path = '/usr/src/app/src/mol2chemfig'
#m2pkg_path = '/usr/src/app/backend/src/mol2chemfig'
pkg = '/mol2chemfig.sty'


def pdfgen(mol):
    tempdir = mkdtemp()
    try:
        os.symlink(m2pkg_path + pkg, tempdir + pkg)
    except OSError:
        pass

    chemfig = mol.render_server()
    width, height = mol.dimensions()

    atomsep = 16
    fixed_extra = 28
    width = round(atomsep * width) + fixed_extra
    height = round(atomsep * height) + fixed_extra
    global width_all, height_all
    width_all = width
    height_all = height
    latex = latex_template % locals()
    
    # Write latex file
    with open(os.path.join(tempdir, latexfn), 'w') as f:
        f.write(latex)

    # Run pdflatex in tempdir without changing global CWD
    cmd = ['pdflatex', '-interaction=nonstopmode', latexfn]
    
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call(cmd, cwd=tempdir, stdout=devnull, stderr=devnull)
            
        pdf_path = os.path.join(tempdir, pdfname)
        if os.path.exists(pdf_path):
            with open(pdf_path, 'rb') as f:
                pdfstring = f.read()
            return True, pdfstring
        else:
            return False, None

    except Exception as e:
        print("Error in pdfgen: " + str(e))
        return False, None
    finally:
        shutil.rmtree(tempdir)


def image_gen(mol):
    tempdir = mkdtemp()
    try:
        os.symlink(m2pkg_path + pkg, tempdir + pkg)
    except OSError:
        pass

    chemfig = mol.render_server()
    width, height = mol.dimensions()

    atomsep = 16
    fixed_extra = 28
    width = round(atomsep * width) + fixed_extra
    height = round(atomsep * height) + fixed_extra
    global width_all, height_all
    width_all = width
    height_all = height
    latex = latex_template % locals()
    
    with open(os.path.join(tempdir, latexfn), 'w') as f:
        f.write(latex)
        
    cmd_latex = ['pdflatex', '-interaction=nonstopmode', latexfn]
    # convert -density 300 molecule.pdf -quality 100 molecule.jpg
    cmd_convert = ['convert', '-density', '300', pdfname, '-quality', '100', image_name]

    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call(cmd_latex, cwd=tempdir, stdout=devnull, stderr=devnull)
            subprocess.call(cmd_convert, cwd=tempdir, stdout=devnull, stderr=devnull)
            
        img_path = os.path.join(tempdir, image_name)
        if os.path.exists(img_path):
            with open(img_path, 'rb') as f:
                image_string = f.read()
            return True, image_string
        else:
            return False, None
    except Exception as e:
        print("Error in image_gen: " + str(e))
        return False, None
    finally:
        shutil.rmtree(tempdir)


def update_pdf(mol):
    tempdir = mkdtemp()
    try:
        os.symlink(m2pkg_path + pkg, tempdir + pkg)
    except OSError:
        pass

    chemfig = "\chemfig {" + mol + "}"
    atomsep = 16
    fixed_extra = 28
    width = width_all
    height = height_all
    latex = latex_template % locals()

    with open(os.path.join(tempdir, latexfn), 'w') as f:
        f.write(latex)
        
    cmd = ['pdflatex', '-interaction=nonstopmode', latexfn]

    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call(cmd, cwd=tempdir, stdout=devnull, stderr=devnull)
            
        pdf_path = os.path.join(tempdir, pdfname)
        if os.path.exists(pdf_path):
            with open(pdf_path, 'rb') as f:
                pdfstring = f.read()
            return True, pdfstring
        else:
            return False, None
    except Exception as e:
        print("Error in update_pdf: " + str(e))
        return False, None
    finally:
        shutil.rmtree(tempdir)


latex_template = r'''
\documentclass{minimal}
\usepackage{xcolor, mol2chemfig}
\usepackage[margin=(margin)spt,papersize={%(width)spt, %(height)spt}]{geometry}

\usepackage[helvet]{sfmath}
\setcrambond{2.5pt}{0.4pt}{1.0pt}
\setbondoffset{1pt}
\setdoublesep{3pt}
\setatomsep{%(atomsep)spt}
\renewcommand{\printatom}[1]{\fontsize{8pt}{10pt}\selectfont{\ensuremath{\mathsf{#1}}}}

\setlength{\parindent}{0pt}
\setlength{\fboxsep}{0pt}
\begin{document}
\vspace*{\fill}
\vspace{-8pt}
\begin{center}
%(chemfig)s
\end{center}
\vspace*{\fill}
\end{document}
'''
