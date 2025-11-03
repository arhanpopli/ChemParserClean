#!/usr/bin/python
# -*- coding: utf-8 -*-

import tempfile
import os
import shutil
import base64

import pubchempy as pcp

from mol2chemfig.processor import process
from mol2chemfig import pdfgen
from mol2chemfig.pdfgen import update_pdf



def smiles_mol_to_chemfig(*args):
    error = None
    all_args = ''
    for arg in args:
        all_args += arg + ' '
    success, result = process(rawargs=all_args, progname='mol2chemfig')
    if success:
        pdfsuccess, pdfresult = pdfgen.pdfgen(result)
        if pdfsuccess:
            encoded = base64.encodestring(pdfresult)
            pdflink = "data:application/pdf;base64,{}".format(encoded)
            # pdflink = "data:image/png;base64,{}".format(encoded)

        else:
            pdflink = 'PDF file cannot be generated'
        try:
            chemfig = result.render_user()
            return chemfig, pdflink.replace('\n', ''), error
        except AttributeError:
            error = "Chemfig cannot be generated."
            return None, None, error
    else:
        error = "Sorry, Chemfig cannot be generated."
        return None, None, error


def convert_mol_format(mol_block, args=None):
    tempdir = tempfile.mkdtemp()
    path_to_file = os.path.join(tempdir, 'molecule.mol')
    with open(path_to_file, 'w') as f:
        f.write(mol_block)
    if args:
        chemfig, pdflink, error = smiles_mol_to_chemfig("-w "
                                                        + args,
                                                        path_to_file)
    else:
        chemfig, pdflink, error = smiles_mol_to_chemfig("-w", path_to_file)

    shutil.rmtree(tempdir)
    return chemfig, pdflink, error


def smiles_to_chemfig_image(*args):
    all_args = ''
    for i in args:
        all_args += i + ' '
    success, result = process(rawargs=all_args, progname='mol2chemfig')
    if success:
        image_success, image_result = pdfgen.image_gen(result)
        if image_success:
            encoded = base64.encodestring(image_result)
            # image_link = "data:image/png;base64,{}".format(encoded)
            image_link = "data:image/svg;base64,{}".format(encoded)

        else:
            image_link = 'image generation foobared'
    try:
        outcome = result.render_user()
        return outcome, image_link
    except AttributeError:
        error = "Chemfig cannot be generated."
        return None, error


def update_chemfig(data):
    pdfsuccess, pdfresult = update_pdf(data)
    if pdfsuccess:
        encoded = base64.encodestring(pdfresult)
        pdflink = "data:application/pdf;base64,{}".format(encoded)
    else:
        pdflink = 'pdf generation foobared'
    return pdflink


def get_smiles(name):
    error = None
    chemical_name = pcp.get_compounds(name, 'name')
    try:
        return error, chemical_name[0].isomeric_smiles
    except IndexError:
        error = "Compound {} was not found in the database."\
                .format(name.upper())
        return error, None


def combine_args(options, angle, indentation, h2):
    angle = " -a " + angle
    indentation = " -d " + indentation
    h2 = ' -y ' + h2
    return " ".join(options) + angle + indentation + h2


def chemfig_to_pdf(chemfig):
    tempdir = tempfile.mkdtemp()
    path_to_file = os.path.join(tempdir, 'reaction.tex')
    reaction_pdf = os.path.join(tempdir, 'reaction.pdf')
    latexcmd = "pdflatex -interaction=nonstopmode -output-directory=%s %s > /dev/null" %(tempdir, path_to_file)
    with open(path_to_file, 'w') as f:
        f.write(chemfig)
    os.system(latexcmd)
    with open(reaction_pdf, 'r') as f:
        pdf_string = f.read()
        encoded_pdf = base64.encodestring(pdf_string)
        pdf_link = "data:application/pdf;base64,{}".format(encoded_pdf)

    # pdf_string = open(reaction_pdf).read()
    # encoded_pdf = base64.encodestring(pdf_string)
    # pdf_link = "data:application/pdf;base64,{}".format(encoded_pdf)
    shutil.rmtree(tempdir)
    return pdf_link
