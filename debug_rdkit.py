from rdkit.Chem.Draw import rdMolDraw2D
d = rdMolDraw2D.MolDraw2DSVG(200, 200)
opts = d.drawOptions()
print(dir(opts))
