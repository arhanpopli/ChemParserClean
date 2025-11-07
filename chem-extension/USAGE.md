# Chemistry Renderer - Extension Usage Guide

## 📌 Updated Syntax: `chem:text:`

The extension now uses **double colons** to clearly define what gets converted to a molecule.

### Syntax
```
chem:TEXT:
```

Everything between the two colons is sent to MoleculeViewer.

### Examples

#### ✅ Chemical Names
```
chem:benzene:
chem:acetone:
chem:1-chloro-benzene:
chem:aspirin:
chem:caffeine:
```

#### ✅ SMILES Notation
```
chem:CCO:          (ethanol)
chem:CC(=O)C:      (acetone)
chem:c1ccccc1:     (benzene)
chem:C1=CC=CC=C1:  (benzene alternative)
```

#### ✅ Complex Names
```
chem:2-methylpropane:
chem:n-butyl acetate:
chem:1,2-dichlorobenzene:
chem:alpha-pinene:
```

## 🖥️ Starting the System

### Windows Users
1. **Double-click `START_ALL.bat`** in the MoleculeViewer folder
2. A black command window will open (leave it running)
3. Your browser will automatically open to http://localhost:5000
4. You can now use the extension in ChatGPT

### Manual Start (if needed)
```bash
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer
python -m flask --app app.api run --host=0.0.0.0 --port=5000
```

Then go to: http://localhost:5000

## 🧪 Testing in ChatGPT

1. Open ChatGPT
2. Type: `Look up this molecule: chem:aspirin:`
3. The extension will replace it with a rendered molecule image
4. Or try: `Show me chem:1-methyl-naphthalene:`

## 🎨 Features

The extension automatically detects if your input is:
- **Chemical name** (like `aspirin`, `caffeine`) → Looks up SMILES
- **SMILES string** (like `CC(=O)O`, `c1ccccc1`) → Renders directly

## ⚙️ Troubleshooting

**"It's not rendering"**
- Make sure the black command window from START_ALL.bat is still open
- Reload the ChatGPT page (F5)
- Check that you're using the correct syntax: `chem:NAME:`

**"I see the text but not the image"**
- Reload the ChatGPT page
- Try a different molecule like `chem:benzene:`
- Open the browser console (F12) to check for errors

**"Port 5000 is already in use"**
- Close other Flask instances or use a different port
- Edit START_ALL.bat to use a different port number

## 📝 Notes

- Everything between `chem:` and `:` is processed
- Names can include numbers, hyphens, and spaces
- SMILES strings must use valid SMILES notation
- Molecules are cached for 24 hours on the server
- Works with ChatGPT and any website with the extension enabled

---

**Quick Start:** Just double-click `START_ALL.bat` and you're ready to go!
