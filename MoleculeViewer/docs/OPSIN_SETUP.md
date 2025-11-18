# OPSIN Setup Instructions

OPSIN (Open Parser for Systematic IUPAC Nomenclature) is required to convert chemical names to SMILES.

## Quick Setup

### Option 1: Automatic Download (Recommended)
```bash
python setup_opsin.py
```

### Option 2: Manual Download
1. Download the OPSIN JAR file:
   - From GitHub: https://github.com/dan2097/opsin/releases
   - Download: `opsin-cli-[version]-jar-with-dependencies.jar`
   
2. Save it to the project root directory as: `opsin-cli.jar`
   ```
   C:\Users\Kapil\Personal\PROJECTS\MoleculeViewer\opsin-cli.jar
   ```

3. Make sure you have Java 8 or higher installed:
   ```bash
   java -version
   ```

## Verify Installation

After downloading, test it:
```bash
cd C:\Users\Kapil\Personal\PROJECTS\MoleculeViewer
java -jar opsin-cli.jar -i "benzene"
```

You should see output like: `c1ccccc1`

## Troubleshooting

- **"Java not found"**: Install Java 8+ from https://www.java.com/
- **"JAR not found"**: Make sure opsin-cli.jar is in the project root
- **"Invalid SMILES"**: The compound name syntax might be incorrect

## Example IUPAC Names

These should work:
- `benzene` → `c1ccccc1`
- `acetic acid` → `CC(=O)O`
- `1-bromo-2,2-dimethylpropane` → `BrCC(C)(C)C`
- `aspirin` → `CC(=O)Oc1ccccc1C(=O)O`
- `caffeine` → `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
