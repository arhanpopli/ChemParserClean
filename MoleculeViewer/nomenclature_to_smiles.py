#!/usr/bin/env python3
"""
Convert chemical nomenclature to SMILES
Called by Node.js server via subprocess
"""

import sys
import json

try:
    import requests
    from urllib.parse import quote

    def nomenclature_to_smiles(nomenclature):
        """Convert chemical name to SMILES using PubChem API"""
        try:
            nomenclature = nomenclature.strip()
            if not nomenclature:
                return {"error": "Empty nomenclature"}

            # Use PubChem API to convert name to SMILES
            # GET /cactus/chemical/name/{name}/smiles
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(nomenclature)}/property/CanonicalSMILES/JSON"

            response = requests.get(url, timeout=10)
            response.raise_for_status()

            data = response.json()

            # Extract SMILES from response - PubChem uses PropertyTable.Properties
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                props = data['PropertyTable']['Properties']
                if len(props) > 0:
                    # Try different SMILES field names (PubChem API variations)
                    smiles = props[0].get('CanonicalSMILES') or props[0].get('ConnectivitySMILES') or props[0].get('IsomericSMILES')
                    if smiles:
                        return {"smiles": smiles}

            return {"error": f"Could not find SMILES for: {nomenclature}"}

        except requests.exceptions.RequestException as e:
            return {"error": f"API error: {str(e)}"}
        except (KeyError, IndexError) as e:
            return {"error": f"Invalid response: {str(e)}"}

    # Main execution
    if __name__ == "__main__":
        nomenclature = sys.argv[1] if len(sys.argv) > 1 else ""
        result = nomenclature_to_smiles(nomenclature)
        print(json.dumps(result))

except ImportError:
    print(json.dumps({"error": "requests library not installed. Run: pip install requests"}))
    sys.exit(1)
