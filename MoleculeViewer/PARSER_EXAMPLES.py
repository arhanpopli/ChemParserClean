"""
Parser Configuration Examples - Copy/Paste Ready

Edit app/config.py and copy one of these configurations.
"""

# ============================================================================
# CONFIGURATION 1: DEFAULT (Recommended)
# ============================================================================
# Try all parsers automatically, best coverage
"""
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = True
PARSER_TIMEOUT = 10
"""

# ============================================================================
# CONFIGURATION 2: OFFLINE MODE
# ============================================================================
# No internet needed, fast response
"""
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = False  # Skip internet lookup
PARSER_TIMEOUT = 5
"""

# ============================================================================
# CONFIGURATION 3: IUPAC ONLY
# ============================================================================
# Only accept systematic IUPAC names
"""
NOMENCLATURE_PARSER = 'opsin'
PARSER_TIMEOUT = 10
"""

# ============================================================================
# CONFIGURATION 4: MAXIMUM COVERAGE
# ============================================================================
# Always try PubChem for comprehensive results
"""
NOMENCLATURE_PARSER = 'pubchem'
PARSER_TIMEOUT = 15
"""

# ============================================================================
# CONFIGURATION 5: FAST & SIMPLE
# ============================================================================
# Common compounds only, instant results
"""
NOMENCLATURE_PARSER = 'chemdoodle'
"""

# ============================================================================
# CONFIGURATION 6: EDUCATION
# ============================================================================
# Common compounds + IUPAC, no internet
"""
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = False
PARSER_TIMEOUT = 5
"""

# ============================================================================
# CONFIGURATION 7: RESEARCH LAB
# ============================================================================
# Try everything, be patient
"""
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = True
PARSER_TIMEOUT = 20
"""

# ============================================================================
# CONFIGURATION 8: CHEMDOODLE + PUBCHEM
# ============================================================================
# Fast local lookup, then comprehensive online
"""
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = False
ENABLE_FALLBACK = False
ENABLE_PUBCHEM = True
PARSER_TIMEOUT = 15
"""

# ============================================================================
# HOW TO APPLY
# ============================================================================
"""
1. Choose a configuration above
2. Copy the lines (without triple quotes)
3. Open app/config.py
4. Find the "Nomenclature Parser Settings" section
5. Replace with your chosen configuration
6. Save the file
7. Stop server (Ctrl+C)
8. Start server: python run_server.py
9. Test: python test_parser_config.py
"""
