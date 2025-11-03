"""
Configuration for molecule rendering options.

This file contains adjustable parameters for fine-tuning the appearance
of chemical structures.
"""

# ============================================================================
# Carbon Label Settings (CH₃, CH₂, CH, C)
# ============================================================================

# Base font size for carbon labels (in pixels)
# This is the font size used when rendering at standard molecule size
# Adjust this value to make CH₃ labels larger or smaller
# Recommended range: 24-40
CARBON_LABEL_FONT_SIZE = 32

# Font scaling mode for carbon labels
# Options:
#   'fixed' - Always use CARBON_LABEL_FONT_SIZE regardless of molecule size
#   'auto'  - Scale font size based on average bond length (recommended)
CARBON_LABEL_SCALING = 'auto'

# When using 'auto' scaling, this multiplier adjusts the font size
# relative to the average bond length in the molecule
# Formula: font_size = avg_bond_length * CARBON_LABEL_SCALE_FACTOR
# Increase for larger labels, decrease for smaller labels
# Recommended range: 0.4 - 0.7
CARBON_LABEL_SCALE_FACTOR = 0.55

# Minimum and maximum font sizes when using auto scaling (in pixels)
# This prevents labels from becoming too small or too large
CARBON_LABEL_MIN_SIZE = 18
CARBON_LABEL_MAX_SIZE = 48

# ============================================================================
# Carbon Label Position Adjustments
# ============================================================================

# X-axis offset from the dummy atom 'E' position (in pixels)
# Positive values move the label right, negative values move it left
# This fine-tunes horizontal positioning
CARBON_LABEL_OFFSET_X = 8

# Y-axis offset from the dummy atom 'E' position (in pixels)
# Positive values move the label down, negative values move it up
# This fine-tunes vertical positioning (adjust for baseline alignment)
CARBON_LABEL_OFFSET_Y = 14

# When using auto scaling, also scale the position offsets
# This maintains consistent positioning relative to label size
CARBON_LABEL_SCALE_OFFSETS = True

# ============================================================================
# Aromatic Circle Settings
# ============================================================================

# Radius scaling factors for aromatic circles (as percentage of inradius)
# Smaller values = smaller circles (more space between circle and bonds)
# Larger values = larger circles (closer to bonds)
AROMATIC_CIRCLE_RADIUS_6 = 0.70  # For 6-membered rings (benzene)
AROMATIC_CIRCLE_RADIUS_5 = 0.68  # For 5-membered rings (pyrrole, etc.)
AROMATIC_CIRCLE_RADIUS_OTHER = 0.63  # For other ring sizes

# ============================================================================
# Advanced Options
# ============================================================================

# Font family for carbon labels
CARBON_LABEL_FONT_FAMILY = 'sans-serif'

# Text anchor for carbon labels (alignment)
# Options: 'start', 'middle', 'end'
CARBON_LABEL_TEXT_ANCHOR = 'middle'

# Color for carbon labels (hex color or named color)
CARBON_LABEL_COLOR = '#191919'  # Black (matches RDKit's dummy atoms)

# ============================================================================
# Nomenclature Parser Settings
# ============================================================================

# Default parser priority for chemical nomenclature to SMILES conversion
# Options:
#   'auto' - Try all parsers in order: ChemDoodle → OPSIN → Fallback → PubChem
#   'chemdoodle' - Use ONLY ChemDoodle database (fast, limited coverage)
#   'opsin' - Use ONLY OPSIN parser (IUPAC names, requires Java)
#   'pubchem' - Use ONLY PubChem API (comprehensive, requires internet)
#   'fallback' - Use ONLY internal fallback dictionary (basic compounds)
NOMENCLATURE_PARSER = 'auto'

# Timeout for parser operations (in seconds)
PARSER_TIMEOUT = 10

# Enable/disable specific parsers when using 'auto' mode
# Set to False to skip a specific parser in auto mode
ENABLE_CHEMDOODLE = True  # ChemDoodle database lookup
ENABLE_OPSIN = True       # OPSIN IUPAC parser
ENABLE_FALLBACK = True    # Internal fallback dictionary
ENABLE_PUBCHEM = True     # PubChem API lookup
