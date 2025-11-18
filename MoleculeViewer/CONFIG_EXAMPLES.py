"""
QUICK CONFIG EXAMPLES - Copy/paste these into app/config.py

Choose one of these configurations based on your needs.
After editing config.py, restart the server for changes to take effect.
"""

# ============================================================================
# EXAMPLE 1: DEFAULT (Current Settings)
# ============================================================================
# Good balance for most molecules
# Labels scale automatically with molecule size
"""
CARBON_LABEL_FONT_SIZE = 32
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.55
CARBON_LABEL_MIN_SIZE = 18
CARBON_LABEL_MAX_SIZE = 48
CARBON_LABEL_OFFSET_X = 8
CARBON_LABEL_OFFSET_Y = 14
CARBON_LABEL_SCALE_OFFSETS = True
"""

# ============================================================================
# EXAMPLE 2: FIXED SIZE - Always 28px
# ============================================================================
# Use this if you want consistent size regardless of molecule
# Good for standardized output
"""
CARBON_LABEL_FONT_SIZE = 28
CARBON_LABEL_SCALING = 'fixed'
CARBON_LABEL_OFFSET_X = 8
CARBON_LABEL_OFFSET_Y = 14
"""

# ============================================================================
# EXAMPLE 3: LARGER LABELS
# ============================================================================
# For presentations or when working with complex molecules
# Makes CH₃ more visible
"""
CARBON_LABEL_FONT_SIZE = 36
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.65
CARBON_LABEL_MIN_SIZE = 24
CARBON_LABEL_MAX_SIZE = 54
CARBON_LABEL_OFFSET_X = 9
CARBON_LABEL_OFFSET_Y = 16
CARBON_LABEL_SCALE_OFFSETS = True
"""

# ============================================================================
# EXAMPLE 4: SMALLER LABELS
# ============================================================================
# For densely packed molecules or to match smaller functional groups
"""
CARBON_LABEL_FONT_SIZE = 24
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.45
CARBON_LABEL_MIN_SIZE = 16
CARBON_LABEL_MAX_SIZE = 36
CARBON_LABEL_OFFSET_X = 7
CARBON_LABEL_OFFSET_Y = 12
CARBON_LABEL_SCALE_OFFSETS = True
"""

# ============================================================================
# EXAMPLE 5: FIXED SIZE WITH ADJUSTED POSITION
# ============================================================================
# If labels appear slightly off-center, adjust offsets
"""
CARBON_LABEL_FONT_SIZE = 30
CARBON_LABEL_SCALING = 'fixed'
CARBON_LABEL_OFFSET_X = 10  # Move 2px right from default
CARBON_LABEL_OFFSET_Y = 16  # Move 2px down from default
CARBON_LABEL_SCALE_OFFSETS = False
"""

# ============================================================================
# EXAMPLE 6: MATCH RDKIT FUNCTIONAL GROUPS EXACTLY
# ============================================================================
# Try to match the size of OH, NO2, etc. as closely as possible
"""
CARBON_LABEL_FONT_SIZE = 30
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.52
CARBON_LABEL_MIN_SIZE = 20
CARBON_LABEL_MAX_SIZE = 44
CARBON_LABEL_OFFSET_X = 8
CARBON_LABEL_OFFSET_Y = 14
CARBON_LABEL_SCALE_OFFSETS = True
"""

# ============================================================================
# EXAMPLE 7: VERY SMALL MOLECULES (like propane, butane)
# ============================================================================
# Optimized for small chain molecules in large canvas
"""
CARBON_LABEL_FONT_SIZE = 32
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.50  # Slightly smaller scaling
CARBON_LABEL_MIN_SIZE = 24  # Higher minimum
CARBON_LABEL_MAX_SIZE = 48
CARBON_LABEL_OFFSET_X = 8
CARBON_LABEL_OFFSET_Y = 14
CARBON_LABEL_SCALE_OFFSETS = True
"""

# ============================================================================
# EXAMPLE 8: LARGE COMPLEX MOLECULES
# ============================================================================
# Optimized for steroids, terpenes, etc.
"""
CARBON_LABEL_FONT_SIZE = 28
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.48  # Smaller scaling
CARBON_LABEL_MIN_SIZE = 16  # Lower minimum
CARBON_LABEL_MAX_SIZE = 40
CARBON_LABEL_OFFSET_X = 7
CARBON_LABEL_OFFSET_Y = 13
CARBON_LABEL_SCALE_OFFSETS = True
"""

# ============================================================================
# HOW TO USE
# ============================================================================
"""
1. Choose one example above that fits your needs
2. Copy the settings (without the triple quotes)
3. Open app/config.py
4. Replace the corresponding lines in the config file
5. Save the file
6. Stop the server (Ctrl+C)
7. Start the server again: python run_server.py
8. Test at http://localhost:5000

If the first try isn't perfect, adjust the values slightly and restart again.

Key parameters to tweak:
- CARBON_LABEL_SCALE_FACTOR: Controls overall size (0.4 = smaller, 0.7 = larger)
- CARBON_LABEL_FONT_SIZE: Base size for fixed mode
- CARBON_LABEL_OFFSET_X/Y: Position fine-tuning
"""
