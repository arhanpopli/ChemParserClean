# CH₃ Label Configuration Guide

This guide explains how to customize the size and positioning of CH₃ labels in your molecule viewer.

## Quick Start

All configuration is in: **`app/config.py`**

Edit this file to customize CH₃ label appearance.

## Common Adjustments

### Make CH₃ Labels Bigger
```python
# Option 1: Increase fixed size
CARBON_LABEL_FONT_SIZE = 36  # Default: 32

# Option 2: Increase scale factor (for auto-scaling)
CARBON_LABEL_SCALE_FACTOR = 0.65  # Default: 0.55
```

### Make CH₃ Labels Smaller
```python
# Option 1: Decrease fixed size
CARBON_LABEL_FONT_SIZE = 24  # Default: 32

# Option 2: Decrease scale factor (for auto-scaling)
CARBON_LABEL_SCALE_FACTOR = 0.45  # Default: 0.55
```

### Use Fixed Size (No Auto-Scaling)
```python
CARBON_LABEL_SCALING = 'fixed'  # Change from 'auto'
CARBON_LABEL_FONT_SIZE = 28  # Set your preferred size
```

### Adjust Position (if labels appear off-center)
```python
# Move label right (+) or left (-)
CARBON_LABEL_OFFSET_X = 10  # Default: 8

# Move label down (+) or up (-)
CARBON_LABEL_OFFSET_Y = 16  # Default: 14

# Should offsets scale with font size?
CARBON_LABEL_SCALE_OFFSETS = True  # Default: True
```

### Set Size Limits for Auto-Scaling
```python
CARBON_LABEL_MIN_SIZE = 20  # Smallest allowed (default: 18)
CARBON_LABEL_MAX_SIZE = 40  # Largest allowed (default: 48)
```

## Configuration Parameters

### Font Size Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CARBON_LABEL_FONT_SIZE` | 32 | Base font size in pixels |
| `CARBON_LABEL_SCALING` | 'auto' | 'auto' or 'fixed' |
| `CARBON_LABEL_SCALE_FACTOR` | 0.55 | Multiplier for auto-scaling |
| `CARBON_LABEL_MIN_SIZE` | 18 | Minimum font size (auto mode) |
| `CARBON_LABEL_MAX_SIZE` | 48 | Maximum font size (auto mode) |

### Position Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CARBON_LABEL_OFFSET_X` | 8 | Horizontal offset in pixels |
| `CARBON_LABEL_OFFSET_Y` | 14 | Vertical offset in pixels |
| `CARBON_LABEL_SCALE_OFFSETS` | True | Scale offsets with font size |

### Style Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CARBON_LABEL_FONT_FAMILY` | 'sans-serif' | Font family |
| `CARBON_LABEL_TEXT_ANCHOR` | 'middle' | Text alignment |
| `CARBON_LABEL_COLOR` | '#191919' | Text color (hex) |

## How Auto-Scaling Works

When `CARBON_LABEL_SCALING = 'auto'`:

1. System calculates average bond length in the molecule
2. Font size = bond_length × `CARBON_LABEL_SCALE_FACTOR`
3. Font size is clamped between `MIN_SIZE` and `MAX_SIZE`
4. If `SCALE_OFFSETS = True`, position offsets also scale proportionally

This ensures CH₃ labels look proportional regardless of molecule size.

## Recommended Settings

### For Most Molecules (Balanced)
```python
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.55
CARBON_LABEL_MIN_SIZE = 18
CARBON_LABEL_MAX_SIZE = 48
```

### For Large Complex Molecules
```python
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.50  # Slightly smaller
CARBON_LABEL_MIN_SIZE = 16
CARBON_LABEL_MAX_SIZE = 42
```

### For Small Simple Molecules
```python
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.60  # Slightly larger
CARBON_LABEL_MIN_SIZE = 20
CARBON_LABEL_MAX_SIZE = 50
```

### For Consistent Size Across All Molecules
```python
CARBON_LABEL_SCALING = 'fixed'
CARBON_LABEL_FONT_SIZE = 30  # Same size everywhere
```

## Testing Your Changes

After editing `app/config.py`:

1. **Stop the server** (Ctrl+C)
2. **Restart the server**:
   ```bash
   python run_server.py
   ```
3. **Test in browser**: http://localhost:5000

The configuration is loaded when the server starts, so restart is required for changes to take effect.

### Test Script

Run this to see how your settings affect different molecule sizes:
```bash
python test_config_scaling.py
```

This generates SVG files you can inspect to verify sizing.

## Troubleshooting

### Labels Too Big on Small Molecules
- Decrease `CARBON_LABEL_MAX_SIZE` to 36 or 40
- Or decrease `CARBON_LABEL_SCALE_FACTOR` to 0.50

### Labels Too Small on Large Molecules
- Increase `CARBON_LABEL_MIN_SIZE` to 22 or 24
- Or increase `CARBON_LABEL_SCALE_FACTOR` to 0.60

### Labels Positioned Incorrectly
- Adjust `CARBON_LABEL_OFFSET_X` and `CARBON_LABEL_OFFSET_Y`
- Try setting `CARBON_LABEL_SCALE_OFFSETS = False` for consistent positioning

### Labels Don't Scale with Molecule Size
- Ensure `CARBON_LABEL_SCALING = 'auto'`
- Check that you restarted the server after config changes

## Examples

### Configuration A: Always 28px
```python
CARBON_LABEL_SCALING = 'fixed'
CARBON_LABEL_FONT_SIZE = 28
```
Result: All CH₃ labels exactly 28px, regardless of molecule size.

### Configuration B: Proportional Scaling
```python
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.55
CARBON_LABEL_MIN_SIZE = 18
CARBON_LABEL_MAX_SIZE = 48
```
Result: CH₃ labels scale with molecule, staying between 18-48px.

### Configuration C: Larger Labels
```python
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.65
CARBON_LABEL_MIN_SIZE = 24
CARBON_LABEL_MAX_SIZE = 54
```
Result: Noticeably larger CH₃ labels, good for presentations.

### Configuration D: Fixed Position Fine-Tuning
```python
CARBON_LABEL_OFFSET_X = 9
CARBON_LABEL_OFFSET_Y = 15
CARBON_LABEL_SCALE_OFFSETS = False  # Don't scale offsets
```
Result: Labels positioned slightly different, position doesn't change with size.

---

**Note**: After any configuration change, you must **restart the server** for changes to take effect.
