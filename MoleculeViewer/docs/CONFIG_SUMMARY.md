# CH₃ Configuration System - Summary

## What I've Created

I've implemented a **configurable system** for CH₃ label sizing and positioning that automatically scales with molecule size.

## Files Created/Modified

### 1. `app/config.py` ⭐ **MAIN CONFIG FILE**
This is where you adjust CH₃ label settings. Key parameters:

- **`CARBON_LABEL_FONT_SIZE`** = 32 (base size in pixels)
- **`CARBON_LABEL_SCALING`** = 'auto' or 'fixed'
- **`CARBON_LABEL_SCALE_FACTOR`** = 0.55 (for auto-scaling)
- **`CARBON_LABEL_MIN_SIZE`** = 18 (minimum when auto-scaling)
- **`CARBON_LABEL_MAX_SIZE`** = 48 (maximum when auto-scaling)
- **`CARBON_LABEL_OFFSET_X`** = 8 (horizontal position adjustment)
- **`CARBON_LABEL_OFFSET_Y`** = 14 (vertical position adjustment)
- **`CARBON_LABEL_SCALE_OFFSETS`** = True (scale position with size)

### 2. `app/chemistry.py` (Modified)
Updated to:
- Calculate average bond length for each molecule
- Scale CH₃ labels based on bond length
- Use config settings for all label properties

### 3. `CONFIG_GUIDE.md` 📖
Complete documentation with:
- Quick start instructions
- Common adjustments
- Recommended settings
- Troubleshooting tips

### 4. `CONFIG_EXAMPLES.py` 📝
Ready-to-use configuration examples:
- Fixed size configurations
- Large label settings
- Small label settings
- Different use cases

### 5. `test_config_scaling.py` 🧪
Test script that:
- Shows current config settings
- Tests multiple molecule sizes
- Generates SVG files for inspection
- Displays font sizes used

## How It Works

### Auto-Scaling Mode (Default)
1. System calculates average bond length in molecule
2. Font size = bond_length × SCALE_FACTOR
3. Size is clamped between MIN and MAX
4. Offsets also scale proportionally

**Result**: Small molecules get appropriately sized labels, large molecules get proportionally sized labels.

### Fixed Mode
- Always uses CARBON_LABEL_FONT_SIZE
- Same size for all molecules
- Good for standardized output

## Quick Start

### To Change Size:

**Option 1: Edit `app/config.py`**
```python
# Make bigger:
CARBON_LABEL_SCALE_FACTOR = 0.65  # (default: 0.55)

# Make smaller:
CARBON_LABEL_SCALE_FACTOR = 0.45  # (default: 0.55)
```

**Option 2: Use Fixed Size**
```python
CARBON_LABEL_SCALING = 'fixed'
CARBON_LABEL_FONT_SIZE = 28  # Your preferred size
```

### To Adjust Position:
```python
CARBON_LABEL_OFFSET_X = 10  # Move right (default: 8)
CARBON_LABEL_OFFSET_Y = 16  # Move down (default: 14)
```

### After Changes:
1. Save `app/config.py`
2. Stop server (Ctrl+C)
3. Restart: `python run_server.py`
4. Test at http://localhost:5000

## Testing Your Changes

```bash
# Test different molecule sizes:
python test_config_scaling.py

# Check generated SVG files in the directory
```

## Current Behavior

With current settings (`CARBON_LABEL_SCALING = 'auto'`, `SCALE_FACTOR = 0.55`):

- **Small molecules** (propane): ~41px font size
- **Medium molecules** (octane): ~41px font size  
- **Large molecules**: Adjusts within 18-48px range
- **Position**: Scales proportionally with size

## Common Adjustments

| Problem | Solution |
|---------|----------|
| Labels too big on small molecules | Decrease `SCALE_FACTOR` to 0.50 or lower `MAX_SIZE` to 42 |
| Labels too small on large molecules | Increase `SCALE_FACTOR` to 0.60 or raise `MIN_SIZE` to 22 |
| Labels positioned wrong | Adjust `OFFSET_X` and `OFFSET_Y` values |
| Want same size everywhere | Set `CARBON_LABEL_SCALING = 'fixed'` |
| Labels vary too much | Narrow the `MIN_SIZE` to `MAX_SIZE` range |

## Files to Reference

1. **`CONFIG_GUIDE.md`** - Complete documentation
2. **`CONFIG_EXAMPLES.py`** - Copy-paste examples
3. **`app/config.py`** - Where to make changes
4. **`test_config_scaling.py`** - Test your settings

## Example Workflow

1. Open `CONFIG_EXAMPLES.py`
2. Find an example that matches your needs
3. Copy those settings
4. Open `app/config.py`
5. Paste settings (replace existing lines)
6. Save file
7. Restart server
8. Test in browser
9. If not perfect, tweak and repeat

## Key Insight

The system now **calculates actual bond lengths** from the molecule's 2D coordinates and scales labels proportionally. This means:

- ✅ Tiny molecules don't have gigantic labels
- ✅ Large molecules don't have tiny labels  
- ✅ Everything stays proportional
- ✅ You can still override with fixed size if needed

## Server Running

Currently running at: **http://localhost:5000**

With configuration:
- Mode: Auto-scaling
- Factor: 0.55
- Range: 18-48px
- Offsets: X=8, Y=14

---

**Remember**: Configuration changes require server restart!
