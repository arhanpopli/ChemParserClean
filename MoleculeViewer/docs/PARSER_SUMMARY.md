# Parser Configuration System - Summary

## What I've Added

A **manual parser selection system** that lets you choose which chemical name parser to use for nomenclature-to-SMILES conversion.

## Configuration Location

**File**: `app/config.py`  
**Section**: "Nomenclature Parser Settings"

## Key Settings

### Main Parser Selection
```python
NOMENCLATURE_PARSER = 'auto'  # or 'chemdoodle', 'opsin', 'pubchem', 'fallback'
```

### Fine-Tuning (for 'auto' mode)
```python
ENABLE_CHEMDOODLE = True   # ChemDoodle database (~300 compounds)
ENABLE_OPSIN = True        # OPSIN IUPAC parser (systematic names)
ENABLE_FALLBACK = True     # Internal dictionary (~50 compounds)
ENABLE_PUBCHEM = True      # PubChem API (millions of compounds)
```

### Timeout
```python
PARSER_TIMEOUT = 10  # seconds to wait for parser response
```

## Parser Options Explained

### 1. AUTO Mode (Default - Recommended)
**Setting**: `NOMENCLATURE_PARSER = 'auto'`

Tries all enabled parsers in order:
1. ChemDoodle Database (instant, common compounds)
2. OPSIN Parser (fast, IUPAC names)
3. Fallback Dictionary (instant, basic compounds)
4. PubChem API (comprehensive, requires internet)

**Best for**: Maximum coverage with automatic fallback

### 2. ChemDoodle Only
**Setting**: `NOMENCLATURE_PARSER = 'chemdoodle'`

- ✓ Instant lookup
- ✓ No dependencies
- ✗ Limited to ~300 compounds
- **Use when**: You only need common compounds (aspirin, caffeine, etc.)

### 3. OPSIN Only
**Setting**: `NOMENCLATURE_PARSER = 'opsin'`

- ✓ Works offline
- ✓ Excellent for IUPAC names
- ✗ Requires Java
- ✗ Doesn't handle common names
- **Use when**: You want only systematic IUPAC nomenclature

### 4. PubChem Only
**Setting**: `NOMENCLATURE_PARSER = 'pubchem'`

- ✓ Most comprehensive (millions of compounds)
- ✗ Requires internet
- ✗ Slower (API call)
- **Use when**: You need maximum coverage and have internet

### 5. Fallback Only
**Setting**: `NOMENCLATURE_PARSER = 'fallback'`

- ✓ Instant
- ✓ No dependencies
- ✗ Very limited (~50 compounds)
- **Use when**: You only need basic organic chemistry compounds

## How to Use

### Step 1: Edit Config
Open `app/config.py` and change:
```python
NOMENCLATURE_PARSER = 'opsin'  # Your choice
```

### Step 2: Restart Server
Stop (Ctrl+C) and restart:
```bash
python run_server.py
```

### Step 3: Test
```bash
python test_parser_config.py
```

## Common Scenarios

### Fast Offline Operation
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = False  # Skip slow API calls
```

### IUPAC Names Only
```python
NOMENCLATURE_PARSER = 'opsin'
```

### Maximum Coverage
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = True
PARSER_TIMEOUT = 15  # Be patient
```

### Quick Response
```python
NOMENCLATURE_PARSER = 'chemdoodle'
# or
NOMENCLATURE_PARSER = 'auto'
ENABLE_PUBCHEM = False
PARSER_TIMEOUT = 5
```

## Test Results

With default settings (`auto` mode, all parsers enabled):
- ✓ aspirin → ChemDoodle Database
- ✓ caffeine → ChemDoodle Database  
- ✓ glucose → OPSIN Parser
- ✓ 2-methylpropane → OPSIN Parser
- ✓ ethanol → OPSIN Parser
- ✓ benzene → ChemDoodle Database
- ✓ adrenaline → OPSIN Parser
- ✓ paracetamol → OPSIN Parser

**Result**: 8/8 success (100%)

## Files Created

1. **`app/config.py`** (updated)
   - Added `NOMENCLATURE_PARSER` setting
   - Added enable/disable flags
   - Added timeout setting

2. **`app/chemistry.py`** (updated)
   - Modified `nomenclature_to_smiles()` function
   - Added parser selection logic
   - Added `parser_override` parameter

3. **`PARSER_CONFIG_GUIDE.md`** 📖
   - Complete documentation
   - Use cases and examples
   - Troubleshooting guide

4. **`PARSER_QUICK_REF.md`** 📝
   - One-page quick reference
   - Common configurations

5. **`test_parser_config.py`** 🧪
   - Test script for your settings
   - Shows which compounds work

## Code Changes

### Before (Automatic Only)
```python
def nomenclature_to_smiles(compound_name):
    # Always tried all parsers in fixed order
    # No way to customize
```

### After (Configurable)
```python
def nomenclature_to_smiles(compound_name, parser_override=None):
    # Respects config.NOMENCLATURE_PARSER
    # Can override per-call
    # Can disable specific parsers
```

## Benefits

1. **Control**: Choose which parsers to use
2. **Speed**: Skip slow parsers (PubChem) when not needed
3. **Offline**: Disable internet-dependent parsers
4. **Focused**: Use only IUPAC parser for systematic names
5. **Flexible**: Fine-tune in AUTO mode

## Current Server Status

✅ **Running** at http://localhost:5000

**Current Configuration**:
- Mode: `auto` (try all parsers)
- All parsers enabled
- Timeout: 10 seconds
- Working perfectly!

## Quick Examples

### Example 1: School/Education (No Internet)
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = False  # No internet needed
```

### Example 2: Research Lab (Everything)
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = True
PARSER_TIMEOUT = 20  # Be patient for rare compounds
```

### Example 3: Production (Fast & Reliable)
```python
NOMENCLATURE_PARSER = 'chemdoodle'  # Known compounds only
```

### Example 4: IUPAC Standards
```python
NOMENCLATURE_PARSER = 'opsin'  # Systematic names only
```

## API Usage

The parser selection also works through the API:

```javascript
// Use configured parser
POST /api/nomenclature-to-smiles
{
  "nomenclature": "aspirin"
}

// Response includes source
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "nomenclature": "aspirin",
  "source": "ChemDoodle Database"
}
```

## Next Steps

1. ✅ Server running with parser config
2. ✅ Test script available
3. ✅ Documentation complete
4. **Try it**: Edit `app/config.py` and test different modes
5. **Test**: Run `python test_parser_config.py` after each change

---

**Remember**: Always **restart the server** after editing `app/config.py`!
