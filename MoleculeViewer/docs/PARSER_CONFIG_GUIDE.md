# Nomenclature Parser Configuration Guide

This guide explains how to manually select which chemical name parser to use.

## Quick Start

Edit **`app/config.py`** and set the `NOMENCLATURE_PARSER` option.

## Parser Options

### 1. AUTO Mode (Default - Recommended)
```python
NOMENCLATURE_PARSER = 'auto'
```
- Tries all parsers in order: ChemDoodle → OPSIN → Fallback → PubChem
- Best coverage, automatically finds the compound
- Use this unless you have a specific reason to limit parsers

### 2. ChemDoodle Only
```python
NOMENCLATURE_PARSER = 'chemdoodle'
```
- Uses only the ChemDoodle database
- **Pros**: Instant lookup, no external dependencies
- **Cons**: Limited to ~300 common compounds
- **Best for**: Common compounds (aspirin, caffeine, glucose, etc.)

### 3. OPSIN Only
```python
NOMENCLATURE_PARSER = 'opsin'
```
- Uses only the OPSIN IUPAC parser
- **Pros**: Excellent for systematic IUPAC names, works offline
- **Cons**: Requires Java, doesn't handle common names
- **Best for**: IUPAC systematic names (2-methylpropane, benzene, etc.)
- **Requires**: Java installed and `opsin-cli.jar` file

### 4. PubChem Only
```python
NOMENCLATURE_PARSER = 'pubchem'
```
- Uses only the PubChem API
- **Pros**: Most comprehensive database (millions of compounds)
- **Cons**: Requires internet connection, slower (API call)
- **Best for**: Obscure compounds, when offline parsers fail
- **Requires**: Internet connection

### 5. Fallback Only
```python
NOMENCLATURE_PARSER = 'fallback'
```
- Uses only the internal fallback dictionary
- **Pros**: Instant, no dependencies
- **Cons**: Very limited (~50 compounds)
- **Best for**: Basic organic chemistry compounds

## Fine-Tuning AUTO Mode

You can disable specific parsers in AUTO mode:

```python
NOMENCLATURE_PARSER = 'auto'

# Disable specific parsers
ENABLE_CHEMDOODLE = True   # Set to False to skip ChemDoodle
ENABLE_OPSIN = True        # Set to False to skip OPSIN
ENABLE_FALLBACK = True     # Set to False to skip fallback dictionary
ENABLE_PUBCHEM = False     # Set to False to skip PubChem (saves time)
```

**Example**: Skip PubChem for faster offline operation:
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_PUBCHEM = False  # Don't wait for API calls
```

## Timeout Setting

Control how long to wait for parsers:

```python
PARSER_TIMEOUT = 10  # seconds (default)
```

- Increase for slow connections: `PARSER_TIMEOUT = 20`
- Decrease for faster response: `PARSER_TIMEOUT = 5`

## Parser Comparison

| Parser | Coverage | Speed | Offline | Requirements |
|--------|----------|-------|---------|--------------|
| **ChemDoodle** | ~300 compounds | Instant | ✓ Yes | None |
| **OPSIN** | IUPAC names | Fast | ✓ Yes | Java |
| **Fallback** | ~50 compounds | Instant | ✓ Yes | None |
| **PubChem** | Millions | Slow | ✗ No | Internet |

## Use Cases

### Research Lab (Comprehensive)
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = True
PARSER_TIMEOUT = 15
```
Try everything, wait if needed.

### Education (Fast, Offline)
```python
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = True
ENABLE_FALLBACK = True
ENABLE_PUBCHEM = False  # No internet needed
PARSER_TIMEOUT = 5
```
Common compounds and IUPAC names only.

### IUPAC Focus
```python
NOMENCLATURE_PARSER = 'opsin'
```
Only accept systematic IUPAC nomenclature.

### Maximum Coverage
```python
NOMENCLATURE_PARSER = 'pubchem'
```
Always query PubChem (most comprehensive).

### Specific Workflow
```python
# Try ChemDoodle first, then PubChem only
NOMENCLATURE_PARSER = 'auto'
ENABLE_CHEMDOODLE = True
ENABLE_OPSIN = False
ENABLE_FALLBACK = False
ENABLE_PUBCHEM = True
```

## Examples

### Example 1: Common Name
**Input**: "aspirin"

- **ChemDoodle**: ✓ Found
- **OPSIN**: ✗ Failed (not IUPAC)
- **Fallback**: ✓ Found
- **PubChem**: ✓ Found

**Result in AUTO mode**: ChemDoodle (first to match)

### Example 2: IUPAC Name
**Input**: "2-acetoxybenzoic acid"

- **ChemDoodle**: ✗ Not in database
- **OPSIN**: ✓ Found
- **Fallback**: ✗ Not in dictionary
- **PubChem**: ✓ Found

**Result in AUTO mode**: OPSIN (second tier)

### Example 3: Obscure Compound
**Input**: "rapamycin"

- **ChemDoodle**: ✗ Not in database
- **OPSIN**: ✗ Failed (not IUPAC)
- **Fallback**: ✗ Not in dictionary
- **PubChem**: ✓ Found

**Result in AUTO mode**: PubChem (last resort)

## Applying Changes

1. Edit `app/config.py`
2. Change `NOMENCLATURE_PARSER` setting
3. **Restart the server** (changes require restart)
4. Test at http://localhost:5000

## Testing Your Configuration

Use the test script:
```bash
python test_parser_config.py
```

This will test your parser settings with various compound names.

## Troubleshooting

### "OPSIN parser could not convert"
- OPSIN requires Java: Install Java JRE
- Check that `opsin-cli.jar` exists in the project root
- Try a proper IUPAC name instead of common name

### "PubChem could not find"
- Check internet connection
- Verify compound name spelling
- Try a more common synonym

### Slow Response
- Decrease `PARSER_TIMEOUT` to 5 seconds
- Set `ENABLE_PUBCHEM = False` to skip API calls
- Use `'chemdoodle'` mode for instant results

### Nothing Works
- Use `NOMENCLATURE_PARSER = 'auto'` (try all parsers)
- Check that Java is installed (for OPSIN)
- Verify internet connection (for PubChem)

---

**Remember**: Configuration changes require **server restart** to take effect!
