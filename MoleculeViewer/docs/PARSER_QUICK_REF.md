# Parser Configuration - Quick Reference

## How to Change Parser

**Edit**: `app/config.py`  
**Find**: `NOMENCLATURE_PARSER`  
**Restart**: Server after changes

## Options

### Try All (Recommended)
```python
NOMENCLATURE_PARSER = 'auto'
```
Best coverage - tries all parsers automatically.

### ChemDoodle Only
```python
NOMENCLATURE_PARSER = 'chemdoodle'
```
Fast, common compounds only (~300).

### OPSIN Only
```python
NOMENCLATURE_PARSER = 'opsin'
```
IUPAC names only. Requires Java.

### PubChem Only
```python
NOMENCLATURE_PARSER = 'pubchem'
```
Most comprehensive. Requires internet.

### Fallback Only
```python
NOMENCLATURE_PARSER = 'fallback'
```
Basic compounds only (~50).

## Disable Specific Parsers

```python
NOMENCLATURE_PARSER = 'auto'

# Skip PubChem (faster, offline)
ENABLE_PUBCHEM = False

# Skip OPSIN (no Java needed)
ENABLE_OPSIN = False
```

## Change Timeout

```python
PARSER_TIMEOUT = 5  # Faster response (default: 10)
PARSER_TIMEOUT = 20  # More patience for slow connections
```

## Test Your Settings

```bash
python test_parser_config.py
```

Shows which compounds work with your configuration.

---

**Remember**: Restart server after editing config!
