"""
SVG Cache Manager - Handles intelligent caching with meaningful filenames
"""

import os
import hashlib
import time
from pathlib import Path
from typing import Dict, Tuple, Optional
import json

# Cache configuration
CACHE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'svg-cache'))
CACHE_EXPIRY_HOURS = 24

# Create cache directory if it doesn't exist
os.makedirs(CACHE_DIR, exist_ok=True)


def sanitize_filename(text: str, max_length: int = 50) -> str:
    """
    Convert text to safe filename format
    Examples:
        'Benzene' -> 'benzene'
        'Acetic Acid' -> 'acetic_acid'
        'CC(=O)O' -> 'smiles_cco'
    """
    # Replace spaces and special chars with underscores
    safe = ''.join(c if c.isalnum() else '_' for c in text.lower())
    # Remove consecutive underscores
    while '__' in safe:
        safe = safe.replace('__', '_')
    # Trim underscores from ends
    safe = safe.strip('_')
    # Limit length
    return safe[:max_length]


def create_cache_key(
    smiles_or_name: str,
    options: Optional[Dict] = None,
    source: str = "moleculeviewer"
) -> str:
    """
    Create a descriptive cache filename
    
    Examples:
        ('C1=CC=CC=C1', {}, 'mol2chemfig') -> 'benzene_default_m2cf'
        ('CC(=O)O', {'aromatic': True}, 'mol2chemfig') -> 'acetic_acid_aromatic_m2cf'
        ('caffeine', {'carbons': True, 'methyls': True}, 'moleculeviewer') -> 'caffeine_carbons_methyls_mv'
    
    Args:
        smiles_or_name: SMILES string or chemical name
        options: Dict of options (aromatic_circles, show_carbons, show_methyls, etc.)
        source: 'mol2chemfig' (m2cf) or 'moleculeviewer' (mv)
    
    Returns:
        Descriptive filename without extension
    """
    
    # Clean up the base name
    base_name = sanitize_filename(smiles_or_name)
    
    # Build option tags
    option_tags = []
    
    if options:
        # Standard moleculeviewer options
        if options.get('fancy_bonds'):
            option_tags.append('fancy')
        if options.get('aromatic_circles'):
            option_tags.append('aromatic')
        if options.get('show_carbons'):
            option_tags.append('carbons')
        if options.get('show_methyls'):
            option_tags.append('methyls')
        if options.get('atom_numbers'):
            option_tags.append('atoms')
        if options.get('recalculate_coordinates'):
            option_tags.append('recalc')
        
        hydrogens = options.get('hydrogens', 'keep')
        if hydrogens != 'keep':
            option_tags.append(f'h_{hydrogens}')
        
        flip = options.get('flip_horizontal') or options.get('flip_vertical')
        if flip:
            option_tags.append('flip')
        
        rotate = options.get('rotate', 0)
        if rotate != 0:
            option_tags.append(f'rot{rotate}')
        
        # Mol2chemfig specific options
        if options.get('angle') is not None:
            option_tags.append(f'angle{options.get("angle", 0)}')
        if options.get('indentation') is not None:
            option_tags.append(f'indent{options.get("indentation", 4)}')
        if options.get('h2') is not None and options.get('h2') != 'default':
            option_tags.append(f'h2_{options.get("h2")}')
    
    # If no options, mark as default
    if not option_tags:
        option_tags = ['default']
    
    # Add source tag
    source_tag = 'm2cf' if 'mol2chemfig' in source.lower() else 'mv'
    
    # Combine all parts
    parts = [base_name] + option_tags + [source_tag]
    
    # Create filename (no extension, will be added when saving)
    cache_key = '_'.join(filter(None, parts))
    
    return cache_key


def get_cache_hash(content: str) -> str:
    """Generate short hash for content (as backup for uniqueness)"""
    return hashlib.md5(content.encode()).hexdigest()[:8]


def save_svg_to_cache(
    svg_content: str,
    smiles_or_name: str,
    options: Optional[Dict] = None,
    source: str = "moleculeviewer",
    base_url: str = "http://localhost:5000"
) -> Tuple[str, str]:
    """
    Save SVG to cache with descriptive filename
    
    Returns:
        (cache_url, local_filepath)
        Example: ('/cache/caffeine_aromatic_m2cf.svg', '/path/to/caffeine_aromatic_m2cf.svg')
    """
    try:
        # Create cache key
        cache_key = create_cache_key(smiles_or_name, options, source)
        
        # Add content hash to ensure uniqueness
        content_hash = get_cache_hash(svg_content)
        filename = f"{cache_key}_{content_hash}.svg"
        filepath = os.path.join(CACHE_DIR, filename)
        
        # Write to cache
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(svg_content)
        
        print(f"✅ Cached SVG: {filename}")
        
        # Return URLs
        base = base_url.rstrip('/')
        cache_url = f"{base}/cache/{filename}"
        
        return cache_url, filepath
    
    except Exception as e:
        print(f"❌ Cache save error: {e}")
        return None, None


def get_cached_svg(
    smiles_or_name: str,
    options: Optional[Dict] = None,
    source: str = "moleculeviewer"
) -> Optional[str]:
    """
    Check if SVG already exists in cache
    
    Returns:
        Cache filename if found, None otherwise
    """
    try:
        cache_key = create_cache_key(smiles_or_name, options, source)
        
        # Look for files matching this key
        for filename in os.listdir(CACHE_DIR):
            if filename.startswith(cache_key + '_') and filename.endswith('.svg'):
                return filename
        
        return None
    
    except Exception as e:
        print(f"Cache lookup error: {e}")
        return None


def cleanup_old_cache():
    """Remove cache files older than CACHE_EXPIRY_HOURS"""
    try:
        now = time.time()
        removed_count = 0
        
        for filename in os.listdir(CACHE_DIR):
            filepath = os.path.join(CACHE_DIR, filename)
            if os.path.isfile(filepath):
                age_hours = (now - os.path.getmtime(filepath)) / 3600
                if age_hours > CACHE_EXPIRY_HOURS:
                    os.remove(filepath)
                    removed_count += 1
                    print(f"Cleaned old cache: {filename}")
        
        if removed_count > 0:
            print(f"🧹 Cache cleanup: Removed {removed_count} old files")
    
    except Exception as e:
        print(f"Cache cleanup error: {e}")


def get_cache_stats() -> Dict:
    """Get statistics about the cache"""
    try:
        files = os.listdir(CACHE_DIR)
        total_size = 0
        
        for filename in files:
            filepath = os.path.join(CACHE_DIR, filename)
            if os.path.isfile(filepath):
                total_size += os.path.getsize(filepath)
        
        return {
            'cache_dir': CACHE_DIR,
            'file_count': len(files),
            'total_size_mb': round(total_size / (1024 * 1024), 2),
            'expiry_hours': CACHE_EXPIRY_HOURS
        }
    
    except Exception as e:
        print(f"Cache stats error: {e}")
        return {}


# Test the cache naming system
if __name__ == '__main__':
    # Test cases
    test_cases = [
        ('C1=CC=CC=C1', {}, 'mol2chemfig'),
        ('C1=CC=CC=C1', {'aromatic_circles': True}, 'mol2chemfig'),
        ('CC(=O)O', {'show_carbons': True, 'show_methyls': True}, 'moleculeviewer'),
        ('caffeine', {'aromatic_circles': True, 'atom_numbers': True}, 'mol2chemfig'),
        ('aspirin', {'fancy_bonds': True, 'rotate': 90}, 'moleculeviewer'),
    ]
    
    print("Cache Key Examples:")
    print("-" * 60)
    for smiles, opts, src in test_cases:
        key = create_cache_key(smiles, opts, src)
        print(f"{smiles:20} + {src:15} = {key}")
