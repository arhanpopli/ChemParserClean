#!/usr/bin/env python3
"""
Cache Deduplication Utility

This script analyzes and deduplicates cache entries across both servers:
- MoleculeViewer (cache/moleculeviewer/)
- Mol2ChemFig (cache/mol2chemfig/)

The issue: When searching by name (e.g., "ethanol"), it caches with that name,
but when searching by SMILES ("CCO"), it creates a new cache entry even though
they're the same molecule.

Solution:
1. Parse all cache filenames and extract metadata
2. For each SVG, attempt to extract SMILES from associated metadata
3. Canonicalize SMILES and identify duplicates
4. Merge duplicate entries, keeping the most recent one
5. Generate a report of space saved

Usage:
    python deduplicate_cache.py --dry-run      # Preview what would be deleted
    python deduplicate_cache.py --execute      # Actually delete duplicates
    python deduplicate_cache.py --report       # Generate detailed report
"""

import os
import json
import hashlib
import argparse
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import re

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. SMILES canonicalization will be limited.")

def canonicalize_smiles(smiles):
    """Canonicalize SMILES string to ensure consistent representation"""
    if not RDKIT_AVAILABLE:
        return smiles.strip()

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles.strip()
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return smiles.strip()

def extract_smiles_from_svg(svg_path):
    """
    Attempt to extract SMILES from SVG metadata or filename pattern
    Returns None if SMILES cannot be determined
    """
    # Try to read SVG content for embedded metadata
    try:
        with open(svg_path, 'r', encoding='utf-8') as f:
            content = f.read(1000)  # Read first 1KB
            # Look for SMILES in metadata comments
            match = re.search(r'SMILES:\s*([^\s<>"]+)', content)
            if match:
                return match.group(1)
    except:
        pass

    return None

def analyze_cache_directory(cache_dir):
    """
    Analyze a cache directory and group files by canonical SMILES

    Returns:
        dict: {canonical_smiles: [list of file paths with metadata]}
    """
    cache_path = Path(cache_dir)
    if not cache_path.exists():
        print(f"Cache directory not found: {cache_dir}")
        return {}, []

    files_by_smiles = defaultdict(list)
    unknown_files = []

    # Get all SVG files
    svg_files = list(cache_path.glob('*.svg'))

    print(f"\nAnalyzing {len(svg_files)} files in {cache_dir}...")

    for svg_file in svg_files:
        file_info = {
            'path': svg_file,
            'size': svg_file.stat().st_size,
            'modified': datetime.fromtimestamp(svg_file.stat().st_mtime),
            'smiles': None,
            'canonical_smiles': None
        }

        # Extract SMILES from SVG
        smiles = extract_smiles_from_svg(svg_file)

        if smiles:
            file_info['smiles'] = smiles
            file_info['canonical_smiles'] = canonicalize_smiles(smiles)
            files_by_smiles[file_info['canonical_smiles']].append(file_info)
        else:
            unknown_files.append(file_info)

    print(f"  - Files with identifiable SMILES: {sum(len(v) for v in files_by_smiles.values())}")
    print(f"  - Files without SMILES metadata: {len(unknown_files)}")

    return files_by_smiles, unknown_files

def find_duplicates(files_by_smiles):
    """
    Find duplicate cache entries for the same molecule

    Returns:
        dict: {canonical_smiles: [list of duplicate files]} (only entries with duplicates)
    """
    duplicates = {}

    for smiles, files in files_by_smiles.items():
        if len(files) > 1:
            duplicates[smiles] = files

    return duplicates

def generate_report(duplicates, unknown_files):
    """Generate a detailed report of duplicates found"""
    print("\n" + "="*80)
    print("CACHE DEDUPLICATION REPORT")
    print("="*80)

    total_duplicates = sum(len(files) - 1 for files in duplicates.values())
    total_wasted_space = sum(
        sum(f['size'] for f in files[1:])  # Sum size of all but the first file
        for files in duplicates.values()
    )

    print(f"\nSummary:")
    print(f"  - Unique molecules with duplicates: {len(duplicates)}")
    print(f"  - Total duplicate files: {total_duplicates}")
    print(f"  - Wasted disk space: {total_wasted_space / 1024:.2f} KB ({total_wasted_space / (1024*1024):.2f} MB)")

    if duplicates:
        print(f"\nDuplicate entries by molecule:")
        print("-" * 80)

        for smiles, files in sorted(duplicates.items(), key=lambda x: len(x[1]), reverse=True):
            print(f"\nCanonical SMILES: {smiles}")
            print(f"Duplicate count: {len(files)} files")
            print(f"Total size: {sum(f['size'] for f in files) / 1024:.2f} KB")
            print(f"Wasted space: {sum(f['size'] for f in files[1:]) / 1024:.2f} KB")

            # Sort files by modification time (newest first)
            files_sorted = sorted(files, key=lambda x: x['modified'], reverse=True)

            print(f"Files (keeping newest):")
            for i, file_info in enumerate(files_sorted):
                status = "KEEP" if i == 0 else "DELETE"
                print(f"  [{status}] {file_info['path'].name}")
                print(f"        Size: {file_info['size']} bytes, Modified: {file_info['modified']}")

    if unknown_files:
        print(f"\n\nFiles without SMILES metadata ({len(unknown_files)} files):")
        print("-" * 80)
        print("These files cannot be automatically deduplicated.")
        for file_info in unknown_files[:10]:  # Show first 10
            print(f"  - {file_info['path'].name} ({file_info['size']} bytes)")
        if len(unknown_files) > 10:
            print(f"  ... and {len(unknown_files) - 10} more")

    print("\n" + "="*80)

def deduplicate_cache(duplicates, dry_run=True):
    """
    Remove duplicate cache entries, keeping the most recent version

    Args:
        duplicates: Dictionary of duplicate files grouped by canonical SMILES
        dry_run: If True, only simulate deletion without actually removing files
    """
    total_deleted = 0
    total_space_freed = 0

    print("\n" + "="*80)
    if dry_run:
        print("DRY RUN MODE - No files will be deleted")
    else:
        print("EXECUTING DEDUPLICATION - Files will be deleted")
    print("="*80)

    for smiles, files in duplicates.items():
        # Sort by modification time (newest first)
        files_sorted = sorted(files, key=lambda x: x['modified'], reverse=True)

        # Keep the newest, delete the rest
        keep_file = files_sorted[0]
        delete_files = files_sorted[1:]

        print(f"\nMolecule: {smiles}")
        print(f"  Keeping: {keep_file['path'].name} (modified {keep_file['modified']})")

        for file_info in delete_files:
            if dry_run:
                print(f"  [DRY RUN] Would delete: {file_info['path'].name} ({file_info['size']} bytes)")
            else:
                try:
                    file_info['path'].unlink()
                    print(f"  [DELETED] {file_info['path'].name} ({file_info['size']} bytes)")
                    total_deleted += 1
                    total_space_freed += file_info['size']
                except Exception as e:
                    print(f"  [ERROR] Failed to delete {file_info['path'].name}: {e}")

    print("\n" + "="*80)
    if dry_run:
        print(f"DRY RUN COMPLETE")
        print(f"Would delete: {sum(len(files) - 1 for files in duplicates.values())} files")
        print(f"Would free: {sum(sum(f['size'] for f in files[1:]) for files in duplicates.values()) / 1024:.2f} KB")
    else:
        print(f"DEDUPLICATION COMPLETE")
        print(f"Deleted: {total_deleted} files")
        print(f"Space freed: {total_space_freed / 1024:.2f} KB ({total_space_freed / (1024*1024):.2f} MB)")
    print("="*80 + "\n")

def main():
    parser = argparse.ArgumentParser(
        description='Deduplicate cache entries across ChemParser servers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python deduplicate_cache.py --dry-run      # Preview what would be deleted
  python deduplicate_cache.py --execute      # Actually delete duplicates
  python deduplicate_cache.py --report       # Generate detailed report only
  python deduplicate_cache.py --cache-dir path/to/cache --dry-run
        """
    )

    parser.add_argument('--dry-run', action='store_true',
                        help='Simulate deduplication without deleting files')
    parser.add_argument('--execute', action='store_true',
                        help='Execute deduplication and delete duplicate files')
    parser.add_argument('--report', action='store_true',
                        help='Generate detailed report only (no deletion)')
    parser.add_argument('--cache-dir', type=str,
                        help='Specific cache directory to analyze (default: analyze all)')

    args = parser.parse_args()

    # Default: dry-run if no mode specified
    if not (args.dry_run or args.execute or args.report):
        args.dry_run = True

    # Determine which cache directories to process
    cache_dirs = []
    if args.cache_dir:
        cache_dirs.append(args.cache_dir)
    else:
        # Default: process both cache directories
        cache_dirs.extend([
            'MoleculeViewer/cache/moleculeviewer',
            'cache/mol2chemfig'
        ])

    print("\n" + "="*80)
    print("CHEMPARSER CACHE DEDUPLICATION UTILITY")
    print("="*80)
    print(f"Mode: {'DRY RUN' if args.dry_run else 'EXECUTE' if args.execute else 'REPORT ONLY'}")
    print(f"RDKit available: {RDKIT_AVAILABLE}")
    print("="*80)

    # Process each cache directory
    all_duplicates = {}
    all_unknown_files = []

    for cache_dir in cache_dirs:
        files_by_smiles, unknown_files = analyze_cache_directory(cache_dir)
        duplicates = find_duplicates(files_by_smiles)

        if duplicates or unknown_files:
            all_duplicates.update(duplicates)
            all_unknown_files.extend(unknown_files)

    # Generate report
    if all_duplicates or all_unknown_files:
        generate_report(all_duplicates, all_unknown_files)
    else:
        print("\nNo duplicates found! Cache is already optimized.")
        return

    # Execute deduplication if requested
    if all_duplicates and (args.execute or args.dry_run) and not args.report:
        deduplicate_cache(all_duplicates, dry_run=args.dry_run)

if __name__ == '__main__':
    main()
