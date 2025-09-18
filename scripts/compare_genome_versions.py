#!/usr/bin/env python3
"""
Compare old and new versions of genomes_annotated.tsv to find changes and missing species.
Reports which species have changed or been removed, and checks if they exist in genomes/ folder.
"""

import sys
import pandas as pd
import argparse
from pathlib import Path

def clean_species_name(species_name):
    """Convert species name to filesystem-safe format (underscore format)"""
    return species_name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')

def load_genomes_file(filepath):
    """Load genomes_annotated.tsv file and return dataframe"""
    try:
        # Read TSV file - assuming species name is first column
        df = pd.read_csv(filepath, sep='\t', header=None)
        # Create a dict with species as key and entire row as value for comparison
        species_dict = {}
        for _, row in df.iterrows():
            species = row.iloc[0]  # First column is species name
            species_dict[species] = row.tolist()
        return species_dict
    except Exception as e:
        print(f"Error loading {filepath}: {e}", file=sys.stderr)
        sys.exit(1)

def check_genome_exists(species_name, genomes_dir="genomes"):
    """Check if genome folder exists for a species"""
    clean_name = clean_species_name(species_name)
    genome_path = Path(genomes_dir) / clean_name
    genome_file = genome_path / "genome.fna.gz"
    return genome_path.exists(), genome_file.exists(), clean_name

def compare_versions(old_file, new_file, genomes_dir="genomes"):
    """Compare two versions of genomes_annotated.tsv"""

    print(f"Loading old version: {old_file}")
    old_data = load_genomes_file(old_file)

    print(f"Loading new version: {new_file}")
    new_data = load_genomes_file(new_file)

    print(f"\nOld version has {len(old_data)} species")
    print(f"New version has {len(new_data)} species")
    print("=" * 80)

    # Track different types of changes
    removed_species = []
    changed_species = []
    unchanged_species = []
    new_species = []

    # Check each species in old version
    for species, old_row in old_data.items():
        if species not in new_data:
            # Species removed in new version
            removed_species.append(species)
        elif old_row != new_data[species]:
            # Species data changed
            changed_species.append(species)
        else:
            # Species unchanged
            unchanged_species.append(species)

    # Check for new species
    for species in new_data:
        if species not in old_data:
            new_species.append(species)

    # Report removed species
    if removed_species:
        print(f"\n🔴 REMOVED SPECIES ({len(removed_species)}):")
        print("-" * 40)
        for species in sorted(removed_species):
            folder_exists, genome_exists, clean_name = check_genome_exists(species, genomes_dir)
            if folder_exists:
                if genome_exists:
                    print(f"  {species}")
                    print(f"    └─ ⚠️  EXISTS in genomes/{clean_name}/ (with genome.fna.gz)")
                else:
                    print(f"  {species}")
                    print(f"    └─ ⚠️  EXISTS in genomes/{clean_name}/ (but NO genome.fna.gz)")
            else:
                print(f"  {species}")
                print(f"    └─ ✓ Not in genomes/{clean_name}/")

    # Report changed species
    if changed_species:
        print(f"\n🟡 CHANGED SPECIES ({len(changed_species)}):")
        print("-" * 40)
        for species in sorted(changed_species):
            folder_exists, genome_exists, clean_name = check_genome_exists(species, genomes_dir)
            print(f"  {species}")

            # Show what changed (comparing key fields)
            old_row = old_data[species]
            new_row = new_data[species]

            # Check if assembly/accession changed (typically column 3 for URL)
            if len(old_row) > 2 and len(new_row) > 2:
                old_url = str(old_row[2]) if len(old_row) > 2 else "N/A"
                new_url = str(new_row[2]) if len(new_row) > 2 else "N/A"

                if old_url != new_url:
                    old_accession = old_url.split('/')[-1] if '/' in old_url else old_url
                    new_accession = new_url.split('/')[-1] if '/' in new_url else new_url
                    print(f"    ├─ Assembly: {old_accession} → {new_accession}")

            if folder_exists:
                if genome_exists:
                    print(f"    └─ ⚠️  EXISTS in genomes/{clean_name}/ (needs update)")
                else:
                    print(f"    └─ ⚠️  EXISTS in genomes/{clean_name}/ (but NO genome.fna.gz)")
            else:
                print(f"    └─ ✓ Not in genomes/{clean_name}/")

    # Report new species
    if new_species:
        print(f"\n🟢 NEW SPECIES ({len(new_species)}):")
        print("-" * 40)
        for species in sorted(new_species):
            folder_exists, genome_exists, clean_name = check_genome_exists(species, genomes_dir)
            if folder_exists:
                if genome_exists:
                    print(f"  {species}")
                    print(f"    └─ ✓ Already in genomes/{clean_name}/")
                else:
                    print(f"  {species}")
                    print(f"    └─ ⚠️  Folder exists but NO genome.fna.gz")
            else:
                print(f"  {species}")
                print(f"    └─ To download: genomes/{clean_name}/")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY:")
    print(f"  Unchanged: {len(unchanged_species)} species")
    print(f"  Changed:   {len(changed_species)} species")
    print(f"  Removed:   {len(removed_species)} species")
    print(f"  New:       {len(new_species)} species")

    # Check which changed/removed species we have locally
    local_affected = []
    for species in changed_species + removed_species:
        folder_exists, genome_exists, clean_name = check_genome_exists(species, genomes_dir)
        if genome_exists:
            local_affected.append((species, clean_name))

    if local_affected:
        print(f"\n⚠️  ACTION REQUIRED: {len(local_affected)} species in genomes/ need attention:")
        for species, clean_name in local_affected:
            status = "REMOVED" if species in removed_species else "CHANGED"
            print(f"  - genomes/{clean_name}/ ({status}: {species})")

def main():
    parser = argparse.ArgumentParser(
        description='Compare old and new versions of genomes_annotated.tsv',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two versions
  python3 compare_genome_versions.py data/genomes_annotated.old.tsv data/genomes_annotated.tsv

  # Use custom genomes directory
  python3 compare_genome_versions.py old.tsv new.tsv --genomes-dir /path/to/genomes
        """
    )

    parser.add_argument('old_file', help='Path to old genomes_annotated.tsv')
    parser.add_argument('new_file', help='Path to new genomes_annotated.tsv')
    parser.add_argument('--genomes-dir', default='genomes',
                       help='Path to genomes directory (default: genomes)')

    args = parser.parse_args()

    # Check if files exist
    if not Path(args.old_file).exists():
        print(f"Error: Old file not found: {args.old_file}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.new_file).exists():
        print(f"Error: New file not found: {args.new_file}", file=sys.stderr)
        sys.exit(1)

    # Run comparison
    compare_versions(args.old_file, args.new_file, args.genomes_dir)

if __name__ == "__main__":
    main()