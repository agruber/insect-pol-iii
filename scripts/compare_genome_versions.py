#!/usr/bin/env python3
"""
Check genomes/ directories against genomes_annotated.tsv to find version mismatches and obsolete directories.
Reports which species have genome version mismatches and which directories are no longer in the database.
"""

import sys
import pandas as pd
import argparse
from pathlib import Path
import gzip

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

def get_genome_version(genome_file):
    """Extract genome version/accession from genome.fna.gz file"""
    try:
        with gzip.open(genome_file, 'rt') as f:
            # Read first line (header)
            first_line = f.readline().strip()
            if first_line.startswith('>'):
                # Try to extract accession from header
                # Typically format is >accession description
                parts = first_line[1:].split()
                if parts:
                    return parts[0]  # First part is usually the accession
        return None
    except Exception as e:
        return None

def get_tblout_genome_version(tblout_file):
    """Extract genome version from tblout.gz file's 'Target file:' line"""
    try:
        with gzip.open(tblout_file, 'rt') as f:
            for line in f:
                if line.startswith('# Target file:'):
                    # Format: # Target file:     genomes/Species_name/GCA_XXXXXX_genomic.fna
                    parts = line.split('/')
                    if len(parts) >= 3:
                        filename = parts[-1].strip()
                        # Extract just the assembly part (GCA_XXXXXX)
                        # Remove _genomic.fna or .fasta_genomic.fna suffix
                        if '_genomic.fna' in filename:
                            assembly = filename.split('_genomic.fna')[0]
                            # Also handle .fasta in the name
                            if '.fasta' in assembly:
                                assembly = assembly.replace('.fasta', '')
                            return assembly
                    break
        return None
    except Exception as e:
        return None

def check_versions(new_file, genomes_dir="genomes"):
    """Check genome directories against the reference file for mismatches and obsolete directories"""

    print(f"Loading reference file: {new_file}")
    reference_data = load_genomes_file(new_file)
    print(f"Reference file has {len(reference_data)} species")
    print("=" * 80)

    # Create a mapping of clean names to species for reverse lookup
    clean_name_to_species = {}
    for species in reference_data:
        clean_name = clean_species_name(species)
        clean_name_to_species[clean_name] = species

    # Get all species directories in genomes/
    genomes_path = Path(genomes_dir)
    if not genomes_path.exists():
        print(f"Error: Genomes directory not found: {genomes_dir}", file=sys.stderr)
        sys.exit(1)

    all_dirs = [d for d in genomes_path.iterdir() if d.is_dir()]
    print(f"Found {len(all_dirs)} directories in {genomes_dir}/\n")

    # Track issues
    version_mismatches = []
    tblout_mismatches = []
    obsolete_directories = []
    missing_genome = []

    # Check each directory in genomes/
    for species_dir in sorted(all_dirs):
        dir_name = species_dir.name
        genome_file = species_dir / "genome.fna.gz"

        # Check if this directory corresponds to a species in reference
        if dir_name not in clean_name_to_species:
            # Directory doesn't match any species in reference
            obsolete_directories.append(dir_name)
            continue

        species_name = clean_name_to_species[dir_name]

        # Check if genome file exists
        if not genome_file.exists():
            missing_genome.append((species_name, dir_name))
            continue

        # Get expected URL/assembly from reference
        species_row = reference_data[species_name]
        if len(species_row) > 2:
            expected_url = str(species_row[2])
            expected_assembly = expected_url.split('/')[-1] if '/' in expected_url else expected_url
        else:
            expected_assembly = None
            continue

        # Look for the original assembly file (e.g., GCA_*_genomic.fna.gz)
        actual_assembly = None
        for genomic_file in species_dir.glob("GC*_genomic.fna.gz"):
            # Extract assembly name from filename
            filename = genomic_file.name
            # Remove _genomic.fna.gz suffix
            actual_assembly = filename.replace('_genomic.fna.gz', '')
            break

        # Check for version mismatch
        if expected_assembly and actual_assembly:
            if expected_assembly != actual_assembly:
                version_mismatches.append((species_name, dir_name, expected_assembly, actual_assembly))
        elif expected_assembly and not actual_assembly:
            # Could not find assembly file - might be an issue
            version_mismatches.append((species_name, dir_name, expected_assembly, "Unknown (no GC*_genomic.fna.gz file found)"))

        # Check tblout.gz files for genome version mismatches
        tblout_files = list(species_dir.glob("*.tblout.gz"))
        for tblout_file in tblout_files:
            tblout_assembly = get_tblout_genome_version(tblout_file)
            if tblout_assembly and expected_assembly:
                if expected_assembly != tblout_assembly:
                    tblout_mismatches.append((species_name, dir_name, tblout_file.name, expected_assembly, tblout_assembly))

    # Report version mismatches
    if version_mismatches:
        print(f"\n🟡 GENOME VERSION MISMATCHES ({len(version_mismatches)}):")
        print("-" * 80)
        for species, dir_name, expected, actual in sorted(version_mismatches):
            print(f"  {species} (genomes/{dir_name}/)")
            print(f"    Expected: {expected}")
            print(f"    Actual:   {actual}")
            print()

    # Report tblout file mismatches
    if tblout_mismatches:
        print(f"\n🟠 TBLOUT FILE VERSION MISMATCHES ({len(tblout_mismatches)}):")
        print("-" * 80)
        for species, dir_name, tblout_name, expected, actual in sorted(tblout_mismatches):
            print(f"  {species} (genomes/{dir_name}/)")
            print(f"    File:     {tblout_name}")
            print(f"    Expected: {expected}")
            print(f"    Actual:   {actual}")
            print()

    # Report obsolete directories
    if obsolete_directories:
        print(f"\n🔴 OBSOLETE DIRECTORIES ({len(obsolete_directories)}):")
        print("-" * 80)
        print("These directories exist in genomes/ but are not in the reference file:")
        for dir_name in sorted(obsolete_directories):
            print(f"  genomes/{dir_name}/")
        print()

    # Report missing genome files
    if missing_genome:
        print(f"\n⚠️  MISSING GENOME FILES ({len(missing_genome)}):")
        print("-" * 80)
        for species, dir_name in sorted(missing_genome):
            print(f"  {species} (genomes/{dir_name}/) - no genome.fna.gz")
        print()

    # Summary
    print("=" * 80)
    print("SUMMARY:")
    print(f"  Total directories checked: {len(all_dirs)}")
    print(f"  Genome version mismatches: {len(version_mismatches)}")
    print(f"  Tblout version mismatches: {len(tblout_mismatches)}")
    print(f"  Obsolete directories:      {len(obsolete_directories)}")
    print(f"  Missing genome files:      {len(missing_genome)}")
    print(f"  OK:                        {len(all_dirs) - len(version_mismatches) - len(obsolete_directories) - len(missing_genome)}")

    total_issues = len(version_mismatches) + len(tblout_mismatches) + len(obsolete_directories)
    if total_issues > 0:
        print(f"\n⚠️  ACTION REQUIRED: {total_issues} issues found")
        if version_mismatches:
            print(f"    - {len(version_mismatches)} genome(s) need updating")
        if tblout_mismatches:
            print(f"    - {len(tblout_mismatches)} tblout file(s) need re-running")
        if obsolete_directories:
            print(f"    - {len(obsolete_directories)} obsolete director(ies) to remove")

def main():
    parser = argparse.ArgumentParser(
        description='Check genome directories for version mismatches and obsolete directories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check genomes against reference file
  python3 compare_genome_versions.py genomes_annotated.tsv

  # Use custom genomes directory
  python3 compare_genome_versions.py genomes_annotated.tsv --genomes-dir /path/to/genomes
        """
    )

    parser.add_argument('reference_file', help='Path to reference genomes_annotated.tsv file')
    parser.add_argument('--genomes-dir', default='genomes',
                       help='Path to genomes directory (default: genomes)')

    args = parser.parse_args()

    # Check if file exists
    if not Path(args.reference_file).exists():
        print(f"Error: Reference file not found: {args.reference_file}", file=sys.stderr)
        sys.exit(1)

    # Run check
    check_versions(args.reference_file, args.genomes_dir)

if __name__ == "__main__":
    main()