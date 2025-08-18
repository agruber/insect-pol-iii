#!/usr/bin/env python3
"""
cmsearch runner script - reads species names from stdin
Generates cmsearch commands for species received via pipe
"""

import sys
import argparse
import os

def generate_cmsearch_commands(species_names, models_file, output_dir="genomes", batch_size=None):
    """Generate cmsearch commands for species list"""
    commands = []
    batch_commands = []
    valid_species = []
    skipped_species = []
    
    # First pass: check which species have genomes available and haven't been processed
    for species in species_names:
        species = species.strip()
        if not species:
            continue
            
        species_dir = f"{output_dir}/{species}"
        genome_symlink = f"{species_dir}/genome.fna.gz"
        cmsearch_output = f"{species_dir}/cmsearch_results.txt"
        cmsearch_gz = f"{species_dir}/cmsearch_results.txt.gz"
        
        # Check if genome symlink exists and points to a valid file
        if os.path.islink(genome_symlink) and os.path.exists(genome_symlink):
            # Check if cmsearch results already exist
            if os.path.exists(cmsearch_output) or os.path.exists(cmsearch_gz):
                skipped_species.append(species)  # Already processed
            else:
                valid_species.append(species)
        else:
            skipped_species.append(species)  # No genome available
    
    # Skip species without genomes (no output)
    if not valid_species:
        return commands
    
    # Second pass: generate commands only for valid species
    for idx, species in enumerate(valid_species):
        species_dir = f"{output_dir}/{species}"
        genome_symlink = f"{species_dir}/genome.fna.gz"
        
        # Use the actual genome filename (without .gz) to preserve version info
        # Read the symlink target to get the actual filename
        if os.path.islink(genome_symlink):
            actual_filename = os.readlink(genome_symlink)
            # Remove .gz extension to get the uncompressed name
            temp_fasta = f"{species_dir}/{actual_filename.replace('.gz', '')}"
        else:
            temp_fasta = f"{species_dir}/genome.fna"
            
        cmsearch_output = f"{species_dir}/cmsearch_results.txt"
        species_commands = [
            f"echo 'Processing {species}...'",
            f"gunzip -c '{genome_symlink}' > '{temp_fasta}'",
            f"cmsearch --tblout '{cmsearch_output}' --noali '{models_file}' '{temp_fasta}'",
            f"gzip '{cmsearch_output}'",
            f"rm '{temp_fasta}'",
            f"echo 'Completed {species}'"
        ]
        
        batch_commands.extend(species_commands)
        
        # If batch_size is specified, create batches
        if batch_size and (idx + 1) % batch_size == 0:
            commands.extend(batch_commands)
            batch_commands = []
    
    # Add remaining commands if not using batches or there are leftover commands
    if not batch_size or batch_commands:
        commands.extend(batch_commands)
    
    return commands

def main():
    parser = argparse.ArgumentParser(description='Generate cmsearch commands from species list input')
    parser.add_argument('models_file', help='Path to the CM models file')
    parser.add_argument('-o', '--output-dir', default='genomes', help='Output directory for genomes (default: genomes)')
    parser.add_argument('-b', '--batch-size', type=int, help='Number of species per batch')
    
    args = parser.parse_args()
    
    # Read species names from stdin
    if sys.stdin.isatty():
        print("Error: This script expects species names from stdin. Use with filter_species.py", file=sys.stderr)
        print("Example: ./filter_species.py data/genomes_annotated.tsv 'Drosophila melanogaster' -f names | ./search_genomes.py data/models/all.cm", file=sys.stderr)
        sys.exit(1)
    
    # Read all input lines
    species_names = []
    for line in sys.stdin:
        line = line.strip()
        if line:
            species_names.append(line)
    
    if not species_names:
        print("Error: No species names received", file=sys.stderr)
        sys.exit(1)
    
    
    # Generate cmsearch commands
    commands = generate_cmsearch_commands(species_names, args.models_file, args.output_dir, args.batch_size)
    
    # Print commands
    for cmd in commands:
        print(cmd)

if __name__ == "__main__":
    main()