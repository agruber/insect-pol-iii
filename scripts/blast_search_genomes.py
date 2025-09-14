#!/usr/bin/env python3
"""
BLAST search runner script - reads species names from stdin
Generates blastn commands for species received via pipe using a query FASTA file
"""

import sys
import argparse
import os

def generate_blast_commands(species_names, query_fasta, output_dir="genomes", batch_size=None, 
                          evalue=1e-5, num_threads=1, blast_type="blastn"):
    """Generate BLAST commands for species list"""
    commands = []
    batch_commands = []
    valid_species = []
    skipped_species = []
    
    # Extract query name from fasta file for output naming
    query_basename = os.path.basename(query_fasta)
    query_name = os.path.splitext(query_basename)[0]  # Remove .fa/.fasta extension
    
    # First pass: check which species have genomes available and haven't been processed
    for species in species_names:
        species = species.strip()
        if not species:
            continue
            
        species_dir = f"{output_dir}/{species}"
        genome_symlink = f"{species_dir}/genome.fna.gz"
        blast_output = f"{species_dir}/{query_name}.blast"
        blast_gz = f"{species_dir}/{query_name}.blast.gz"
        
        # Check if genome symlink exists and points to a valid file
        if os.path.islink(genome_symlink) and os.path.exists(genome_symlink):
            # Check if BLAST results already exist
            if os.path.exists(blast_output) or os.path.exists(blast_gz):
                skipped_species.append(species)  # Already processed
            else:
                valid_species.append(species)
        else:
            skipped_species.append(species)  # No genome available
    
    # Skip species without genomes (no output)
    if not valid_species:
        return commands
    
    # Report progress to stderr
    print(f"# Found {len(valid_species)} species to process", file=sys.stderr)
    if skipped_species:
        print(f"# Skipping {len(skipped_species)} species (no genome or already processed)", file=sys.stderr)
    
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
            
        blast_output = f"{species_dir}/{query_name}.blast"
        
        # Build BLAST command with appropriate output format
        # Format 6 is tabular output
        blast_cmd = (f"{blast_type} -query '{query_fasta}' -subject '{temp_fasta}' "
                    f"-outfmt '6 qseqid sseqid pident length mismatch gapopen "
                    f"qstart qend sstart send evalue bitscore qlen slen sstrand' "
                    f"-evalue {evalue} -num_threads {num_threads}")
        
        species_commands = [
            f"echo 'Processing {species} with BLAST...'",
            f"gunzip -c '{genome_symlink}' > '{temp_fasta}'",
            f"{blast_cmd} | ./scripts/process_blast_hits.py > '{blast_output}'",
            f"gzip '{blast_output}'",
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

def validate_blast_tool(blast_type):
    """Check if BLAST tool is available"""
    import shutil
    if not shutil.which(blast_type):
        print(f"Error: {blast_type} not found in PATH", file=sys.stderr)
        print(f"Please install BLAST+ tools: sudo apt install ncbi-blast+", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Generate BLAST search commands from species list input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script generates BLAST commands to search for sequences in genome files.
It reads species names from stdin and creates commands to:
1. Decompress genome files
2. Run BLAST search with the query sequences
3. Compress results and clean up

Examples:
  # Basic usage with a query FASTA file
  echo "Drosophila_melanogaster" | ./blast_search_genomes.py data/seqs/query.fa
  
  # With species filter
  ./filter_species.py data/genomes.tsv 'Insecta' | ./blast_search_genomes.py data/seqs/noeCR34335.fa
  
  # With custom e-value and threads
  ./filter_species.py data/genomes.tsv 'Diptera' | ./blast_search_genomes.py data/seqs/query.fa -e 1e-10 -t 4
  
  # Using tblastn for protein queries
  echo "Apis_mellifera" | ./blast_search_genomes.py protein_query.fa --blast-type tblastn

Output format:
  The script generates .blast.gz files with tabular BLAST output (format 6):
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand
        """
    )
    
    parser.add_argument('query_fasta', help='Path to the query FASTA file')
    parser.add_argument('-o', '--output-dir', default='genomes', 
                       help='Output directory for genomes (default: genomes)')
    parser.add_argument('-b', '--batch-size', type=int, 
                       help='Number of species per batch')
    parser.add_argument('-e', '--evalue', type=float, default=1e-5,
                       help='E-value threshold for BLAST (default: 1e-5)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Number of threads for BLAST (default: 1)')
    parser.add_argument('--blast-type', choices=['blastn', 'tblastn', 'blastx', 'tblastx'], 
                       default='blastn',
                       help='Type of BLAST to use (default: blastn)')
    
    args = parser.parse_args()
    
    # Check if query file exists
    if not os.path.exists(args.query_fasta):
        print(f"Error: Query file '{args.query_fasta}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Check if BLAST tool is available
    validate_blast_tool(args.blast_type)
    
    # Read species names from stdin
    if sys.stdin.isatty():
        print("Error: This script expects species names from stdin. Use with filter_species.py", file=sys.stderr)
        print(f"Example: ./filter_species.py data/genomes.tsv 'Insecta' | ./blast_search_genomes.py {args.query_fasta}", file=sys.stderr)
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
    
    print(f"# Received {len(species_names)} species names", file=sys.stderr)
    print(f"# Using query file: {args.query_fasta}", file=sys.stderr)
    print(f"# BLAST type: {args.blast_type}, E-value: {args.evalue}, Threads: {args.threads}", file=sys.stderr)
    
    # Generate BLAST commands
    commands = generate_blast_commands(
        species_names, 
        args.query_fasta, 
        args.output_dir, 
        args.batch_size,
        args.evalue,
        args.threads,
        args.blast_type
    )
    
    # Print commands
    for cmd in commands:
        print(cmd)

if __name__ == "__main__":
    main()