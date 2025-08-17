#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

def load_taxonomy_data(taxonomy_file):
    """Load taxonomy data using simple file reading for speed."""
    
    taxonomy_map = {}
    name_to_uid = {}
    
    with open(taxonomy_file, 'r') as f:
        # Skip header
        next(f)
        
        for line_num, line in enumerate(f, 2):
            line = line.strip()
            if not line:
                continue
                
            # Skip bacterial domain entries
            if 'in domain Bacteria' in line:
                continue
                
            parts = [p.strip() for p in line.split('|')]
            if len(parts) < 4:
                continue
                
            uid, parent_uid, name, rank = parts[:4]
            
            # Clean up empty parent_uid
            parent_uid = parent_uid if parent_uid else None
            rank = rank.lower() if rank else 'no rank'
            
            taxonomy_map[uid] = {
                'parent_uid': parent_uid,
                'name': name,
                'rank': rank
            }
            
            # Case-insensitive name lookup
            name_to_uid[name.lower()] = uid
            
            if line_num % 100000 == 0:
                print(f"Processed {line_num} lines...", file=sys.stderr)
    
    return taxonomy_map, name_to_uid

def get_taxonomic_lineage(species_name, taxonomy_map, name_to_uid):
    """Get taxonomic lineage for a species."""
    
    species_uid = name_to_uid.get(species_name.lower())
    if not species_uid:
        return None
    
    # Target ranks
    target_ranks = ['phylum', 'class', 'order', 'suborder', 'infraorder', 
                   'superfamily', 'family', 'genus', 'subgenus', 
                   'species group', 'species subgroup']
    
    taxonomy = {rank.capitalize(): 'NA' for rank in target_ranks}
    taxonomy['Species'] = species_name
    
    # Traverse lineage
    current_uid = species_uid
    visited = set()
    max_depth = 50  # Prevent infinite loops
    depth = 0
    
    while current_uid and current_uid in taxonomy_map and current_uid not in visited and depth < max_depth:
        visited.add(current_uid)
        info = taxonomy_map[current_uid]
        rank = info['rank']
        
        if rank in target_ranks:
            taxonomy[rank.capitalize()] = info['name']
        
        current_uid = info['parent_uid']
        depth += 1
    
    return taxonomy

def main():
    parser = argparse.ArgumentParser(description='Add taxonomic information to genome data')
    parser.add_argument('genomes_file', help='Input TSV file with genome data')
    parser.add_argument('taxonomy_file', help='Open Tree of Life taxonomy.tsv file')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.genomes_file).exists():
        sys.stderr.write(f"Error: Genomes file '{args.genomes_file}' not found\n")
        sys.exit(1)
    
    if not Path(args.taxonomy_file).exists():
        sys.stderr.write(f"Error: Taxonomy file '{args.taxonomy_file}' not found\n")
        sys.exit(1)
    
    if args.verbose:
        print(f"Loading taxonomy data from {args.taxonomy_file}...", file=sys.stderr)
    
    # Load taxonomy
    taxonomy_map, name_to_uid = load_taxonomy_data(args.taxonomy_file)
    
    if args.verbose:
        print(f"Loaded {len(taxonomy_map)} taxonomic entries", file=sys.stderr)
        print(f"Processing genomes from {args.genomes_file}...", file=sys.stderr)
    
    # Set up output
    output_file = open(args.output, 'w') if args.output else sys.stdout
    
    try:
        matched = 0
        total = 0
        
        with open(args.genomes_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                species_name = parts[0]
                total += 1
                
                taxonomy = get_taxonomic_lineage(species_name, taxonomy_map, name_to_uid)
                
                if taxonomy:
                    matched += 1
                    tax_info = [taxonomy.get(rank, 'NA') for rank in 
                               ['Phylum', 'Class', 'Order', 'Suborder', 'Infraorder',
                                'Superfamily', 'Family', 'Genus', 'Subgenus', 
                                'Species group', 'Species subgroup']]
                    output_line = line + '\t' + '\t'.join(tax_info)
                else:
                    # No match found
                    output_line = line + '\t' + '\t'.join(['NA'] * 11)
                
                print(output_line, file=output_file)
        
        if args.verbose:
            print(f"Matched {matched}/{total} species ({matched/total*100:.1f}%)", file=sys.stderr)
    
    finally:
        if args.output:
            output_file.close()

if __name__ == "__main__":
    main()