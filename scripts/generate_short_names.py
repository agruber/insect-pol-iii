#!/usr/bin/env python3
"""
Generate unique short names for species from genomes_annotated.tsv
Format: First 3 letters of genus + underscore + first 3+ letters of species
Extends species name length to ensure uniqueness.
"""

import sys
import pandas as pd
from collections import defaultdict

def clean_species_name(name):
    """Clean species name by removing problematic characters."""
    return name.replace(' ', '_').replace('(', '').replace(')', '').replace('[', '').replace(']', '')

def generate_short_name(genus, species, min_species_len=3):
    """Generate short name: first 3 letters of genus + _ + species letters."""
    genus_short = genus[:3].capitalize()
    species_short = species[:min_species_len].lower()
    return f"{genus_short}_{species_short}"

def make_unique_short_names(species_list):
    """Generate unique short names - first come, first served with shortest names."""
    name_mapping = {}
    used_names = set()
    
    for species in species_list:
        parts = species.strip().split()
        if len(parts) < 2:
            continue
            
        genus = parts[0]
        species_epithet = parts[1]
        genus_short = genus[:3].capitalize()
        
        # Try progressively longer species names until we find one that's not used
        for species_len in range(3, len(species_epithet) + 1):
            species_short = species_epithet[:species_len].lower()
            candidate_name = f"{genus_short}_{species_short}"
            
            if candidate_name not in used_names:
                name_mapping[species] = candidate_name
                used_names.add(candidate_name)
                break
        else:
            # Fallback: use full name with number suffix if needed
            base_name = f"{genus_short}_{species_epithet.lower()}"
            if base_name not in used_names:
                name_mapping[species] = base_name
                used_names.add(base_name)
            else:
                counter = 1
                while f"{base_name}_{counter}" in used_names:
                    counter += 1
                final_name = f"{base_name}_{counter}"
                name_mapping[species] = final_name
                used_names.add(final_name)
    
    return name_mapping

def main():
    if len(sys.argv) != 2:
        print("Usage: generate_short_names.py <genomes_annotated.tsv>", file=sys.stderr)
        sys.exit(1)
    
    tsv_file = sys.argv[1]
    
    try:
        # Read TSV file
        df = pd.read_csv(tsv_file, sep='\t', header=None)
        species_list = df.iloc[:, 0].unique()  # Get unique species from column 1
        
        # Generate unique short names
        name_mapping = make_unique_short_names(species_list)
        
        # Output mapping
        print("Full_Species_Name\tShort_Name")
        for full_name in sorted(name_mapping.keys()):
            short_name = name_mapping[full_name]
            print(f"{full_name}\t{short_name}")
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()