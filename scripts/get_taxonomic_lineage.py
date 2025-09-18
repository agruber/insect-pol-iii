#!/usr/bin/env python3
"""
Script to extract taxonomic lineage from genomes_annotated.tsv
"""

import sys
import argparse

def get_taxonomic_lineage(species_name, tsv_file):
    """
    Extract taxonomic lineage for a given species from the TSV file.
    Returns a formatted lineage string.
    """
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 14 and fields[0] == species_name:
                    # Extract taxonomic levels (skipping NA values)
                    taxonomic_levels = []

                    # Define the taxonomic hierarchy (0-indexed)
                    level_names = [
                        fields[6],  # Suborder (Brachycera) - index 6
                        fields[7],  # Infraorder (Muscomorpha) - index 7
                        fields[8],  # Superfamily (Ephydroidea) - index 8
                        fields[9],  # Family (Drosophilidae) - index 9
                        fields[10], # Genus (Drosophila) - index 10
                        fields[11], # Subgenus (Sophophora) - index 11
                        fields[12], # Species group (melanogaster group) - index 12
                        fields[13]  # Subgroup (melanogaster subgroup) - index 13
                    ]

                    # Add non-NA levels to the lineage
                    for level in level_names:
                        if level and level != 'NA':
                            taxonomic_levels.append(level)

                    # Add the species name at the end
                    taxonomic_levels.append(species_name)

                    # Join with ' > ' separator
                    return ' > '.join(taxonomic_levels)

        return f"Taxonomic lineage not found for {species_name}"

    except FileNotFoundError:
        return f"Error: Could not find file {tsv_file}"
    except Exception as e:
        return f"Error: {str(e)}"

def main():
    parser = argparse.ArgumentParser(description='Get taxonomic lineage for a species')
    parser.add_argument('species_name', help='Species name to look up')
    parser.add_argument('--tsv_file', default='data/genomes_annotated.tsv',
                       help='Path to genomes_annotated.tsv file')

    args = parser.parse_args()

    lineage = get_taxonomic_lineage(args.species_name, args.tsv_file)
    print(lineage)

if __name__ == '__main__':
    main()