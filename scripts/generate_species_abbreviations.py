#!/usr/bin/env python3
"""
Generate species abbreviations from genomes_annotated.tsv
Creates a TSV file with species name, underscore name, and unique abbreviation
"""

import sys
import os

def main():
    unique_abbrev = {}

    # Read the genomes_annotated.tsv file
    with open("data/genomes_annotated.tsv", 'r') as f:
        for line in f:
            line = line.strip()
            columns = line.split('\t')

            phylum = columns[3]

            # Only process Arthropoda
            if phylum == 'Arthropoda':
                species = columns[0]
                species_underscore = species.replace(' ', '_')

                # Generate abbreviation like the Perl script
                words = species.split()
                if len(words) >= 2:
                    jj = 3  # characters from first word
                    kk = 3  # characters from second word

                    abbrev2 = words[0][:jj] + "_" + words[1][:kk]

                    # Make abbreviation unique
                    while abbrev2 in unique_abbrev:
                        kk += 1
                        if kk > len(words[1]):
                            kk = len(words[1])
                            jj += 1
                        abbrev2 = words[0][:jj] + "_" + words[1][:kk]

                    unique_abbrev[abbrev2] = True

                    # Output TSV line
                    print(f"{species}\t{species_underscore}\t{abbrev2}")

if __name__ == "__main__":
    main()