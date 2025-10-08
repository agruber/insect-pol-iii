#!/usr/bin/env python3
"""
Compare previous noeCR34335 results from results_summary.csv with current results/ folder data
"""

import sys
import csv
from pathlib import Path
import argparse
from collections import defaultdict

def load_taxonomy(tsv_file="data/genomes_annotated.tsv"):
    """Load taxonomy information from genomes_annotated.tsv"""
    taxonomy = {}

    tsv_path = Path(tsv_file)
    if not tsv_path.exists():
        print(f"Warning: Taxonomy file not found: {tsv_file}", file=sys.stderr)
        return taxonomy

    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            species = parts[0]
            phylum = parts[3] if parts[3] != 'NA' else ''
            class_ = parts[4] if parts[4] != 'NA' else ''
            order = parts[5] if parts[5] != 'NA' else ''
            suborder = parts[6] if parts[6] != 'NA' else ''
            infraorder = parts[7] if parts[7] != 'NA' else ''
            superfamily = parts[8] if parts[8] != 'NA' else ''
            family = parts[9] if parts[9] != 'NA' else ''

            # Build lineage string
            lineage_parts = []
            for taxon in [phylum, class_, order, suborder, infraorder, superfamily, family]:
                if taxon:
                    lineage_parts.append(taxon)

            taxonomy[species] = ' > '.join(lineage_parts)

    return taxonomy

def clean_species_name(species_name):
    """Convert species name to filesystem-safe format"""
    return species_name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')

def load_previous_results(csv_file):
    """Load previous results from results_summary.csv"""
    species_hits = defaultdict(list)

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            species = row['Species']
            sequence_name = row['Sequence_Name']

            # Store details for each hit
            hit_info = {
                'name': sequence_name,
                'length': int(row['Length']) if row['Length'] else 0,
                'noe_evalue': row['noe E-value'] if row['noe E-value'] else 'N/A',
                'contig': row['Contig'],
                'start': row['Start'],
                'end': row['End'],
                'strand': row['Strand']
            }
            species_hits[species].append(hit_info)

    return species_hits

def load_current_results(results_dir="results"):
    """Load current results from results/ folder"""
    results_path = Path(results_dir)
    species_hits = defaultdict(list)

    if not results_path.exists():
        print(f"Warning: Results directory {results_dir} does not exist", file=sys.stderr)
        return species_hits

    # For each species directory in results/
    for species_dir in results_path.iterdir():
        if not species_dir.is_dir():
            continue

        # Convert directory name to species name (underscore to space)
        species_name = species_dir.name.replace('_', ' ')

        # Check for lncrna.fa file
        lncrna_file = species_dir / "lncrna.fa"
        if not lncrna_file.exists():
            continue

        # Parse lncrna.fa to get current hits
        with open(lncrna_file, 'r') as f:
            current_seq = None
            current_header = None

            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if any
                    if current_seq and current_header:
                        # Parse header format: >ID|noeCR34335|contig|start|end|strand
                        parts = current_header[1:].split('|')
                        if len(parts) >= 6:
                            hit_info = {
                                'name': parts[0],
                                'length': len(current_seq),
                                'contig': parts[2],
                                'start': parts[3],
                                'end': parts[4],
                                'strand': parts[5]
                            }
                            species_hits[species_name].append(hit_info)

                    # Start new sequence
                    current_header = line
                    current_seq = ""
                else:
                    if current_seq is not None:
                        current_seq += line

            # Don't forget the last sequence
            if current_seq and current_header:
                parts = current_header[1:].split('|')
                if len(parts) >= 6:
                    hit_info = {
                        'name': parts[0],
                        'length': len(current_seq),
                        'contig': parts[2],
                        'start': parts[3],
                        'end': parts[4],
                        'strand': parts[5]
                    }
                    species_hits[species_name].append(hit_info)

    # Also check for noe.scores files to get e-values if available
    for species_dir in results_path.iterdir():
        if not species_dir.is_dir():
            continue

        species_name = species_dir.name.replace('_', ' ')
        noe_scores = species_dir / "noe.scores"

        if noe_scores.exists() and species_name in species_hits:
            evalues = {}
            with open(noe_scores, 'r') as f:
                header = f.readline()  # Skip header
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            seq_id = parts[0].split('|')[0]
                            evalue = parts[1]
                            evalues[seq_id] = evalue

            # Update e-values for hits
            for hit in species_hits[species_name]:
                if hit['name'] in evalues:
                    hit['noe_evalue'] = evalues[hit['name']]
                else:
                    hit['noe_evalue'] = 'N/A'

    return species_hits

def compare_coordinates(prev_hits, curr_hits):
    """Compare coordinates between previous and current hits"""
    # Create coordinate signatures for comparison
    prev_coords = set()
    curr_coords = set()

    for hit in prev_hits:
        # Create coordinate signature: contig:start-end:strand
        coord_sig = f"{hit['contig']}:{hit['start']}-{hit['end']}:{hit['strand']}"
        prev_coords.add(coord_sig)

    for hit in curr_hits:
        coord_sig = f"{hit['contig']}:{hit['start']}-{hit['end']}:{hit['strand']}"
        curr_coords.add(coord_sig)

    # Find overlaps and differences
    same_coords = prev_coords & curr_coords
    only_prev = prev_coords - curr_coords
    only_curr = curr_coords - prev_coords

    return {
        'same': same_coords,
        'only_previous': only_prev,
        'only_current': only_curr,
        'exact_match': len(prev_coords) == len(curr_coords) and prev_coords == curr_coords
    }

def compare_results(previous, current, taxonomy=None):
    """Compare previous and current results"""

    if taxonomy is None:
        taxonomy = {}

    all_species = set(previous.keys()) | set(current.keys())

    # Categories of comparison
    same_count_same_coords = []
    same_count_diff_coords = []
    more_hits = []
    fewer_hits = []
    no_current_data = []
    new_species = []

    coordinate_details = {}  # Store detailed coordinate comparisons

    for species in sorted(all_species):
        prev_hits = previous.get(species, [])
        curr_hits = current.get(species, [])

        if species not in previous:
            # New species in current results
            new_species.append((species, len(curr_hits)))
        elif species not in current:
            # No current data for this species
            no_current_data.append((species, len(prev_hits)))
        else:
            # Compare hit counts
            prev_count = len(prev_hits)
            curr_count = len(curr_hits)

            if curr_count == prev_count:
                # Same count - now check coordinates
                coord_comparison = compare_coordinates(prev_hits, curr_hits)
                coordinate_details[species] = coord_comparison

                if coord_comparison['exact_match']:
                    same_count_same_coords.append((species, curr_count))
                else:
                    same_count_diff_coords.append((species, curr_count, coord_comparison))
            elif curr_count > prev_count:
                coord_comparison = compare_coordinates(prev_hits, curr_hits)
                coordinate_details[species] = coord_comparison
                more_hits.append((species, prev_count, curr_count))
            else:
                coord_comparison = compare_coordinates(prev_hits, curr_hits)
                coordinate_details[species] = coord_comparison
                fewer_hits.append((species, prev_count, curr_count))

    # Print comparison results
    print("=" * 80)
    print("NOECR34335 RESULTS COMPARISON")
    print("=" * 80)
    print(f"\nTotal species in previous results: {len(previous)}")
    print(f"Total species in current results: {len(current)}")
    print()

    # Species with same number of hits AND same coordinates
    if same_count_same_coords:
        print(f"✅ SAME COUNT & SAME COORDINATES ({len(same_count_same_coords)} species):")
        print("-" * 50)
        for species, count in same_count_same_coords:
            print(f"  {species}: {count} hits (identical coordinates)")

    # Species with same number of hits BUT different coordinates
    if same_count_diff_coords:
        print(f"\n⚠️  SAME COUNT & DIFFERENT COORDINATES ({len(same_count_diff_coords)} species):")
        print("-" * 50)
        for species, count, coord_comp in same_count_diff_coords:
            same_coord_count = len(coord_comp['same'])
            only_prev_count = len(coord_comp['only_previous'])
            only_curr_count = len(coord_comp['only_current'])

            print(f"  {species}: {count} hits")
            print(f"    ├─ Same coordinates: {same_coord_count}")
            print(f"    ├─ Only in previous: {only_prev_count}")
            print(f"    └─ Only in current:  {only_curr_count}")

            # Show specific coordinate differences if there are any
            if only_prev_count > 0:
                for coord in sorted(coord_comp['only_previous']):
                    print(f"      - Lost: {coord}")
            if only_curr_count > 0:
                for coord in sorted(coord_comp['only_current']):
                    print(f"      + Gained: {coord}")

    # Species with more hits
    if more_hits:
        print(f"\n📈 MORE HITS IN CURRENT ({len(more_hits)} species):")
        print("-" * 40)
        for species, prev, curr in more_hits:
            print(f"  {species}: {prev} → {curr} hits (+{curr-prev})")

    # Species with fewer hits
    if fewer_hits:
        print(f"\n📉 FEWER HITS IN CURRENT ({len(fewer_hits)} species):")
        print("-" * 40)
        for species, prev, curr in fewer_hits:
            print(f"  {species}: {prev} → {curr} hits (-{prev-curr})")

    # New species in current results
    if new_species:
        print(f"\n🆕 NEW SPECIES IN CURRENT ({len(new_species)} species):")
        print("-" * 40)
        for species, count in new_species:
            lineage = taxonomy.get(species, 'Unknown taxonomy')
            print(f"  {species}: {count} hits")
            print(f"    Lineage: {lineage}")

    # Species with no current data
    if no_current_data:
        print(f"\n⚠️  NO CURRENT DATA ({len(no_current_data)} species):")
        print("-" * 40)
        for species, prev_count in no_current_data:
            clean_name = clean_species_name(species)
            lineage = taxonomy.get(species, 'Unknown taxonomy')
            print(f"  {species}: had {prev_count} hits (check results/{clean_name}/)")
            print(f"    Lineage: {lineage}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY:")
    print(f"  Same count & coordinates: {len(same_count_same_coords)} species")
    print(f"  Same count, diff coords:  {len(same_count_diff_coords)} species")
    print(f"  More hits:                {len(more_hits)} species")
    print(f"  Fewer hits:               {len(fewer_hits)} species")
    print(f"  New species:              {len(new_species)} species")
    print(f"  No current data:          {len(no_current_data)} species")

    # Calculate total hits
    prev_total = sum(len(hits) for hits in previous.values())
    curr_total = sum(len(hits) for hits in current.values())
    print(f"\n  Total previous hits: {prev_total}")
    print(f"  Total current hits:  {curr_total}")

    if curr_total > prev_total:
        print(f"  Change: +{curr_total - prev_total} hits overall")
    elif curr_total < prev_total:
        print(f"  Change: -{prev_total - curr_total} hits overall")
    else:
        print(f"  Change: No change in total hits")

def main():
    parser = argparse.ArgumentParser(
        description='Compare previous noeCR34335 results with current results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare using default paths
  python3 compare_noecr_results.py

  # Specify custom CSV file
  python3 compare_noecr_results.py --csv old_results.csv

  # Specify custom results directory
  python3 compare_noecr_results.py --results-dir /path/to/results

  # Verbose mode to see detailed hit information
  python3 compare_noecr_results.py --verbose
        """
    )

    parser.add_argument('--csv', default='results_summary.csv',
                       help='Path to previous results CSV file (default: results_summary.csv)')
    parser.add_argument('--results-dir', default='results',
                       help='Path to current results directory (default: results)')
    parser.add_argument('--verbose', action='store_true',
                       help='Show detailed hit information')

    args = parser.parse_args()

    # Check if CSV file exists
    if not Path(args.csv).exists():
        print(f"Error: CSV file not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    # Load taxonomy information
    print(f"Loading taxonomy information...")
    taxonomy = load_taxonomy()

    # Load previous results
    print(f"Loading previous results from {args.csv}...")
    previous = load_previous_results(args.csv)

    # Load current results
    print(f"Loading current results from {args.results_dir}/...")
    current = load_current_results(args.results_dir)

    # Compare results
    compare_results(previous, current, taxonomy)

    # Verbose output if requested
    if args.verbose:
        print("\n" + "=" * 80)
        print("DETAILED COMPARISON:")

        # Show details for species with different counts
        for species in sorted(set(previous.keys()) | set(current.keys())):
            prev_hits = previous.get(species, [])
            curr_hits = current.get(species, [])

            if len(prev_hits) != len(curr_hits) and species in previous and species in current:
                print(f"\n{species}:")
                print(f"  Previous hits ({len(prev_hits)}):")
                for hit in prev_hits:
                    print(f"    - {hit['name']} ({hit['length']}bp, {hit['contig']}:{hit['start']}-{hit['end']})")

                print(f"  Current hits ({len(curr_hits)}):")
                for hit in curr_hits:
                    print(f"    - {hit['name']} ({hit['length']}bp, {hit['contig']}:{hit['start']}-{hit['end']})")

if __name__ == "__main__":
    main()