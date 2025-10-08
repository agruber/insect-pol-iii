#!/usr/bin/env python3
"""
Generate a taxonomic tree showing lncRNA hit distribution.
Shows hierarchy from Phylum down to Family level with hit rates.
"""

import sys
import os
from pathlib import Path
from collections import defaultdict

def normalize_name(species_name):
    """Convert species name to directory name format"""
    return species_name.replace(' ', '_')

def check_lncrna_hit(species_name, results_dir):
    """Check if species has non-empty lncrna.fa file"""
    species_dir = results_dir / normalize_name(species_name)
    lncrna_file = species_dir / 'lncrna.fa'

    if lncrna_file.exists() and lncrna_file.stat().st_size > 0:
        # Check if file has actual sequences (not just empty)
        with open(lncrna_file, 'r') as f:
            content = f.read().strip()
            # Has content and has at least one sequence header
            return len(content) > 0 and '>' in content
    return False

def build_taxonomy_trees(tsv_file, results_dir):
    """
    Build two taxonomy trees: one with hits, one with all species.
    Returns (hit_tree, total_tree)
    """
    # Tree structure: phylum -> class -> order -> suborder -> infraorder -> superfamily -> family -> [species]
    hit_tree = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))))
    total_tree = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))))

    with open(tsv_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            species = parts[0]
            phylum = parts[3] if parts[3] != 'NA' else 'Unknown_Phylum'
            class_ = parts[4] if parts[4] != 'NA' else 'Unknown_Class'
            order = parts[5] if parts[5] != 'NA' else 'Unknown_Order'
            suborder = parts[6] if parts[6] != 'NA' else None
            infraorder = parts[7] if parts[7] != 'NA' else None
            superfamily = parts[8] if parts[8] != 'NA' else None
            family = parts[9] if parts[9] != 'NA' else None

            # Check if species has lncRNA hit
            has_hit = check_lncrna_hit(species, results_dir)

            # Always add to total tree
            total_tree[phylum][class_][order][suborder][infraorder][superfamily][family].append(species)

            # Only add to hit tree if species has a hit
            if has_hit:
                hit_tree[phylum][class_][order][suborder][infraorder][superfamily][family].append(species)

    return hit_tree, total_tree

def count_species_recursive(node):
    """Recursively count species in a tree node"""
    if isinstance(node, list):
        return len(node)
    elif isinstance(node, dict):
        return sum(count_species_recursive(v) for v in node.values())
    return 0

def print_tree(hit_tree, total_tree):
    """Print the taxonomy tree with hit/total species counts, starting at Brachycera level"""

    # Navigate directly to Brachycera within Arthropoda > Insecta > Diptera
    if 'Arthropoda' not in hit_tree:
        print("No data found for Arthropoda")
        return

    hit_arthropoda = hit_tree['Arthropoda']
    total_arthropoda = total_tree['Arthropoda']

    if 'Insecta' not in hit_arthropoda:
        print("No data found for Insecta")
        return

    hit_insecta = hit_arthropoda['Insecta']
    total_insecta = total_arthropoda['Insecta']

    if 'Diptera' not in hit_insecta:
        print("No data found for Diptera")
        return

    hit_diptera = hit_insecta['Diptera']
    total_diptera = total_insecta['Diptera']

    # Find Brachycera suborder
    suborders = hit_diptera

    # Brachycera is a suborder, so we need to navigate through the structure
    # The structure is: order -> suborder -> infraorder -> superfamily -> family
    # Since we're already at Diptera (order level), we can access suborders directly
    suborders_sorted = sorted(
        [s for s in suborders.keys() if s is not None],
        key=lambda s: count_species_recursive(suborders[s]),
        reverse=True
    )

    for i_sub, suborder in enumerate(suborders_sorted):
        if suborder != 'Brachycera':
            continue

        sub_hit = count_species_recursive(suborders[suborder])
        sub_total = count_species_recursive(total_diptera[suborder])

        print(f"Brachycera ({sub_hit}/{sub_total})")

        hit_infraorders = suborders[suborder]
        total_infraorders = total_diptera[suborder]

        # Process from infraorder level down
        # Get all infraorders from total tree (includes those with 0 hits)
        all_infraorders = [i for i in total_infraorders.keys() if i is not None]

        # Sort by hit count (0 if not in hit tree)
        infraorders_sorted = sorted(
            all_infraorders,
            key=lambda i: count_species_recursive(hit_infraorders.get(i, {})),
            reverse=True
        )

        for i_infra, infraorder in enumerate(infraorders_sorted):
            infra_hit = count_species_recursive(hit_infraorders.get(infraorder, {}))
            infra_total = count_species_recursive(total_infraorders[infraorder])
            is_last_infra = (i_infra == len(infraorders_sorted) - 1)
            infra_prefix = "└──" if is_last_infra else "├──"
            infra_continuation = "    " if is_last_infra else "│   "

            print(f"{infra_prefix} {infraorder} ({infra_hit}/{infra_total})")

            hit_superfamilies = hit_infraorders.get(infraorder, {})
            total_superfamilies = total_infraorders[infraorder]

            process_superfamily_level(
                hit_superfamilies, total_superfamilies,
                infra_continuation, False
            )

def process_infraorder_level(hit_infraorders, total_infraorders, prefix, is_last_parent):
    """Process infraorder level and below"""
    infraorders_sorted = sorted(
        [i for i in hit_infraorders.keys() if i is not None],
        key=lambda i: count_species_recursive(hit_infraorders[i]),
        reverse=True
    )

    # Handle direct superfamilies/families under previous level (no infraorder)
    if None in hit_infraorders:
        process_superfamily_level(
            hit_infraorders[None], total_infraorders[None],
            prefix, False
        )

    # Show infraorders
    for i_infra, infraorder in enumerate(infraorders_sorted):
        infra_hit = count_species_recursive(hit_infraorders[infraorder])
        infra_total = count_species_recursive(total_infraorders[infraorder])
        is_last_infra = (i_infra == len(infraorders_sorted) - 1)
        infra_prefix = "└──" if is_last_infra else "├──"
        infra_continuation = "    " if is_last_infra else "│   "

        print(f"│   {prefix}{infra_prefix} {infraorder} ({infra_hit}/{infra_total})")

        hit_superfamilies = hit_infraorders[infraorder]
        total_superfamilies = total_infraorders[infraorder]

        process_superfamily_level(
            hit_superfamilies, total_superfamilies,
            prefix + infra_continuation, False
        )

def process_superfamily_level(hit_superfamilies, total_superfamilies, prefix, is_last_parent):
    """Process superfamily level and below"""
    # Get all superfamilies from total tree (includes those with 0 hits)
    all_superfamilies = [s for s in total_superfamilies.keys() if s is not None]

    # Sort by hit count (0 if not in hit tree)
    superfamilies_sorted = sorted(
        all_superfamilies,
        key=lambda s: count_species_recursive(hit_superfamilies.get(s, {})),
        reverse=True
    )

    # Handle direct families under previous level (no superfamily)
    if None in total_superfamilies:
        process_family_level(
            hit_superfamilies.get(None, {}), total_superfamilies[None],
            prefix, False
        )

    # Show superfamilies
    for i_sf, superfamily in enumerate(superfamilies_sorted):
        sf_hit = count_species_recursive(hit_superfamilies.get(superfamily, {}))
        sf_total = count_species_recursive(total_superfamilies[superfamily])
        is_last_sf = (i_sf == len(superfamilies_sorted) - 1)
        sf_prefix = "└──" if is_last_sf else "├──"
        sf_continuation = "    " if is_last_sf else "│   "

        print(f"{prefix}{sf_prefix} {superfamily} ({sf_hit}/{sf_total})")

        hit_families = hit_superfamilies.get(superfamily, {})
        total_families = total_superfamilies[superfamily]

        process_family_level(
            hit_families, total_families,
            prefix + sf_continuation, False
        )

def process_family_level(hit_families, total_families, prefix, is_last_parent):
    """Process family level (deepest level shown)"""
    # Get all families from total tree (includes those with 0 hits)
    all_families = [f for f in total_families.keys() if f is not None]

    # Sort by hit count (0 if not in hit tree)
    families_sorted = sorted(
        all_families,
        key=lambda f: count_species_recursive(hit_families.get(f, {})),
        reverse=True
    )

    # Show families
    for i_fam, family in enumerate(families_sorted):
        fam_hit = count_species_recursive(hit_families.get(family, {}))
        fam_total = count_species_recursive(total_families[family])
        is_last_fam = (i_fam == len(families_sorted) - 1)
        fam_prefix = "└──" if is_last_fam else "├──"

        print(f"{prefix}{fam_prefix} {family} ({fam_hit}/{fam_total})")

def main():
    # Paths
    base_dir = Path(__file__).parent.parent
    tsv_file = base_dir / 'data' / 'genomes_annotated.tsv'
    results_dir = base_dir / 'results'

    if not tsv_file.exists():
        print(f"Error: {tsv_file} not found", file=sys.stderr)
        sys.exit(1)

    if not results_dir.exists():
        print(f"Error: {results_dir} not found", file=sys.stderr)
        sys.exit(1)

    # Build trees
    hit_tree, total_tree = build_taxonomy_trees(tsv_file, results_dir)

    # Calculate totals
    total_hits = count_species_recursive(hit_tree)
    total_species = count_species_recursive(total_tree)

    print(f"\nnoeCR34335 lncRNA Distribution Tree")
    print(f"{'='*50}")
    print(f"Species with lncRNA hits: {total_hits}/{total_species}")
    print(f"Format: Taxon (hits/total)\n")

    # Print tree
    print_tree(hit_tree, total_tree)
    print()

if __name__ == '__main__':
    main()
