#!/usr/bin/env python3
"""
Update file status cache for fast explore script startup
Scans genomes directory and stores status counts in JSON file
"""

import os
import json
import sys
from pathlib import Path

def clean_species_name(species_name):
    """Convert species name to filesystem-safe format"""
    return species_name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')

def scan_file_status(genomes_dir="genomes"):
    """Scan genomes directory for file status"""
    status_data = {}
    
    if not os.path.exists(genomes_dir):
        print(f"Warning: {genomes_dir} directory not found")
        return status_data
    
    print(f"Scanning {genomes_dir} directory for file status...")
    
    total_dirs = 0
    for species_name in os.listdir(genomes_dir):
        species_dir = os.path.join(genomes_dir, species_name)
        if os.path.isdir(species_dir):
            total_dirs += 1
    
    processed = 0
    for species_name in sorted(os.listdir(genomes_dir)):
        species_dir = os.path.join(genomes_dir, species_name)
        if os.path.isdir(species_dir):
            processed += 1
            if processed % 200 == 0 or processed == total_dirs:
                print(f"  Processed {processed}/{total_dirs} species directories...")
            
            # Check for genome, cmsearch results, and trnascan results
            genome_file = os.path.join(species_dir, "genome.fna.gz")
            search_results = os.path.join(species_dir, "cmsearch_results.txt.gz")
            trnascan_results = os.path.join(species_dir, "trnascan_results.gff.gz")
            
            status_data[species_name] = {
                'has_genome': os.path.exists(genome_file) and os.path.islink(genome_file),
                'has_search': os.path.exists(search_results),
                'has_trnascan': os.path.exists(trnascan_results)
            }
    
    return status_data

def save_status_cache(status_data, cache_file="data/file_status_cache.json"):
    """Save status data to cache file"""
    # Create data directory if it doesn't exist
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    
    # Add timestamp to cache
    import time
    cache_data = {
        'timestamp': time.time(),
        'status': status_data
    }
    
    with open(cache_file, 'w') as f:
        json.dump(cache_data, f, indent=2)
    
    print(f"Status cache saved to {cache_file}")

def main():
    genomes_dir = sys.argv[1] if len(sys.argv) > 1 else "genomes"
    cache_file = sys.argv[2] if len(sys.argv) > 2 else "data/file_status_cache.json"
    
    print("File Status Cache Update")
    print("="*40)
    print(f"Genomes directory: {genomes_dir}")
    print(f"Cache file: {cache_file}")
    print()
    
    # Scan file status
    status_data = scan_file_status(genomes_dir)
    
    # Calculate summary
    total_species = len(status_data)
    total_downloaded = sum(1 for s in status_data.values() if s['has_genome'])
    total_searched = sum(1 for s in status_data.values() if s['has_search'])
    total_trnascan = sum(1 for s in status_data.values() if s['has_trnascan'])
    
    print(f"\nSummary:")
    print(f"  Total species directories: {total_species}")
    print(f"  Downloaded genomes: {total_downloaded}")
    print(f"  cmsearch results: {total_searched}")
    print(f"  tRNAscan results: {total_trnascan}")
    
    # Save cache
    save_status_cache(status_data, cache_file)
    print(f"\nCache updated successfully!")

if __name__ == "__main__":
    main()