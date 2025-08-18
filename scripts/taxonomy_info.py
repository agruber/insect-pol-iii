#!/usr/bin/env python3
"""
Taxonomic information database with English and German names
Loads data from taxonomy_info.json file
"""

import json
import os

# Global variable to cache the loaded data
_TAXONOMY_INFO = None

def _load_taxonomy_info():
    """Load taxonomy information from JSON file"""
    global _TAXONOMY_INFO
    
    if _TAXONOMY_INFO is None:
        # Try to find the JSON file relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        json_path = os.path.join(script_dir, '..', 'data', 'taxonomy_info.json')
        
        # Fallback paths
        fallback_paths = [
            os.path.join(script_dir, 'taxonomy_info.json'),
            os.path.join(os.getcwd(), 'data', 'taxonomy_info.json'),
            os.path.join(os.getcwd(), 'taxonomy_info.json')
        ]
        
        # Try to load from the primary path first
        for path in [json_path] + fallback_paths:
            if os.path.exists(path):
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        _TAXONOMY_INFO = json.load(f)
                    return _TAXONOMY_INFO
                except (json.JSONDecodeError, IOError) as e:
                    print(f"Warning: Could not load taxonomy info from {path}: {e}")
                    continue
        
        # If no file found, return empty dict
        print("Warning: No taxonomy_info.json file found, using empty taxonomy database")
        _TAXONOMY_INFO = {}
    
    return _TAXONOMY_INFO

def get_taxonomy_info(taxon_name):
    """Get English/German names and description for a taxon"""
    taxonomy_data = _load_taxonomy_info()
    return taxonomy_data.get(taxon_name, {
        'english': taxon_name,
        'german': taxon_name,
        'description': 'No information available'
    })

def format_taxon_display(taxon_name, species_count=None):
    """Format taxon name with additional information for display"""
    info = get_taxonomy_info(taxon_name)
    
    # Format: "Scientific name (English/German) (count)"
    display = f"{taxon_name}"
    if info['english'] != taxon_name or info['german'] != taxon_name:
        display += f" ({info['english']}/{info['german']})"
    
    if species_count is not None:
        display += f" ({species_count})"
    
    return display

def get_taxon_description(taxon_name):
    """Get description for a taxon"""
    info = get_taxonomy_info(taxon_name)
    return info['description']