#!/usr/bin/env python3
"""
Interactive tree-based species explorer for genomes_annotated.tsv
Navigate with arrow keys, expand/collapse with +/-, select with Enter
"""

import sys
import pandas as pd
import argparse
from collections import defaultdict
import os
from taxonomy_info import format_taxon_display, get_taxon_description

try:
    import termios
    import tty
    import select
except ImportError:
    print("This script requires a Unix-like system with termios support")
    sys.exit(1)

class TreeNode:
    def __init__(self, name, level, parent=None):
        self.name = name
        self.level = level
        self.parent = parent
        self.children = {}
        self.species_indices = set()
        self.expanded = False
        self.is_species = False
        self.downloaded_count = 0
        self.searched_count = 0
        
    def add_species(self, species_idx, taxonomy_path):
        """Add a species following the taxonomy path"""
        if not taxonomy_path:
            self.species_indices.add(species_idx)
            self.is_species = True
            return
        
        next_level = taxonomy_path[0]
        if next_level not in self.children:
            self.children[next_level] = TreeNode(next_level, self.level + 1, self)
        
        self.children[next_level].add_species(species_idx, taxonomy_path[1:])
    
    def get_species_count(self):
        """Get total number of species in this subtree"""
        count = len(self.species_indices)
        for child in self.children.values():
            count += child.get_species_count()
        return count
    
    def update_status_counts(self, downloaded, searched):
        """Update download and search status counts for this node only"""
        self.downloaded_count += downloaded
        self.searched_count += searched
    
    def get_status_counts(self):
        """Get total downloaded and searched counts for this subtree"""
        downloaded = self.downloaded_count
        searched = self.searched_count
        
        for child in self.children.values():
            child_downloaded, child_searched = child.get_status_counts()
            downloaded += child_downloaded
            searched += child_searched
        
        return downloaded, searched

class SpeciesTreeExplorer:
    def __init__(self, tsv_file, start_taxon=None):
        self.df = pd.read_csv(tsv_file, sep='\t', header=None)
        self.levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        self.root = TreeNode("Root", -1)
        self.current_row = 0
        self.display_nodes = []
        self.selected_species = set()
        self.start_taxon = start_taxon or "Arthropoda"
        
        # Build the tree
        self.build_tree()
        self.check_file_status()
        self.auto_expand_to_start()
        self.update_display()
    
    def clean_species_name(self, species_name):
        """Convert species name to filesystem-safe format"""
        return species_name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')
    
    def build_tree(self):
        """Build the taxonomic tree from the dataframe"""
        print("Building taxonomic tree...")
        
        for idx, row in self.df.iterrows():
            # Extract taxonomy path: phylum, class, order, family, genus, species
            taxonomy = []
            
            # Add taxonomic levels (columns 3-10: phylum to genus)
            for col_idx in [3, 4, 5, 9, 10]:  # phylum, class, order, family, genus
                value = row.iloc[col_idx]
                if pd.notna(value) and value != 'NA':
                    taxonomy.append(str(value))
                else:
                    taxonomy.append(f"Unknown_{self.levels[len(taxonomy)]}")
            
            # Add species name (column 0)
            taxonomy.append(row.iloc[0])
            
            # Filter data based on start taxon if specified
            if self.start_taxon and self.start_taxon != "Root":
                # Check if this species belongs to the start taxon
                if not self._species_belongs_to_taxon(taxonomy, self.start_taxon):
                    continue
                # Trim taxonomy to start from the specified taxon
                taxonomy = self._trim_taxonomy_from_start(taxonomy, self.start_taxon)
            
            self.root.add_species(idx, taxonomy)
    
    def _trim_taxonomy_from_start(self, taxonomy_path, start_taxon):
        """Trim taxonomy path to start from the specified taxon"""
        try:
            start_index = taxonomy_path.index(start_taxon)
            # Return taxonomy starting from the specified taxon
            return taxonomy_path[start_index:]
        except ValueError:
            # If start taxon not found, return original path
            return taxonomy_path
    
    def _species_belongs_to_taxon(self, taxonomy_path, target_taxon):
        """Check if a species belongs to the target taxon"""
        return target_taxon in taxonomy_path
    
    def check_file_status(self):
        """Check which species have genomes downloaded and cmsearch results"""
        print("Checking genome and search status...")
        
        # Count all files in genomes directory
        total_downloaded = 0
        total_searched = 0
        
        if os.path.exists("genomes"):
            for species_name in os.listdir("genomes"):
                species_dir = f"genomes/{species_name}"
                if os.path.isdir(species_dir):
                    genome_file = f"{species_dir}/genome.fna.gz"
                    search_results = f"{species_dir}/cmsearch_results.txt.gz"
                    
                    # Check if genome is downloaded
                    has_genome = os.path.exists(genome_file) and os.path.islink(genome_file)
                    
                    # Check if cmsearch results exist
                    has_search = os.path.exists(search_results)
                    
                    if has_genome:
                        total_downloaded += 1
                    if has_search:
                        total_searched += 1
                    
                    # Mark in tree if either file exists
                    if has_genome or has_search:
                        self._mark_species_status(species_name, has_genome, has_search)
        
        print(f"Found {total_downloaded} downloaded genomes, {total_searched} search results")
    
    def _mark_species_status(self, species_name, has_genome, has_search):
        """Mark a species as downloaded/searched by finding it in the tree"""
        # Convert cleaned species name back to find in dataframe
        for idx, row in self.df.iterrows():
            if self.clean_species_name(row.iloc[0]) == species_name:
                # Found the species, now mark it in the tree
                self._update_species_status(idx, has_genome, has_search)
                break
    
    def _update_species_status(self, species_idx, has_genome, has_search):
        """Update status for a specific species in the tree"""
        row = self.df.iloc[species_idx]
        taxonomy = []
        
        # Build taxonomy path same as in build_tree
        for col_idx in [3, 4, 5, 9, 10]:  # phylum, class, order, family, genus
            value = row.iloc[col_idx]
            if pd.notna(value) and value != 'NA':
                taxonomy.append(str(value))
            else:
                taxonomy.append(f"Unknown_{self.levels[len(taxonomy)]}")
        
        taxonomy.append(row.iloc[0])  # species name
        
        # Apply same filtering as build_tree
        if self.start_taxon and self.start_taxon != "Root":
            if not self._species_belongs_to_taxon(taxonomy, self.start_taxon):
                return
            taxonomy = self._trim_taxonomy_from_start(taxonomy, self.start_taxon)
        
        # Find the species node and update its status
        current_node = self.root
        for taxon in taxonomy[:-1]:  # Navigate to parent of species
            if taxon in current_node.children:
                current_node = current_node.children[taxon]
            else:
                return  # Species not in tree
        
        # Update the species node (leaf node only, parent counts are calculated in get_status_counts)
        species_name = taxonomy[-1]
        if species_name in current_node.children:
            species_node = current_node.children[species_name]
            # Only update if not already counted (avoid double counting)
            if species_node.downloaded_count == 0 and species_node.searched_count == 0:
                downloaded_increment = 1 if has_genome else 0
                searched_increment = 1 if has_search else 0
                species_node.update_status_counts(downloaded_increment, searched_increment)
    
    def auto_expand_to_start(self):
        """Auto-expand tree to show the start taxon"""
        if self.start_taxon and self.start_taxon != "Root":
            # If we filtered to a specific taxon, expand the root to show its children
            self.root.expanded = True
        else:
            # Default behavior - find and expand to show the start taxon
            if not self._find_and_expand_taxon(self.root, self.start_taxon):
                # If start taxon not found, show available top-level taxa
                available_taxa = sorted(self.root.children.keys())
                print(f"Warning: Start taxon '{self.start_taxon}' not found.")
                print(f"Available top-level taxa: {', '.join(available_taxa)}")
                print(f"Starting from root instead.\n")
                self.start_taxon = None
    
    def _find_and_expand_taxon(self, node, target_taxon):
        """Recursively find and expand to target taxon"""
        # Check if target taxon is in this subtree
        if target_taxon in node.children:
            node.expanded = True
            # Set current row to the target taxon
            return True
        
        # Check children recursively
        for child in node.children.values():
            if self._find_and_expand_taxon(child, target_taxon):
                node.expanded = True
                return True
        
        return False
    
    def update_display(self):
        """Update the list of nodes to display"""
        self.display_nodes = []
        self._add_node_to_display(self.root, 0)
        
        # Only reset current_row during initial setup
        if not hasattr(self, '_initial_setup_done'):
            # Find the start taxon in display and set current row (only if not filtered)
            if self.start_taxon and self.start_taxon == "Root":
                for i, (node, indent) in enumerate(self.display_nodes):
                    if node.name == self.start_taxon:
                        self.current_row = i
                        break
            else:
                # If we filtered the tree, start at the first item
                self.current_row = 0
            self._initial_setup_done = True
        else:
            # Keep current_row within bounds after tree updates
            if self.current_row >= len(self.display_nodes):
                self.current_row = max(0, len(self.display_nodes) - 1)
    
    def _add_node_to_display(self, node, indent):
        """Recursively add nodes to display list"""
        if node.level >= 0:  # Skip root
            self.display_nodes.append((node, indent))
        
        if node.expanded or node.level < 0:
            for child_name in sorted(node.children.keys()):
                self._add_node_to_display(node.children[child_name], indent + 1)
    
    def get_key(self):
        """Get a single key press"""
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
            
            # Handle escape sequences (arrow keys)
            if ch == '\x1b':
                ch2 = sys.stdin.read(1)
                if ch2 == '[':
                    ch3 = sys.stdin.read(1)
                    if ch3 == 'A':
                        return 'UP'
                    elif ch3 == 'B':
                        return 'DOWN'
                    elif ch3 == 'C':
                        return 'RIGHT'
                    elif ch3 == 'D':
                        return 'LEFT'
            
            return ch
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    
    def draw_screen(self):
        """Draw the current tree view"""
        os.system('clear')
        
        print("Species Tree Explorer")
        if self.start_taxon and self.start_taxon != "Root":
            from taxonomy_info import get_taxonomy_info
            info = get_taxonomy_info(self.start_taxon)
            print(f"Filtered to: {self.start_taxon} ({info['english']}/{info['german']})")
        print("=" * 80)
        print("Navigation: ↑/↓ arrows, +/- expand/collapse, Enter select, 'i' info, 's' save, 'q' quit")
        
        # Show what will be saved
        current_node, _ = self.display_nodes[self.current_row] if self.display_nodes else (None, 0)
        if self.selected_species:
            print(f"Will save: {len(self.selected_species)} manually selected species")
        elif current_node and current_node.is_species:
            print(f"Will save: 1 species ('{current_node.name}')")
        elif current_node:
            count = self._get_all_species_count_in_subtree(current_node)
            print(f"Will save: {count} species from '{current_node.name}'")
        else:
            print("Will save: Nothing selected")
        print("-" * 80)
        
        # Display visible portion of tree
        start_row = max(0, self.current_row - 8)
        end_row = min(len(self.display_nodes), self.current_row + 12)
        
        for i in range(start_row, end_row):
            node, indent = self.display_nodes[i]
            prefix = "  " * indent
            
            # Current row indicator
            cursor = "→ " if i == self.current_row else "  "
            
            # Expansion indicator
            if node.children:
                expand_char = "- " if node.expanded else "+ "
            else:
                expand_char = "  "
            
            # Selection indicator for species
            selected = "✓ " if node.is_species and any(idx in self.selected_species for idx in node.species_indices) else ""
            
            # Format taxon name with translations and status counts
            if not node.is_species:
                total_count = node.get_species_count()
                downloaded, searched = node.get_status_counts()
                display_name = format_taxon_display(node.name, total_count)
                if downloaded > 0 or searched > 0:
                    display_name += f" [📥{downloaded} 🔍{searched}]"
            else:
                # For species nodes, show individual status
                species_name = self.clean_species_name(node.name)
                species_dir = f"genomes/{species_name}"
                has_genome = os.path.exists(f"{species_dir}/genome.fna.gz")
                has_search = os.path.exists(f"{species_dir}/cmsearch_results.txt.gz")
                
                status_icons = ""
                if has_genome:
                    status_icons += "📥"
                if has_search:
                    status_icons += "🔍"
                
                display_name = node.name
                if status_icons:
                    display_name += f" [{status_icons}]"
            
            line = f"{cursor}{prefix}{expand_char}{selected}{display_name}"
            
            # Highlight current row
            if i == self.current_row:
                print(f"\033[7m{line}\033[0m")  # Reverse video
            else:
                print(line)
        
        # Show description of current taxon
        if self.display_nodes:
            current_node, _ = self.display_nodes[self.current_row]
            if not current_node.is_species:
                description = get_taxon_description(current_node.name)
                print("-" * 80)
                print(f"Info: {description}")
        
        print("-" * 80)
        print("Controls: ↑/↓ navigate, + expand, - collapse, Enter select/deselect, 'i' info, 's' save, 'q' quit")
    
    def handle_key(self, key):
        """Handle keyboard input"""
        if key == 'UP' and self.current_row > 0:
            self.current_row -= 1
        elif key == 'DOWN' and self.current_row < len(self.display_nodes) - 1:
            self.current_row += 1
        elif key == '+':
            node, _ = self.display_nodes[self.current_row]
            if node.children and not node.expanded:
                node.expanded = True
                self.update_display()
        elif key == '-':
            node, _ = self.display_nodes[self.current_row]
            if node.expanded:
                node.expanded = False
                self.update_display()
        elif key == '\r' or key == '\n':  # Enter
            node, _ = self.display_nodes[self.current_row]
            if node.is_species:
                # Toggle selection for this species
                for idx in node.species_indices:
                    if idx in self.selected_species:
                        self.selected_species.remove(idx)
                    else:
                        self.selected_species.add(idx)
            else:
                # Toggle expansion
                node.expanded = not node.expanded
                self.update_display()
        elif key == 's':
            self.save_selection()
        elif key == 'i':
            self.show_detailed_info()
        elif key == 'q':
            return False
        
        return True
    
    def save_selection(self):
        """Save currently selected taxon or manually selected species to files"""
        current_node, _ = self.display_nodes[self.current_row]
        
        # Determine what to save and filename
        if self.selected_species:
            # Save manually selected species
            species_indices = self.selected_species
            filename_base = "manually_selected"
            save_type = "manually selected species"
        elif current_node.is_species:
            # Save current species
            species_indices = current_node.species_indices
            filename_base = self.clean_species_name(current_node.name)
            save_type = f"species '{current_node.name}'"
        else:
            # Save all species in current taxonomic group
            species_indices = self._get_all_species_in_subtree(current_node)
            filename_base = self.clean_species_name(current_node.name)
            save_type = f"all species in '{current_node.name}'"
        
        if not species_indices:
            input("No species to save. Press Enter to continue...")
            return
        
        selected_df = self.df.loc[list(species_indices)]
        
        # Create filenames based on selection
        tsv_file = f"{filename_base}.txt"
        names_file = f"{filename_base}_names.txt"
        
        # Save full TSV format
        with open(tsv_file, 'w') as f:
            for _, row in selected_df.iterrows():
                f.write('\t'.join(map(str, row.values)) + '\n')
        
        # Save cleaned names format
        with open(names_file, 'w') as f:
            for _, row in selected_df.iterrows():
                f.write(self.clean_species_name(row.iloc[0]) + '\n')
        
        print(f"\nSaved {len(selected_df)} {save_type} to:")
        print(f"  - {tsv_file} (TSV format)")
        print(f"  - {names_file} (names only)")
        input("Press Enter to continue...")
    
    def _get_all_species_in_subtree(self, node):
        """Get all species indices in this node and all its children"""
        species_indices = set(node.species_indices)
        
        for child in node.children.values():
            species_indices.update(self._get_all_species_in_subtree(child))
        
        return species_indices
    
    def _get_all_species_count_in_subtree(self, node):
        """Get count of all species in this node and all its children"""
        return len(self._get_all_species_in_subtree(node))
    
    def show_detailed_info(self):
        """Show detailed information about the current taxon"""
        if not self.display_nodes:
            return
        
        current_node, _ = self.display_nodes[self.current_row]
        
        os.system('clear')
        print("Detailed Taxonomic Information")
        print("=" * 80)
        
        if current_node.is_species:
            # Show species information
            print(f"Species: {current_node.name}")
            
            # Get the first species entry to show additional data
            if current_node.species_indices:
                idx = list(current_node.species_indices)[0]
                row = self.df.iloc[idx]
                
                print(f"Assembly Type: {row.iloc[1]}")
                print(f"Download URL: {row.iloc[2]}")
                print("\nTaxonomic Classification:")
                print(f"  Phylum: {row.iloc[3]}")
                print(f"  Class: {row.iloc[4]}")
                print(f"  Order: {row.iloc[5]}")
                print(f"  Family: {row.iloc[9]}")
                print(f"  Genus: {row.iloc[10]}")
        else:
            # Show taxonomic group information
            from taxonomy_info import get_taxonomy_info
            info = get_taxonomy_info(current_node.name)
            
            print(f"Taxon: {current_node.name}")
            print(f"English: {info['english']}")
            print(f"German: {info['german']}")
            print(f"Species count: {current_node.get_species_count()}")
            
            # Show status counts
            downloaded, searched = current_node.get_status_counts()
            print(f"📥 Downloaded genomes: {downloaded}")
            print(f"🔍 Completed searches: {searched}")
            if current_node.get_species_count() > 0:
                download_percent = (downloaded / current_node.get_species_count()) * 100
                search_percent = (searched / current_node.get_species_count()) * 100
                print(f"Download progress: {download_percent:.1f}%")
                print(f"Search progress: {search_percent:.1f}%")
            
            print(f"\nDescription:")
            print(f"  {info['description']}")
            
            # Show some example species
            if current_node.species_indices:
                example_indices = list(current_node.species_indices)[:5]
                print(f"\nExample species:")
                for idx in example_indices:
                    species_name = self.df.iloc[idx, 0]
                    print(f"  • {species_name}")
                
                if len(current_node.species_indices) > 5:
                    print(f"  ... and {len(current_node.species_indices) - 5} more")
        
        print("\n" + "=" * 80)
        input("Press Enter to return...")
    
    def run(self):
        """Main interactive loop"""
        try:
            while True:
                self.draw_screen()
                key = self.get_key()
                if not self.handle_key(key):
                    break
        except KeyboardInterrupt:
            pass
        
        print("\nGoodbye!")

def main():
    parser = argparse.ArgumentParser(description='Interactive tree-based species explorer')
    parser.add_argument('tsv_file', help='Path to the genomes TSV file')
    parser.add_argument('-s', '--start', default='Arthropoda', 
                       help='Starting taxon to focus on (default: Arthropoda)')
    
    args = parser.parse_args()
    
    try:
        explorer = SpeciesTreeExplorer(args.tsv_file, args.start)
        explorer.run()
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()