#!/usr/bin/env python3
"""
Colorize alignment sequences for HTML display with ClustalW-like formatting
"""

import sys
import io
import colorsys
from Bio import AlignIO

def adjust_color(color, factor=0.5):
    """ Lighten the given color by the specified factor. """
    r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
    # Convert RGB to HLS, modify lightness, and convert back to RGB
    h, l, s = colorsys.rgb_to_hls(r/255.0, g/255.0, b/255.0)
    lighter = colorsys.hls_to_rgb(h, min(1, l * (1 + factor)), s=s)
    # Convert RGB back to hex
    return f"#{int(lighter[0]*255):02x}{int(lighter[1]*255):02x}{int(lighter[2]*255):02x}"

def colorize_sequence(sequence):
    colors = {
        'A': '#64f73f',
        'C': '#ffb340',
        'G': '#eb413c',
        'T': '#3c88ee'
    }
    colored_seq = ""
    for nt in sequence:
        char = nt.upper()
        if char in colors:
            color = colors[char]
            bold = "font-weight:bold;"
            if nt.islower():
                bold = ""
                color = adjust_color(color, factor=0.2)  # Lighten the color for lowercase letters
            colored_seq += f"<span style='{bold}background-color:{color};'>{nt}</span>"
        else:
            colored_seq += nt  # Non-ATCG characters are left unchanged
    return colored_seq

def format_sequence_name(seq_id):
    """Format sequence names - abbreviate noeCR34335, use ncRNA type for others"""
    parts = seq_id.split('|')
    if len(parts) >= 2:
        species_part = parts[0]
        gene_name = parts[1]

        if gene_name == 'noeCR34335':
            # For noeCR34335, just use the abbreviated species name
            return species_part
        else:
            # For others, use the ncRNA type (gene name)
            return gene_name
    else:
        # Fallback to original ID
        return seq_id

def main():
    # Read alignment from stdin into a string
    alignment_data = sys.stdin.read()

    # Use StringIO to create a file-like object from the string data
    file_like_data = io.StringIO(alignment_data)

    # Read the alignment using Biopython (assume FASTA format with gaps)
    alignment = AlignIO.read(file_like_data, "fasta")

    # Format sequence names and find maximum length for padding
    formatted_names = {}
    for record in alignment:
        formatted_names[record.id] = format_sequence_name(record.id)

    max_name_len = max(len(name) for name in formatted_names.values())

    # Colorize each sequence in the alignment and print in blocks
    block_width = 150
    num_blocks = (max(len(record.seq) for record in alignment) + block_width - 1) // block_width

    for block in range(num_blocks):
        for record in alignment:
            start = block * block_width
            end = start + block_width
            # Slice the sequence for the current block
            sequence_block = record.seq[start:end]
            # Colorize the sliced block
            colorized_sequence_block = colorize_sequence(sequence_block)
            # Get formatted name
            formatted_name = formatted_names[record.id]
            # Print each sequence name with aligned sequence block
            print(f"{formatted_name.ljust(max_name_len)}  {colorized_sequence_block}")
        # No vertical padding between blocks

if __name__ == "__main__":
    main()