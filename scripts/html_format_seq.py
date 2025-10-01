#!/usr/bin/env python3
"""
Format FASTA sequences with score information for HTML display
"""

import sys
import re
from Bio import SeqIO

def load_scores(score_file):
    """Load scores from TSV file"""
    scores = {}
    try:
        with open(score_file, 'r') as f:
            header = f.readline().strip()  # Skip header
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        seq_id = parts[0]
                        pol2_score = parts[2]
                        pol3_score = parts[3]
                        scores[seq_id] = (pol2_score, pol3_score)
    except:
        pass
    return scores

def load_noe_scores(noe_score_file):
    """Load noe e-values from TSV file"""
    noe_scores = {}
    try:
        with open(noe_score_file, 'r') as f:
            header = f.readline().strip()  # Skip header
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        # Parse the full sequence ID
                        full_id = parts[0]
                        evalue = parts[1]
                        # Extract just the short ID (before first |)
                        seq_id = full_id.split('|')[0]
                        noe_scores[seq_id] = evalue
    except:
        pass
    return noe_scores

def load_lncrna_stats(stats_file):
    """Load lncRNA sequence statistics from TSV file"""
    stats = {}
    try:
        with open(stats_file, 'r') as f:
            header = f.readline().strip()  # Skip header
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        # Parse the full sequence ID and stats
                        full_id = parts[0]
                        # Extract just the short ID (before first |)
                        seq_id = full_id.split('|')[0]
                        stats[seq_id] = {
                            'has_5prime_motif': parts[1],
                            'has_3prime_motif': parts[2],
                            'longest_polyt': parts[3],
                            'trailing_t': parts[4],
                            'motif_5prime_type': parts[5] if len(parts) > 5 else 'none',
                            'motif_3prime_type': parts[6] if len(parts) > 6 else 'none'
                        }
    except:
        pass
    return stats

def load_incomplete_ids(incomplete_file):
    """Load incomplete sequence IDs from file"""
    incomplete_ids = set()
    try:
        with open(incomplete_file, 'r') as f:
            for line in f:
                if line.strip():
                    incomplete_ids.add(line.strip())
    except:
        pass
    return incomplete_ids

def load_extended_ids(extended_file):
    """Load extended sequence IDs from file"""
    extended_ids = set()
    try:
        with open(extended_file, 'r') as f:
            for line in f:
                if line.strip():
                    extended_ids.add(line.strip())
    except:
        pass
    return extended_ids

def get_score_color(score_str):
    """Get color class based on score value"""
    try:
        score = float(score_str)
        if score >= 0.8:
            return "score-high"
        elif score >= 0.6:
            return "score-medium"
        else:
            return "score-low"
    except:
        return "score-low"

def format_evalue(evalue_str):
    """Format e-value for display"""
    try:
        evalue = float(evalue_str)
        if evalue == 0:
            return "0"
        elif evalue < 1e-100:
            return f"{evalue:.0e}"
        else:
            return f"{evalue:.1e}"
    except:
        return evalue_str

def color_nucleotides(seq):
    """Color each nucleotide with background colors using ClustalX color scheme"""
    colored_seq = ""

    # Use ClustalX color scheme
    for nt in seq.upper():
        if nt == 'A':
            colored_seq += f'<span style="background-color:#64f73f;">{nt}</span>'
        elif nt == 'T':
            colored_seq += f'<span style="background-color:#3c88ee;">{nt}</span>'
        elif nt == 'C':
            colored_seq += f'<span style="background-color:#ffb340;">{nt}</span>'
        elif nt == 'G':
            colored_seq += f'<span style="background-color:#eb413c;">{nt}</span>'
        else:
            colored_seq += nt  # Default for any non-standard nucleotide

    return colored_seq

def wrap_sequence_html(seq, line_length=150):
    """Wrap the sequence in lines with HTML formatting"""
    # Show the full sequence, wrapped at line_length
    wrapped_seq = ""
    for i in range(0, len(seq), line_length):
        segment = seq[i:i+line_length]
        colored_segment = color_nucleotides(segment)

        # If this is the last segment and it's less than 12 nt, add padding
        if i + line_length >= len(seq) and len(segment) < 12:
            # Add 6 non-breaking spaces for padding
            colored_segment += '&nbsp;' * 6

        wrapped_seq += colored_segment
        if i + line_length < len(seq):
            wrapped_seq += '\n'
    return wrapped_seq

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 html_format_seq.py <score_file> [incomplete_file] [noe_score_file] [stats_file] [extended_file]", file=sys.stderr)
        sys.exit(1)

    score_file = sys.argv[1]
    scores = load_scores(score_file)

    # Load incomplete sequence IDs if provided
    incomplete_ids = set()
    if len(sys.argv) >= 3:
        incomplete_file = sys.argv[2]
        incomplete_ids = load_incomplete_ids(incomplete_file)

    # Load noe scores if provided
    noe_scores = {}
    if len(sys.argv) >= 4:
        noe_score_file = sys.argv[3]
        noe_scores = load_noe_scores(noe_score_file)

    # Load lncRNA stats if provided
    lncrna_stats = {}
    if len(sys.argv) >= 5:
        stats_file = sys.argv[4]
        lncrna_stats = load_lncrna_stats(stats_file)

    # Load extended sequence IDs if provided
    extended_ids = set()
    if len(sys.argv) >= 6:
        extended_file = sys.argv[5]
        extended_ids = load_extended_ids(extended_file)

    # Read FASTA from stdin
    for record in SeqIO.parse(sys.stdin, "fasta"):
        # Extract sequence ID (before the |)
        seq_id = record.id.split('|')[0]
        id_parts = record.id.split('|')

        # Extract coordinates and info
        if len(id_parts) >= 6:
            gene_name = id_parts[1] if len(id_parts) > 1 else "noeCR34335"
            chrom = id_parts[2] if len(id_parts) > 2 else ""
            start = id_parts[3] if len(id_parts) > 3 else ""
            end = id_parts[4] if len(id_parts) > 4 else ""
            strand = id_parts[5] if len(id_parts) > 5 else ""

            # Start sequence entry container
            print('<div class="sequence-entry">')

            # Build sequence header with basic info only
            header_parts = [f"{seq_id} | {chrom}:{start}-{end} ({strand}) | {len(str(record.seq))} nt"]

            # Add NOE e-value if available
            if seq_id in noe_scores:
                noe_evalue = noe_scores[seq_id]
                formatted_evalue = format_evalue(noe_evalue)
                header_parts.append(f'lncRNA:noe consensus e-value: {formatted_evalue}')

            print(f'<div class="sequence-header">{" | ".join(header_parts)}</div>')

            # Build inline metadata
            metadata_items = []

            if seq_id in lncrna_stats:
                stats = lncrna_stats[seq_id]

                # 5' motif
                if stats['has_5prime_motif'] == '1':
                    motif_5prime_type = stats.get('motif_5prime_type', 'GCGGT')
                    metadata_items.append(f'<span class="metadata-item">5\' motif: {motif_5prime_type}</span>')
                else:
                    metadata_items.append('<span class="metadata-item">5\' motif: ✗</span>')

                # 3' motif
                if stats['has_3prime_motif'] == '1':
                    motif_3prime_type = stats.get('motif_3prime_type', 'ATCGC')
                    metadata_items.append(f'<span class="metadata-item">3\' motif: {motif_3prime_type}</span>')
                else:
                    metadata_items.append('<span class="metadata-item">3\' motif: ✗</span>')

                # Poly-T info
                metadata_items.append(f'<span class="metadata-item">Internal max. Poly-T: {stats["longest_polyt"]}nt</span>')
                metadata_items.append(f'<span class="metadata-item">Trailing-T: {stats["trailing_t"]}nt</span>')

                # PSE score for metadata line
                if seq_id in scores:
                    pol2_score, pol3_score = scores[seq_id]
                    metadata_items.append(f'<span class="metadata-item">PSE: {pol3_score}</span>')

            # Add extended marker if needed
            if seq_id in extended_ids:
                metadata_items.append('<span class="metadata-item" style="color:green;font-weight:bold;">EXTENDED</span>')

            # Add incomplete marker if needed
            if seq_id in incomplete_ids:
                metadata_items.append('<span class="metadata-item" style="color:red;font-weight:bold;">INCOMPLETE</span>')

            if metadata_items:
                print(f'<div class="sequence-metadata">{"".join(metadata_items)}</div>')

            # Print colored sequence
            colored_sequence = wrap_sequence_html(str(record.seq))
            print(f'<div class="sequence">{colored_sequence}</div>')

            print('</div>')  # Close sequence-entry

        else:
            # Fallback for malformed headers
            print('<div class="sequence-entry">')

            header_parts = [f"{seq_id} | {len(str(record.seq))} nt"]

            if seq_id in noe_scores:
                noe_evalue = noe_scores[seq_id]
                formatted_evalue = format_evalue(noe_evalue)
                header_parts.append(f'lncRNA:noe consensus e-value: {formatted_evalue}')

            print(f'<div class="sequence-header">{" | ".join(header_parts)}</div>')

            metadata_items = []
            if seq_id in extended_ids:
                metadata_items.append('<span class="metadata-item" style="color:green;font-weight:bold;">EXTENDED</span>')
            if seq_id in incomplete_ids:
                metadata_items.append('<span class="metadata-item" style="color:red;font-weight:bold;">INCOMPLETE</span>')

            if metadata_items:
                print(f'<div class="sequence-metadata">{"".join(metadata_items)}</div>')

            colored_sequence = wrap_sequence_html(str(record.seq))
            print(f'<div class="sequence">{colored_sequence}</div>')

            print('</div>')  # Close sequence-entry

if __name__ == "__main__":
    main()