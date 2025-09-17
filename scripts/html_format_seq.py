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

def get_score_color(score_str):
    """Get color class based on score value"""
    try:
        score = float(score_str)
        if score >= 0.8:
            return "label-green"
        elif score >= 0.6:
            return "label-yellow"
        else:
            return "label-red"
    except:
        return "label-red"

def color_nucleotides(seq):
    """Color each nucleotide with specific background colors"""
    colored_seq = ""
    colors = {
        'A': '#64f73f',
        'C': '#ffb340',
        'G': '#eb413c',
        'T': '#3c88ee'
    }

    # Color each nucleotide by wrapping it in a span with the appropriate background color
    for nt in seq.upper():
        color = colors.get(nt, 'grey')  # Default to grey for any non-standard nucleotide
        colored_seq += f"<span style='background-color:{color};'>{nt}</span>"

    return colored_seq

def wrap_sequence_html(seq, line_length=150):
    """Wrap the sequence in lines with HTML formatting"""
    wrapped_seq = ""
    for i in range(0, len(seq), line_length):
        # Apply coloring to the current line segment
        colored_segment = color_nucleotides(seq[i:i+line_length])
        wrapped_seq += colored_segment + '<br>'
    return wrapped_seq

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 html_format_seq.py <score_file> [incomplete_file] [noe_score_file]", file=sys.stderr)
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

            # Print header with coordinate info and Pol III score
            header_text = f'{seq_id} {chrom}:{start}-{end} ({strand}) {len(str(record.seq))}nt'

            if seq_id in scores:
                pol2_score, pol3_score = scores[seq_id]
                pol3_class = get_score_color(pol3_score)
                header_text += f' | Pol III PWM match score: <span class="{pol3_class}">{pol3_score}</span>'

            # Add noe e-value if available
            if seq_id in noe_scores:
                noe_evalue = noe_scores[seq_id]
                header_text += f' | noe e-value: {noe_evalue}'

            # Add incomplete marker if sequence is incomplete
            if seq_id in incomplete_ids:
                header_text += ' | <span style="color:red;font-weight:bold;">INCOMPLETE (no poly-T termination)</span>'

            print(f'<div style="margin-bottom:4px;"><b>{header_text}</b></div>')

            # Print colored sequence wrapped in lines
            colored_sequence = wrap_sequence_html(str(record.seq))
            print(colored_sequence)
            print("<br>")
        else:
            # Fallback for malformed headers
            header_text = f">{record.id}"
            if seq_id in scores:
                pol2_score, pol3_score = scores[seq_id]
                pol3_class = get_score_color(pol3_score)
                header_text += f" | Pol III PWM match score: <span class=\"{pol3_class}\">{pol3_score}</span>"

            # Add noe e-value if available
            if seq_id in noe_scores:
                noe_evalue = noe_scores[seq_id]
                header_text += f' | noe e-value: {noe_evalue}'

            # Add incomplete marker if sequence is incomplete
            if seq_id in incomplete_ids:
                header_text += ' | <span style="color:red;font-weight:bold;">INCOMPLETE (no poly-T termination)</span>'

            print(f'<div style="margin-bottom:4px;"><b>{header_text}</b></div>')

            # Print colored sequence
            colored_sequence = wrap_sequence_html(str(record.seq))
            print(colored_sequence)
            print("<br>")

if __name__ == "__main__":
    main()