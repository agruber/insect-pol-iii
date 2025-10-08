#!/usr/bin/env python3
"""
Analyze RNA secondary structure from RNAfold output.
Takes a FASTA file, runs RNAfold, and processes the structure.

Usage:
    python3 analyze_structure.py input.fa
    cat input.fa | python3 analyze_structure.py
"""

import sys
import subprocess
import argparse
from typing import List, Tuple, Optional


def classify_structure(structure: str) -> str:
    """Classify dot-bracket structure type"""
    if not structure:
        return "unpaired"

    open_count = structure.count('(')
    close_count = structure.count(')')

    # Must have equal opening and closing brackets
    if open_count != close_count:
        return "invalid"

    # If no brackets, it's just unpaired
    if open_count == 0:
        return "unpaired"

    # Check if structure is balanced
    if not is_balanced_structure(structure):
        return "invalid"

    # Check for multiple separate structures
    if has_separate_structures(structure):
        return "invalid"

    # Calculate max nesting depth
    depth = 0
    max_depth = 0
    for char in structure:
        if char == '(':
            depth += 1
            max_depth = max(max_depth, depth)
        elif char == ')':
            depth -= 1

    # Check for bulges
    has_bulges = has_bulges_in_structure(structure)

    # Classify based on depth and bulges
    if max_depth == 1:
        return "hairpin_with_bulges" if has_bulges else "hairpin"
    else:
        return "nested_with_bulges" if has_bulges else "nested"


def has_bulges_in_structure(structure: str) -> bool:
    """Check for bulges (dots at depth > 0)"""
    depth = 0
    dots_at_depth = set()

    for char in structure:
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
        elif char == '.':
            dots_at_depth.add(depth)

    # If we have dots at depth > 0, there are bulges
    return any(d > 0 for d in dots_at_depth)


def has_separate_structures(structure: str) -> bool:
    """Check for separate structures (disconnected hairpins)"""
    depth = 0
    found_complete_structure = False

    for i, char in enumerate(structure):
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
            if depth == 0:
                if found_complete_structure:
                    return True
                found_complete_structure = True
                # Check if there are more opening brackets ahead
                if '(' in structure[i+1:]:
                    return True

    return False


def is_balanced_structure(structure: str) -> bool:
    """Check if structure has proper bracket matching"""
    depth = 0
    for char in structure:
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def check_3prime_structure(region: str) -> bool:
    """Check 3' end structure rules"""
    # Rule 1: Simply 20 or fewer unpaired nucleotides (no brackets)
    if '(' not in region and ')' not in region:
        return len(region) <= 20

    # Rule 2: Check total unpaired nucleotides and structure type
    unpaired_count = region.count('.')
    structure_type = classify_structure(region)

    # Allow up to 20 nt unpaired AND/OR nested_with_bulges structure
    return (unpaired_count <= 20 and
            structure_type in ["hairpin", "hairpin_with_bulges", "nested_with_bulges"])


def check_region_rules(region: str, side: str) -> bool:
    """Check region rules for 5' or 3' end"""
    if side == "5prime":
        # 5' end: max 5 nt unpaired region allowed
        unpaired_count = region.count('.')
        return unpaired_count <= 5
    elif side == "3prime":
        # 3' end: allow hairpin structure with bulges OR simply 20 unpaired nucleotides
        return check_3prime_structure(region)

    return False


def find_base_pairs(structure: str) -> List[Tuple[int, int]]:
    """Find all base pairs using a stack"""
    stack = []
    pairs = []

    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                open_pos = stack.pop()
                pairs.append((open_pos, i))

    # Sort pairs by opening position
    pairs.sort(key=lambda x: x[0])
    return pairs


def process_structure(structure: str) -> str:
    """Process RNA structure and mark valid regions with 'x'"""
    chars = list(structure)
    pairs = find_base_pairs(structure)

    if not pairs:
        return structure

    # Start with the first base pair (leftmost opening bracket)
    first_pair = pairs[0]
    current_open, current_close = first_pair
    final_open = current_open
    final_close = current_close

    # Try to expand inward by finding nested base pairs
    expanded = True
    while expanded:
        expanded = False

        # Look for the next innermost pair
        for candidate_open, candidate_close in pairs:
            # Check if this pair is nested within our current final pair
            if candidate_open > final_open and candidate_close < final_close:

                # Check if there are any other pairs between current and candidate
                has_intermediate = False
                for other_open, other_close in pairs:
                    if (other_open > final_open and other_close < final_close and
                        other_open < candidate_open and other_close > candidate_close):
                        has_intermediate = True
                        break

                if not has_intermediate:
                    # Extract the 5' and 3' regions for rule checking
                    fivep_start = final_open + 1
                    fivep_end = candidate_open - 1
                    threep_start = candidate_close + 1
                    threep_end = final_close - 1

                    fivep_region = ""
                    threep_region = ""

                    if fivep_end >= fivep_start:
                        fivep_region = ''.join(chars[fivep_start:fivep_end + 1])

                    if threep_end >= threep_start:
                        threep_region = ''.join(chars[threep_start:threep_end + 1])

                    # Check both regions using the rules
                    fivep_valid = check_region_rules(fivep_region, "5prime")
                    threep_valid = check_region_rules(threep_region, "3prime")

                    # If both rules are satisfied, we can expand
                    if fivep_valid and threep_valid:
                        final_open = candidate_open
                        final_close = candidate_close
                        expanded = True
                        break  # Found our next level, continue expanding
                    else:
                        break  # Stop expanding - rules violated

    # Create result array
    result = list(chars)

    # Mark everything between the final expanded pair with 'x'
    for i in range(final_open + 1, final_close):
        result[i] = 'x'

    return ''.join(result)


def run_rnafold(fasta_input: str) -> List[Tuple[str, str, str]]:
    """
    Run RNAfold on FASTA input.
    Returns list of (header, sequence, structure) tuples.
    """
    try:
        # Run RNAfold with --noPS flag to avoid PostScript output
        result = subprocess.run(
            ['RNAfold', '--noPS'],
            input=fasta_input,
            capture_output=True,
            text=True,
            check=True
        )

        # Parse RNAfold output
        lines = result.stdout.strip().split('\n')
        results = []

        i = 0
        while i < len(lines):
            if lines[i].startswith('>'):
                header = lines[i]
                if i + 2 < len(lines):
                    sequence = lines[i + 1]
                    # Structure line contains structure and energy, extract just structure
                    structure_line = lines[i + 2]
                    # Structure is before the energy value in parentheses
                    structure = structure_line.split()[0]
                    results.append((header, sequence, structure))
                    i += 3
                else:
                    i += 1
            else:
                i += 1

        return results

    except subprocess.CalledProcessError as e:
        print(f"Error running RNAfold: {e}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: RNAfold not found. Please install ViennaRNA package.", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze RNA secondary structure using RNAfold',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From FASTA file
  python3 analyze_structure.py input.fa

  # From stdin
  cat input.fa | python3 analyze_structure.py

  # With specific input
  python3 analyze_structure.py -i sequences.fa -o structures.txt
        """
    )

    parser.add_argument('input', nargs='?',
                       help='Input FASTA file (default: stdin)')
    parser.add_argument('-i', '--input-file',
                       help='Input FASTA file (alternative to positional arg)')
    parser.add_argument('-o', '--output',
                       help='Output file (default: stdout)')

    args = parser.parse_args()

    # Determine input source
    input_file = args.input_file or args.input

    if input_file:
        try:
            with open(input_file, 'r') as f:
                fasta_input = f.read()
        except IOError as e:
            print(f"Error reading input file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Read from stdin
        fasta_input = sys.stdin.read()

    if not fasta_input.strip():
        print("Error: No input provided", file=sys.stderr)
        sys.exit(1)

    # Run RNAfold and get results
    results = run_rnafold(fasta_input)

    # Determine output destination
    if args.output:
        try:
            outfile = open(args.output, 'w')
        except IOError as e:
            print(f"Error opening output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        outfile = sys.stdout

    # Process each structure and output
    try:
        for header, sequence, structure in results:
            processed_structure = process_structure(structure)
            outfile.write(f"{header}\n")
            outfile.write(f"{structure}\n")
            outfile.write(f"{processed_structure}\n")
    finally:
        if args.output:
            outfile.close()


if __name__ == "__main__":
    main()
