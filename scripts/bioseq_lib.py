#!/usr/bin/env python3
"""
Bioinformatics Sequence Processing Library

A comprehensive library for common bioinformatics operations including:
- FASTA sequence processing
- GFF3 file handling
- Genomic coordinate operations
- K-mer analysis and sequence operations

Author: Claude Code
"""

import sys
import gzip
import argparse
from typing import Dict, List, Tuple, Optional, TextIO, Union, Iterator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter


# =============================================================================
# File I/O Utilities
# =============================================================================

def smart_open(filename: Optional[str], mode: str = 'r') -> TextIO:
    """
    Open a file with automatic gzip detection and stdin/stdout handling.
    
    Args:
        filename: File path or None for stdin/stdout
        mode: File mode ('r', 'w', etc.)
        
    Returns:
        File handle
    """
    if filename is None:
        return sys.stdin if 'r' in mode else sys.stdout
    
    if filename.endswith('.gz'):
        if 'r' in mode:
            return gzip.open(filename, 'rt')
        else:
            return gzip.open(filename, 'wt')
    else:
        return open(filename, mode)


def progress_reporter(current: int, interval: int = 1000, message: str = "Processed") -> None:
    """
    Report progress to stderr at regular intervals.
    
    Args:
        current: Current count
        interval: Reporting interval
        message: Progress message prefix
    """
    if current % interval == 0:
        print(f"{message} {current}...", file=sys.stderr)


# =============================================================================
# FASTA Processing
# =============================================================================

def load_fasta_sequences(filename: str) -> Dict[str, SeqRecord]:
    """
    Load FASTA sequences into a dictionary for fast random access.
    
    Args:
        filename: FASTA file path (supports .gz compression)
        
    Returns:
        Dictionary mapping sequence ID to SeqRecord
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    print(f"Loading sequences from {filename}...", file=sys.stderr)
    
    try:
        with smart_open(filename, 'r') as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        
        seq_count = len(sequences)
        total_length = sum(len(seq.seq) for seq in sequences.values())
        print(f"Loaded {seq_count} sequences ({total_length:,} bp total)", file=sys.stderr)
        
        return sequences
        
    except Exception as e:
        print(f"Error loading FASTA file: {e}", file=sys.stderr)
        raise


def parse_fasta_from_handle(handle: TextIO) -> Iterator[SeqRecord]:
    """
    Parse FASTA sequences from an open file handle.
    
    Args:
        handle: Open file handle
        
    Yields:
        SeqRecord objects
    """
    return SeqIO.parse(handle, 'fasta')


def write_fasta_sequence(handle: TextIO, seq_id: str, sequence: Union[str, Seq], 
                        lowercase: bool = False, line_wrap: Optional[int] = None) -> None:
    """
    Write a single FASTA sequence to file handle.
    
    Args:
        handle: Output file handle
        seq_id: Sequence identifier
        sequence: Sequence string or Bio.Seq object
        lowercase: Convert to lowercase
        line_wrap: Wrap lines at this width (None = no wrapping)
    """
    handle.write(f">{seq_id}\n")
    
    seq_str = str(sequence)
    if lowercase:
        seq_str = seq_str.lower()
    
    if line_wrap and line_wrap > 0:
        for i in range(0, len(seq_str), line_wrap):
            handle.write(seq_str[i:i+line_wrap] + '\n')
    else:
        handle.write(seq_str + '\n')


def get_sequence_lengths(sequences: Dict[str, SeqRecord]) -> Dict[str, int]:
    """
    Extract sequence lengths from a sequence dictionary.
    
    Args:
        sequences: Dictionary of SeqRecord objects
        
    Returns:
        Dictionary mapping sequence ID to length
    """
    return {seq_id: len(record.seq) for seq_id, record in sequences.items()}


# =============================================================================
# GFF3 Processing
# =============================================================================

class GFF3Feature:
    """Represents a single GFF3 feature with easy attribute access."""
    
    def __init__(self, seqid: str, source: str, feature_type: str, start: int, end: int,
                 score: str = '.', strand: str = '.', phase: str = '.', attributes: str = '.'):
        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self._parsed_attributes = None
    
    @classmethod
    def from_line(cls, line: str) -> Optional['GFF3Feature']:
        """
        Create GFF3Feature from a GFF3 line.
        
        Args:
            line: GFF3 format line
            
        Returns:
            GFF3Feature object or None if line is invalid
        """
        if line.startswith('#') or not line.strip():
            return None
        
        parts = line.strip().split('\t')
        if len(parts) != 9:
            print(f"Warning: Malformed GFF3 line (expected 9 columns): {line}", file=sys.stderr)
            return None
        
        try:
            return cls(
                seqid=parts[0],
                source=parts[1], 
                feature_type=parts[2],
                start=int(parts[3]),
                end=int(parts[4]),
                score=parts[5],
                strand=parts[6],
                phase=parts[7],
                attributes=parts[8]
            )
        except (ValueError, IndexError) as e:
            print(f"Warning: Invalid GFF3 line: {line}, error: {e}", file=sys.stderr)
            return None
    
    def to_line(self) -> str:
        """Convert feature back to GFF3 format line."""
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{self.attributes}"
    
    def parse_attributes(self) -> Dict[str, str]:
        """
        Parse GFF3 attributes into a dictionary.
        
        Returns:
            Dictionary of attribute key-value pairs
        """
        if self._parsed_attributes is not None:
            return self._parsed_attributes
        
        self._parsed_attributes = {}
        if self.attributes and self.attributes != '.':
            for attr in self.attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    self._parsed_attributes[key] = value
        
        return self._parsed_attributes
    
    def get_attribute(self, key: str, default: Optional[str] = None) -> Optional[str]:
        """Get a specific attribute value."""
        return self.parse_attributes().get(key, default)
    
    def is_valid_coordinates(self, seq_lengths: Optional[Dict[str, int]] = None) -> bool:
        """
        Check if feature coordinates are valid.
        
        Args:
            seq_lengths: Optional sequence length constraints
            
        Returns:
            True if coordinates are valid
        """
        if self.start < 1 or self.start > self.end:
            return False
        
        if seq_lengths and self.seqid in seq_lengths:
            return self.end <= seq_lengths[self.seqid]
        
        return True


def parse_gff3_file(filename: Optional[str]) -> Iterator[GFF3Feature]:
    """
    Parse GFF3 file and yield features.
    
    Args:
        filename: GFF3 file path or None for stdin
        
    Yields:
        GFF3Feature objects
    """
    with smart_open(filename, 'r') as handle:
        for line in handle:
            feature = GFF3Feature.from_line(line)
            if feature:
                yield feature


# =============================================================================
# Genomic Coordinate Operations
# =============================================================================

def extract_genomic_sequence(feature: GFF3Feature, sequences: Dict[str, SeqRecord]) -> Optional[Seq]:
    """
    Extract genomic sequence for a GFF3 feature with strand handling.
    
    Args:
        feature: GFF3 feature
        sequences: Dictionary of genomic sequences
        
    Returns:
        Extracted sequence or None if extraction fails
    """
    if feature.seqid not in sequences:
        print(f"Warning: Sequence '{feature.seqid}' not found in genome", file=sys.stderr)
        return None
    
    genome_seq = sequences[feature.seqid].seq
    seq_length = len(genome_seq)
    
    # Validate coordinates
    if not feature.is_valid_coordinates({feature.seqid: seq_length}):
        print(f"Warning: Invalid coordinates {feature.start}-{feature.end} for sequence '{feature.seqid}' (length {seq_length})", file=sys.stderr)
        return None
    
    # Extract sequence (convert from 1-based GFF3 to 0-based Python)
    extracted_seq = genome_seq[feature.start-1:feature.end]
    
    # Handle strand orientation
    if feature.strand == '-':
        extracted_seq = extracted_seq.reverse_complement()
    elif feature.strand not in ['+', '.']:
        print(f"Warning: Unknown strand '{feature.strand}', treating as forward", file=sys.stderr)
    
    return extracted_seq


def calculate_upstream_region(feature: GFF3Feature, upstream_distance: int, 
                            seq_lengths: Dict[str, int], strict: bool = False) -> Optional[GFF3Feature]:
    """
    Calculate upstream region coordinates for a feature.
    
    Args:
        feature: Original GFF3 feature
        upstream_distance: Distance upstream to extract
        seq_lengths: Dictionary of sequence lengths
        strict: If True, only return regions that are exactly upstream_distance long
        
    Returns:
        New GFF3Feature for upstream region or None if invalid
    """
    if feature.seqid not in seq_lengths:
        print(f"Warning: No length information for sequence '{feature.seqid}'", file=sys.stderr)
        return None
    
    seq_length = seq_lengths[feature.seqid]
    
    # Calculate upstream coordinates based on strand
    if feature.strand == '+':
        # For + strand: upstream is before the start
        upstream_end = feature.start - 1
        upstream_start = max(1, upstream_end - upstream_distance + 1)
    elif feature.strand == '-':
        # For - strand: upstream is after the end  
        upstream_start = feature.end + 1
        upstream_end = min(seq_length, upstream_start + upstream_distance - 1)
    else:
        print(f"Warning: Unknown strand '{feature.strand}'", file=sys.stderr)
        return None
    
    # Validate region
    if upstream_start > upstream_end or upstream_start < 1 or upstream_end > seq_length:
        return None
    
    # Check strict mode: only return if region is exactly the requested length
    actual_length = upstream_end - upstream_start + 1
    if strict and actual_length != upstream_distance:
        return None
    
    # Create upstream feature
    return GFF3Feature(
        seqid=feature.seqid,
        source=feature.source,
        feature_type=feature.type,
        start=upstream_start,
        end=upstream_end,
        strand=feature.strand,
        attributes='.'
    )


# =============================================================================
# K-mer Analysis and Conservation
# =============================================================================

def extract_kmers(sequence: str, k: int, start_pos: int = 0, end_pos: Optional[int] = None) -> Dict[str, List[int]]:
    """
    Extract k-mers from a sequence region and their positions.
    
    Args:
        sequence: Input sequence
        k: K-mer length
        start_pos: Start position (0-based)
        end_pos: End position (0-based, None for end of sequence)
        
    Returns:
        Dictionary mapping k-mer to list of positions
    """
    if end_pos is None:
        end_pos = len(sequence)
    
    sequence = sequence.upper()
    kmers = {}
    
    for i in range(start_pos, end_pos - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:  # Skip k-mers with ambiguous nucleotides
            if kmer not in kmers:
                kmers[kmer] = []
            kmers[kmer].append(i)
    
    return kmers



def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two sequences.
    
    Args:
        seq1, seq2: Sequences to compare
        
    Returns:
        Hamming distance (number of differing positions)
    """
    if len(seq1) != len(seq2):
        return float('inf')
    
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))




# =============================================================================
# Command Line Utilities
# =============================================================================

def create_standard_parser(description: str, epilog: Optional[str] = None) -> argparse.ArgumentParser:
    """
    Create a standard argument parser with common formatting.
    
    Args:
        description: Program description
        epilog: Optional epilog text with examples
        
    Returns:
        Configured ArgumentParser
    """
    return argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog
    )


def add_io_arguments(parser: argparse.ArgumentParser, input_help: str, output_help: str,
                    require_input: bool = False) -> None:
    """
    Add standard input/output arguments to parser.
    
    Args:
        parser: ArgumentParser to modify
        input_help: Help text for input argument
        output_help: Help text for output argument
        require_input: Whether input argument is required
    """
    if require_input:
        parser.add_argument('input', help=input_help)
    else:
        parser.add_argument('input', nargs='?', help=input_help + ' (use stdin if not provided)')
    
    parser.add_argument('-o', '--output', help=output_help + ' (use stdout if not provided)')


# =============================================================================
# Validation and Error Handling
# =============================================================================

def validate_file_exists(filename: str) -> bool:
    """
    Check if a file exists and is readable.
    
    Args:
        filename: File path to check
        
    Returns:
        True if file exists and is readable
    """
    try:
        with smart_open(filename, 'r') as f:
            return True
    except (FileNotFoundError, PermissionError):
        return False


def validate_sequence_format(sequence: str, allowed_chars: str = "ATCGN-") -> bool:
    """
    Validate sequence contains only allowed characters.
    
    Args:
        sequence: Sequence to validate
        allowed_chars: String of allowed characters
        
    Returns:
        True if sequence is valid
    """
    return all(char.upper() in allowed_chars for char in sequence)


# =============================================================================
# Data Structures
# =============================================================================



if __name__ == "__main__":
    print("Bioinformatics Sequence Processing Library")
    print("This is a library module - import it in your scripts")
    print(f"Available functions: {[name for name in globals() if not name.startswith('_')]}")