#!/usr/bin/env python3
"""
Process BLAST output to keep only the best non-overlapping hits
Removes hits that are contained within longer hits
For equally long hits, keeps the one with the lowest e-value
"""

import sys
import argparse
from typing import List, Tuple, Dict
from dataclasses import dataclass

@dataclass
class BlastHit:
    """Represents a BLAST hit with all relevant information"""
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    sstrand: str
    original_line: str
    
    def __post_init__(self):
        """Ensure start <= end for subject coordinates"""
        if self.sstart > self.send:
            self.sstart, self.send = self.send, self.sstart
    
    @classmethod
    def from_line(cls, line: str):
        """Create BlastHit from BLAST output line"""
        fields = line.strip().split('\t')
        if len(fields) < 15:
            return None
        
        try:
            return cls(
                qseqid=fields[0],
                sseqid=fields[1],
                pident=float(fields[2]),
                length=int(fields[3]),
                mismatch=int(fields[4]),
                gapopen=int(fields[5]),
                qstart=int(fields[6]),
                qend=int(fields[7]),
                sstart=int(fields[8]),
                send=int(fields[9]),
                evalue=float(fields[10]),
                bitscore=float(fields[11]),
                qlen=int(fields[12]),
                slen=int(fields[13]),
                sstrand=fields[14] if len(fields) > 14 else 'plus',
                original_line=line.strip()
            )
        except (ValueError, IndexError):
            return None
    
    def overlaps_with(self, other: 'BlastHit') -> bool:
        """Check if this hit overlaps with another hit on the same sequence"""
        if self.sseqid != other.sseqid:
            return False
        
        # Check for overlap in subject coordinates
        return not (self.send < other.sstart or self.sstart > other.send)
    
    def contains(self, other: 'BlastHit') -> bool:
        """Check if this hit completely contains another hit"""
        if self.sseqid != other.sseqid:
            return False
        
        return (self.sstart <= other.sstart and self.send >= other.send)
    
    def subject_length(self) -> int:
        """Return the length of the hit on the subject sequence"""
        return abs(self.send - self.sstart) + 1

def process_blast_hits(hits: List[BlastHit]) -> List[BlastHit]:
    """
    Process BLAST hits using greedy algorithm with prioritization:
    1. Prioritize full-length 100% identity matches first
    2. Then sort by length (longest first)  
    3. Then by e-value (lowest first) for ties
    4. Keep the first hit, discard all overlapping hits
    5. Continue with next non-overlapping hit
    """
    if not hits:
        return []
    
    # Sort with priority:
    # 1. Full-length 100% identity matches first
    # 2. Then by length (descending)
    # 3. Then by e-value (ascending)
    def sort_key(h):
        is_full_length_perfect = (h.pident == 100.0 and h.length == h.qlen)
        return (not is_full_length_perfect, -h.subject_length(), h.evalue)
    
    hits.sort(key=sort_key)
    
    kept_hits = []
    
    for hit in hits:
        # Check if this hit overlaps with any already kept hit
        overlaps = False
        for kept_hit in kept_hits:
            if hit.overlaps_with(kept_hit):
                overlaps = True
                break
        
        # If no overlap, keep this hit
        if not overlaps:
            kept_hits.append(hit)
    
    # Sort output by sequence ID and position for consistent output
    kept_hits.sort(key=lambda h: (h.sseqid, h.sstart))
    
    return kept_hits

def main():
    parser = argparse.ArgumentParser(
        description='Process BLAST output to keep only the best non-overlapping hits',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script processes BLAST tabular output to remove redundant hits:

1. Removes hits that are completely contained within longer hits
2. For overlapping hits of equal length, keeps the one with lowest e-value
3. Maintains the best representative hit for each genomic region

Examples:
  # Process BLAST output from file
  python3 process_blast_hits.py < results.blast > filtered_results.blast
  
  # In a pipeline
  zcat genome/species/query.blast.gz | python3 process_blast_hits.py > filtered.blast
  
  # With minimum e-value filtering
  zcat genome/species/query.blast.gz | python3 process_blast_hits.py -e 1e-10 > filtered.blast

Input format expected (BLAST format 6 with extended fields):
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand

The script reads from stdin and writes to stdout, preserving the original format.
        """
    )
    
    parser.add_argument('-e', '--evalue', type=float, default=None,
                       help='Filter hits with e-value worse than threshold before processing')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print processing statistics to stderr')
    
    args = parser.parse_args()
    
    # Read and parse all hits
    hits = []
    lines_read = 0
    
    for line in sys.stdin:
        line = line.strip()
        lines_read += 1
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            continue
        
        hit = BlastHit.from_line(line)
        if hit is None:
            if args.verbose:
                print(f"Warning: Could not parse line {lines_read}: {line}", file=sys.stderr)
            continue
        
        # Apply e-value filter if specified
        if args.evalue is not None and hit.evalue > args.evalue:
            continue
        
        hits.append(hit)
    
    if args.verbose:
        print(f"Read {len(hits)} valid hits from {lines_read} lines", file=sys.stderr)
    
    # Process hits to remove redundancy
    filtered_hits = process_blast_hits(hits)
    
    if args.verbose:
        removed = len(hits) - len(filtered_hits)
        print(f"Kept {len(filtered_hits)} hits, removed {removed} redundant hits", file=sys.stderr)
        if len(hits) > 0:
            print(f"Reduction: {removed/len(hits)*100:.1f}%", file=sys.stderr)
    
    # Output filtered hits
    for hit in filtered_hits:
        print(hit.original_line)

if __name__ == "__main__":
    main()