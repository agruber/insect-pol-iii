#!/usr/bin/env python3
"""
Transfer files from Google Cloud instance to local machine
Takes a file with remote paths and copies them locally preserving directory structure
"""

import sys
import os
import subprocess
from pathlib import Path

def copy_files_from_instance(transfer_file, instance_name, zone, local_base_dir="genomes"):
    """
    Copy files from Google Cloud instance based on transfer list
    
    Args:
        transfer_file: File containing remote paths (one per line)
        instance_name: Name of the Google Cloud instance
        zone: Zone where the instance is located
        local_base_dir: Local directory to copy files to
    """
    
    if not os.path.exists(transfer_file):
        print(f"Error: Transfer file '{transfer_file}' not found")
        return False
    
    # Read the transfer list
    with open(transfer_file, 'r') as f:
        remote_paths = [line.strip() for line in f if line.strip()]
    
    if not remote_paths:
        print("No paths found in transfer file")
        return False
    
    print(f"Found {len(remote_paths)} files to transfer")
    
    # Create local base directory if it doesn't exist
    os.makedirs(local_base_dir, exist_ok=True)
    
    success_count = 0
    error_count = 0
    
    for remote_path in remote_paths:
        try:
            # Extract the relative path from the remote path
            # Expected format: /home/andreas/insect-pol-iii/genomes/Species_name/filename
            path_parts = Path(remote_path).parts
            
            # Find 'genomes' in the path and take everything after it
            try:
                genomes_index = path_parts.index('genomes')
                relative_path = Path(*path_parts[genomes_index+1:])  # Skip 'genomes' itself
            except ValueError:
                print(f"Warning: 'genomes' not found in path {remote_path}, using full filename")
                relative_path = Path(remote_path).name
            
            # Create local directory structure
            local_file_path = Path(local_base_dir) / relative_path
            local_file_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Build gcloud scp command
            cmd = [
                'gcloud', 'compute', 'scp',
                f'{instance_name}:{remote_path}',
                str(local_file_path),
                f'--zone={zone}'
            ]
            
            print(f"Copying {remote_path} -> {local_file_path}")
            
            # Execute the copy command
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                success_count += 1
                print(f"  ✓ Success")
            else:
                error_count += 1
                print(f"  ✗ Error: {result.stderr.strip()}")
                
        except Exception as e:
            error_count += 1
            print(f"  ✗ Exception copying {remote_path}: {e}")
    
    print(f"\nTransfer complete: {success_count} successful, {error_count} errors")
    return error_count == 0

def main():
    if len(sys.argv) < 2:
        print("Usage: transfer_files.py <transfer_file> [instance_name] [zone] [local_dir]")
        print("\nExample:")
        print("  transfer_files.py TRANSFER")
        print("  transfer_files.py TRANSFER instance-20250817-211253 us-central1-c")
        print("  transfer_files.py TRANSFER instance-20250817-211253 us-central1-c results/")
        sys.exit(1)
    
    transfer_file = sys.argv[1]
    instance_name = sys.argv[2] if len(sys.argv) > 2 else "instance-20250817-211253"
    zone = sys.argv[3] if len(sys.argv) > 3 else "us-central1-c"
    local_dir = sys.argv[4] if len(sys.argv) > 4 else "genomes"
    
    print(f"Transfer configuration:")
    print(f"  Transfer file: {transfer_file}")
    print(f"  Instance: {instance_name}")
    print(f"  Zone: {zone}")
    print(f"  Local directory: {local_dir}")
    print()
    
    success = copy_files_from_instance(transfer_file, instance_name, zone, local_dir)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()