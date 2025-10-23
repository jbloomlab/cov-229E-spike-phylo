"""
Correct collection dates in metadata file based on a JSON corrections file.

Author: Sheri Harari
"""

import json
import pandas as pd
import argparse
from datetime import datetime

def format_date_with_unknown_precision(date_str):
    """
    Format dates with XX for unknown month/day
    - "2017" becomes "2017-XX-XX"
    - "2017-06" becomes "2017-06-XX"
    - "2017-06-15" stays "2017-06-15"
    """
    date_str = str(date_str).strip()
    
    # Count components to determine precision
    parts = date_str.split('-')
    
    if len(parts) == 1:  # Year only
        return f"{parts[0]}-XX-XX"
    elif len(parts) == 2:  # Year-month
        return f"{parts[0]}-{parts[1]}-XX"
    elif len(parts) == 3:  # Full date
        # Validate it's a proper date
        try:
            date_obj = pd.to_datetime(date_str)
            return date_obj.strftime('%Y-%m-%d')
        except:
            print(f"Warning: Could not parse date '{date_str}', keeping as-is")
            return date_str
    else:
        print(f"Warning: Unexpected date format '{date_str}', keeping as-is")
        return date_str

def correct_metadata_dates(metadata_file, corrections_file, output_file, 
                          strain_col='strain', date_col='date'):
    """
    Update dates in metadata TSV based on corrections in JSON file
    
    Parameters:
    -----------
    metadata_file : str
        Path to input metadata TSV file
    corrections_file : str
        Path to JSON file with date corrections (accession: date pairs)
    output_file : str
        Path to save corrected metadata
    strain_col : str
        Name of column containing strain/accession IDs
    date_col : str
        Name of column containing collection dates
    """
    
    # Load corrections
    print(f"Loading corrections from {corrections_file}...")
    with open(corrections_file, 'r') as f:
        date_corrections = json.load(f)
    print(f"Found {len(date_corrections)} accessions to correct\n")
    
    # Load metadata
    print(f"Loading metadata from {metadata_file}...")
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"Loaded {len(metadata)} rows\n")
    
    # Check if columns exist
    if strain_col not in metadata.columns:
        raise ValueError(f"Column '{strain_col}' not found in metadata. Available columns: {list(metadata.columns)}")
    if date_col not in metadata.columns:
        raise ValueError(f"Column '{date_col}' not found in metadata. Available columns: {list(metadata.columns)}")
    
    # Track updates
    updated_count = 0
    not_found = []
    
    # Update dates
    print("Updating dates...")
    for accession, corrected_date in date_corrections.items():
        mask = metadata[strain_col] == accession
        
        if mask.any():
            # Format date with XX for unknown precision
            formatted_date = format_date_with_unknown_precision(corrected_date)
            
            old_date = metadata.loc[mask, date_col].values[0]
            metadata.loc[mask, date_col] = formatted_date
            
            print(f"  ✓ {accession}: {old_date} → {formatted_date}")
            updated_count += 1
        else:
            not_found.append(accession)
            print(f"  ✗ {accession}: not found in metadata")
    
    # Report results
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Updated: {updated_count} dates")
    print(f"  Not found: {len(not_found)} accessions")
    
    if not_found and len(not_found) <= 10:
        print(f"  Missing accessions: {', '.join(not_found)}")
    
    # Save
    metadata.to_csv(output_file, sep='\t', index=False)
    print(f"\nCorrected metadata saved to: {output_file}")
    
    return metadata

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Correct collection dates in metadata based on JSON corrections file"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Input metadata TSV file"
    )
    parser.add_argument(
        "--corrections",
        required=True,
        help="JSON file with date corrections (accession: date pairs)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output corrected metadata TSV file"
    )
    parser.add_argument(
        "--strain-col",
        default="strain",
        help="Name of column containing strain/accession IDs (default: strain)"
    )
    parser.add_argument(
        "--date-col",
        default="date",
        help="Name of column containing dates (default: date)"
    )
    
    args = parser.parse_args()
    
    correct_metadata_dates(
        metadata_file=args.metadata,
        corrections_file=args.corrections,
        output_file=args.output,
        strain_col=args.strain_col,
        date_col=args.date_col
    )