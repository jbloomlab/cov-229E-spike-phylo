"""
Classify sequences into RBD classes based on defining mutations.

Author: Sheri Harari
"""

import pandas as pd
import argparse
from Bio import AlignIO
import sys

# RBD class defining mutations
LINEAGE_DEFINITIONS = {
    'class 2': ['318Y', '321R', '324V', '352V', '353A', '404L', '407I', '408N', '409S'],
    'class 3': ['316R', '321R', '324V', '310L', '314V', '352V', '353G', '357F', '350F', 
                '404L', '407L', '408N', '409S'],
    'class 4': ['310L', '311R', '314V', '316R', '318Y', '321R', '324V', '349Q', '350F', 
                '352V', '353G', '356K', '357F', '401M', '404L', '407L', '408N', '409S', '410H'],
    'class 5': ['309E', '310L', '310H', '312R', '314P', '318Y', '321R', '324V', '349Q', 
                '350F', '352V', '353G', '355V', '356K', '357F', '401M', '404L', '407L', 
                '408N', '410H', '411N'],
    'class 6': ['309E', '310L', '310H', '312R', '314P', '318F', '321R', '324V', '349Q', 
                '349R', '350F', '352V', '353G', '353D', '353S', '354N', '355V', '356K', 
                '357L', '358A', '358P', '358D', '401M', '404L', '405V', '406N', '406G', 
                '407H', '407Y', '410H', '411N']
}

def get_mutations(ref_seq, query_seq):
    """
    Get list of mutations between reference and query sequence.
    Returns mutations in format: 'A123B' (original AA, position, new AA)
    """
    mutations = []
    for pos, (ref_aa, query_aa) in enumerate(zip(ref_seq, query_seq), start=1):
        if ref_aa != query_aa and ref_aa != "-" and query_aa != "-":
            mutations.append(f"{ref_aa}{pos}{query_aa}")
    return mutations

def classify_rbd(mutations, strain_name):
    """
    Classify sequence into RBD class based on mutations.
    
    Parameters:
    -----------
    mutations : list
        List of mutations in format 'A123B'
    strain_name : str
        Name of strain (for logging)
    
    Returns:
    --------
    tuple: (rbd_class, matched_mutations_dict)
    """
    # Simplify mutations to position+newAA format (e.g., 'A310P' -> '310P')
    simplified_mutations = set()
    for mut in mutations:
        simplified = mut[1:]  # Remove original AA
        simplified_mutations.add(simplified)
    
    # Count matches for each lineage
    lineage_matches = {}
    for lineage, defining_muts in LINEAGE_DEFINITIONS.items():
        matched = list(set(defining_muts) & simplified_mutations)
        lineage_matches[lineage] = {
            'count': len(matched),
            'mutations': matched
        }
    
    # Check if all lineages have zero matches
    all_zero = all(v['count'] == 0 for v in lineage_matches.values())
    
    if all_zero:
        return 'unclassified', {}
    
    # Find lineage with most matches (prefer earlier class in ranking if tied)
    # Ranking order: class 2, class 3, class 4, class 5, class 6
    preferred_order = ['class 2', 'class 3', 'class 4', 'class 5', 'class 6']
    max_lineage = max(
        lineage_matches.keys(),
        key=lambda k: (lineage_matches[k]['count'], -preferred_order.index(k))
    )
    
    return max_lineage, lineage_matches[max_lineage]['mutations']

def add_rbd_classification(alignment_file, metadata_file, output_file, 
                          reference_strain, log_file=None):
    """
    Add RBD class column to metadata based on alignment.
    
    Parameters:
    -----------
    alignment_file : str
        Path to protein alignment FASTA file
    metadata_file : str
        Path to metadata TSV file
    output_file : str
        Path to save updated metadata
    reference_strain : str
        Name of reference strain in alignment
    log_file : str, optional
        Path to save classification details
    """
    
    # Read alignment
    print(f"Reading alignment from {alignment_file}...")
    alignment_dict = {}
    for record in AlignIO.read(alignment_file, "fasta"):
        alignment_dict[str(record.id)] = str(record.seq)
    
    print(f"Found {len(alignment_dict)} sequences in alignment")
    
    # Get reference sequence
    if reference_strain not in alignment_dict:
        raise ValueError(f"Reference strain '{reference_strain}' not found in alignment")
    
    ref_sequence = alignment_dict[reference_strain]
    print(f"Using reference strain: {reference_strain}")
    print(f"Reference strain will be automatically classified as class 1\n")
    
    # Read metadata
    print(f"Reading metadata from {metadata_file}...")
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"Found {len(metadata)} strains in metadata")
    
    # Check if 'strain' column exists
    if 'strain' not in metadata.columns:
        raise ValueError(f"'strain' column not found in metadata. Available columns: {list(metadata.columns)}")
    
    # Prepare logging
    log_lines = []
    log_lines.append("RBD Classification Results")
    log_lines.append("=" * 80)
    log_lines.append(f"\nReference strain: {reference_strain}")
    log_lines.append(f"Reference strain is automatically classified as: class 1")
    log_lines.append(f"\nDefining mutations for each class:")
    log_lines.append(f"  class 1: Reference strain (no defining mutations)")
    for lineage in ['class 2', 'class 3', 'class 4', 'class 5', 'class 6']:
        muts = LINEAGE_DEFINITIONS[lineage]
        log_lines.append(f"  {lineage}: {', '.join(muts)}")
    log_lines.append(f"\nClass ranking order: class 2 > class 3 > class 4 > class 5 > class 6")
    log_lines.append("\n" + "=" * 80)
    log_lines.append("\nClassification results:\n")
    
    # Classify each strain
    classifications = []
    class_counts = {'class 1': 0, 'class 2': 0, 'class 3': 0, 
                   'class 4': 0, 'class 5': 0, 'class 6': 0}
    
    for idx, row in metadata.iterrows():
        strain = row['strain']
        
        # Check if this is the reference strain
        if strain == reference_strain:
            classifications.append('class 1')
            class_counts['class 1'] += 1
            log_lines.append(f"{strain}: class 1 (reference strain - automatically assigned)")
            continue
        
        if strain not in alignment_dict:
            print(f"Warning: {strain} not found in alignment, assigning 'unclassified'")
            classifications.append('unclassified')
            log_lines.append(f"{strain}: unclassified (not in alignment)")
            continue
        
        # Get mutations
        query_seq = alignment_dict[strain]
        mutations = get_mutations(ref_sequence, query_seq)
        
        # Classify
        rbd_class, matched_muts = classify_rbd(mutations, strain)
        classifications.append(rbd_class)
        
        if rbd_class in class_counts:
            class_counts[rbd_class] += 1
        
        # Log details
        if matched_muts:
            log_lines.append(
                f"{strain}: {rbd_class} "
                f"(matched {len(matched_muts)} mutations: {', '.join(sorted(matched_muts))})"
            )
        else:
            log_lines.append(f"{strain}: {rbd_class} (no defining mutations matched)")
    
    # Add RBD class column to metadata
    metadata['rbd_class'] = classifications
    
    # Summary
    log_lines.append("\n" + "=" * 80)
    log_lines.append("\nSummary:")
    log_lines.append(f"  Total strains classified: {len(metadata)}")
    for rbd_class in ['class 1', 'class 2', 'class 3', 'class 4', 'class 5', 'class 6']:
        count = class_counts[rbd_class]
        log_lines.append(f"  {rbd_class}: {count} strains")
    
    unclassified_count = classifications.count('unclassified')
    if unclassified_count > 0:
        log_lines.append(f"  unclassified: {unclassified_count} strains")
    
    # Save metadata
    metadata.to_csv(output_file, sep='\t', index=False)
    print(f"\nUpdated metadata saved to: {output_file}")
    
    # Save log
    if log_file:
        with open(log_file, 'w') as f:
            f.write('\n'.join(log_lines))
        print(f"Classification details saved to: {log_file}")
    else:
        # Print to stdout if no log file specified
        print('\n'.join(log_lines))
    
    print(f"\nClassification complete!")
    print(f"  {sum(class_counts.values())} strains classified")
    for rbd_class in ['class 1', 'class 2', 'class 3', 'class 4', 'class 5', 'class 6']:
        count = class_counts.get(rbd_class, 0)
        if count > 0:
            print(f"  {rbd_class}: {count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify sequences into RBD classes based on defining mutations"
    )
    parser.add_argument(
        "--alignment",
        required=True,
        help="Input protein alignment FASTA file"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Input metadata TSV file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output metadata TSV file with rbd_class column"
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Name of reference strain in alignment"
    )
    parser.add_argument(
        "--log",
        help="Optional log file for classification details"
    )
    
    args = parser.parse_args()
    
    add_rbd_classification(
        alignment_file=args.alignment,
        metadata_file=args.metadata,
        output_file=args.output,
        reference_strain=args.reference,
        log_file=args.log
    )