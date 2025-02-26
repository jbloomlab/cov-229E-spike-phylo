# Description:
# Main snakemake file that contains the rules, the output files, and 
# configures the order of rule execution.

# Author:
# Caleb Carr

# Imports
import pandas as pd 
import os
from snakemake.utils import min_version
min_version("7.25.4")

# Initialize config file global variable
configfile: "Configure/config.yml"


# Rules included from Rules/ directory 
include: "Rules/download_and_process_accessions.smk"
include: "Rules/extract_ORFs.smk"
include: "Rules/align_sequences.smk"
include: "Rules/augur_processing.smk"
include: "Rules/calculate_sequence_variation.smk"

GENES = [
    "spike",
]

# Output files
rule all:
    """ 
    This rule controls all the output files expected from the workflow. 
    """
    input:
        # ORF extraction log
        expand("Results/{gene}/ORF_sequences/ORF_extraction_verification.txt", gene=GENES),
        # Auspice
        expand("auspice/{gene}.json", gene=GENES),
        # Natural variation
        expand("Results/{gene}/Alignments/{gene}_natural_variation.csv", gene=GENES)


# Make DAG of rule execution
rule make_DAG:
    """
    This rule produces a DAG graph of rule execution order.
    """
    shell:
        """
        mkdir -p Results/DAG/
        snakemake --rulegraph | dot -Tpdf > Results/DAG/dag.pdf
        """
