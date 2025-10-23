# Description:
# Downloads accessions and process/extract sequences and metadata

# Author:
# Caleb Carr

# modified by:
# Sheri Harari

rule download_and_process_accessions:
    """
    This rule runs the download_NCBI_sequences.py script which
    downloads all genbank files from a list of accessions and 
    extracts metadata. 
    """
    input:
        accession_list = config["Accession_list"],
    params:
        genome_size_threshold_lower = config["Genome_size_threshold_lower"],
        genome_size_threshold_upper = config["Genome_size_threshold_upper"],
        max_frac_N = config["max_frac_N"]["genes"],
        accesstions_to_exclude = config["Accessions_to_exclude"],
        accessions_from_pre_study = config['Accessions_from_pre_study'],
        accessions_to_include = config["Accessions_to_include"],
        remove_duplicates = config["Remove_duplicates"],
        asia = config["Asia"],
        oceania = config["Oceania"],
        africa = config["Africa"],
        europe = config["Europe"],
        south_america = config["South America"],
        north_america = config["North America"],
        genbank_features = config["Genbank_features"],
        host_grouper = config["Host_grouper"],
        gene_names = lambda wildcards: config["Gene_groupings"][wildcards.gene],
        desired_gene = lambda wildcards: wildcards.gene,
    output:
        fasta_sequences = "Results/{gene}/nucleotide.fasta",
        metadata = "Results/{gene}/metadata.tsv",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/download_and_process_accessions.txt",
    script:
        "../Scripts/download_NCBI_sequences.py"

rule correct_metadata_dates:
    """
    This rule corrects collection dates in metadata based on a JSON corrections file.
    """
    input:
        metadata = "Results/{gene}/metadata.tsv",
        corrections = config["date_corrections"],
    output:
        corrected_metadata = "Results/{gene}/metadata_corrected.tsv",
    params:
        strain_col = "accession",  
        date_col = "date",      
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/correct_metadata_dates.txt",
    shell:
        "python Scripts/correct_metadata_dates.py "
        "--metadata {input.metadata} "
        "--corrections {input.corrections} "
        "--output {output.corrected_metadata} "
        "--strain-col {params.strain_col} "
        "--date-col {params.date_col} &> {log}"


rule create_new_references:
    """
    This rule runs a script to create new references.
    """
    input:
        config["Reference_genbank"],
    params:
        gene_name = lambda wildcards: config["Gene_groupings"][wildcards.gene][0],
    output:
        new_reference = "Results/{gene}/{gene}_reference.gb",
    conda:
        "../environment.yml",
    log:
        "Results/Logs/{gene}/create_new_references.txt",
    shell:
        "python Scripts/newreference.py "
        "--reference {input} "
        "--output-genbank {output.new_reference} "
        "--new-name spike_full_length "
        "--gene {params.gene_name} &> {log}"
