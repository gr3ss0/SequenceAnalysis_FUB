import glob
DIAMOND_DB_FILE = config["diamond_db"]
QUERY_FILE=config["protien_sequences"]

QUERY_PROTEINS = [line[1:].strip().split()[0] for line in open(QUERY_FILE) if line.startswith(">")]
rule diamond_build:
    input:
        protein_pool="results/annotation/core/combined_protein_CDS.fasta"
    output:
        db="DIAMOND_DB_FILE"+".dmb",
    conda:
        "../envs/diamond.yaml"
    log:
        "logs/prot_phylogeny/build.log"
    threads: 16
    shell:
        "diamond makedb --in {input.protein_pool} -d {output.db_file} > {log} 2>&1"


rule diamond_find_align_orthologs:
    input:
        db_instance = DIAMOND_DB_FILE,
        POI="results/protein_queries/{protein}.faa",
        sample_proteins = expand("results/annotation/{protein}/{protein}.faa", sample=SAMPLES_LONG.index),
    output:
        alignment= "results/protein_queries/alignment.tsv"
    log:
        "logs/prot_phylogeny/blast_{protein}.log"
    conda:
        "../envs/diamond.yaml" # Environment containing blast/python
    threads: 32
    params:
        format = "--outfmt 6 qseqid sseqid pident",
        queries=QUERY_FILE
    shell:
        "diamond blastp -d {input.db_instance} -q {params.queries} -o {output.alignment} {params.format}"

checkpoint extract_homologs:
    input:
        tsv = rules.diamond_find_align_orthologs.output.alignment,
        protein_pool=rules.diamond_build.input.protein_pool,
    output:
        splitted=directory("results/protein_queries/fasta"),
    params:
        threshold = "50"
    run:
        "../scripts/extract_homologs.py"

rule mafft_msa:
    input:
        homos = "results/protein_queries/fasta/{query}.fasta"
    output:
        msa =  "results/protein_queries/fasta/{query}.aln"
    log:
        "logs/prot_phylogeny/mafft_{protein}.log"
    threads: 2
    conda:
        "../envs/diamond.yaml"
    shell:
        "mafft --thread {threads} --auto {input.homos} > {output.msa} 2> {log}"


rule tree_per_protein:
    input:
        alignment = rules.mafft_msa.output.msa,
    output:
        out_dir = directory("results/protein_queries/trees/{query}")
    log:
        "logs/prot_phylogeny/iqtree_{query}.log"
    threads: 4
    params:
        prefix = "{output.out_dir}/{query}"
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        iqtree \
            -s {input.alignment} \
            -pre {params.prefix} \
            -nt {threads} \
            > {log} 2>&1
        """

rule visualize_tree:
    input:
        "results/protein_queries/trees/{query}/{query}.treefile"
    output:
        "results/protein_queries/trees/{query}/{query}.png"
    conda:
        "../envs/phylo.yaml"
    script:
        "../scripts/visualize_tree.py"

rule run_prot_module:
    input:
        expand("results/protein_queries/trees/{query}/{query}.png", query=QUERY_PROTEINS),