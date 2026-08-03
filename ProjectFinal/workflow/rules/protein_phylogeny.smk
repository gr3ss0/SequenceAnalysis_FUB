DIAMOND_DB_FILE = config["diamond_db"]

rule diamond_build:
    input:
        protein_pool=rules.combine_prokka_proteins.output.protein_pool
    output:
        db_file=DIAMOND_DB_FILE + ".dmnd"
    conda:
        "../envs/diamond.yaml"
    log:
        "logs/prot_phylogeny/build.log"
    params:
        db=DIAMOND_DB_FILE
    threads: 16
    shell:
        """
        diamond makedb \
            --in {input.protein_pool} \
            -d {params.db} \
            > {log} 2>&1
        """


rule diamond_find_align_orthologs:
    input:
        db_instance = rules.diamond_build.output.db_file, 
        queries = QUERY_FILE
    output:
        alignment= "results/protein_queries/alignment.tsv"
    log:
        "logs/prot_phylogeny/blast_protein.log"
    conda:
        "../envs/diamond.yaml" 
    threads: 32
    params:
        format = "--outfmt 6 qseqid sseqid pident"
    shell:
        "diamond blastp -d {input.db_instance} -q {input.queries} -o {output.alignment} {params.format}"

rule extract_homologs:
    input:
        tsv = rules.diamond_find_align_orthologs.output.alignment,
        protein_pool=rules.combine_prokka_proteins.output.protein_pool,
        queries=QUERY_FILE
    output:
        fasta_files = expand("results/protein_queries/fasta/{query}.fasta", query=QUERY_PROTEINS),
    conda:
        "../envs/diamond.yaml"
    log:
        "logs/prot_phylogeny/extract_homologs.log"
    params:
        threshold = 50
    script:
        "../scripts/extract_homologs.py"

rule mafft_msa:
    input:
        homos = "results/protein_queries/fasta/{query}.fasta"
    output:
        msa =  "results/protein_queries/fasta/{query}.aln"
    log:
        "logs/prot_phylogeny/mafft_{query}.log"
    threads: 2
    conda:
        "../envs/diamond.yaml"
    shell:
        "mafft --thread {threads} --auto {input.homos} > {output.msa} 2> {log}"


rule tree_per_protein:
    input:
        alignment = rules.mafft_msa.output.msa,
    output:
        tree = "results/protein_queries/trees/{query}/{query}.treefile"
    log:
        "logs/prot_phylogeny/iqtree_{query}.log"
    threads: 4
    params:
        prefix=lambda wildcards: f"results/protein_queries/trees/{wildcards.query}/{wildcards.query}"
    conda:
        "../envs/phylo.yaml"
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
        rules.tree_per_protein.output.tree
    output:
        "results/protein_queries/trees/{query}/{query}.png"
    conda:
        "../envs/phylo.yaml"
    log:
        "logs/prot_phylogeny/visualize_tree_{query}.log"
    script:
        "../scripts/visualize_tree.py"

#rule run_prot_module:
    #input:
        #expand("results/protein_queries/trees/{query}/{query}.png", query=QUERY_PROTEINS),