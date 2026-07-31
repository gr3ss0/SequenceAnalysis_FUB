import glob

def get_db_const(wildcards):
    return {
        "name":config["diamond_db"] + "{sample}.dmnd",
        "prefix": config["diamond_db"] + "{sample}"
    }

DIAMOND_DB_FILE = config["diamond_db"] + ".dmnd"
DIAMOND_DB_PREFIX = DIAMOND_DB_FILE.replace(".dmnd", "") # Simple string manipulation
QUERY_FILE=config["protien_sequences"]

QUERY_PROTEINS = [line[1:].strip().split('|')[1] for line in open(QUERY_FILE) if line.startswith(">")]
#print(QUERY_PROTEINS)

## deprecated
# rule diamond_build:
#     input:
#         protein_pool="results/annotation/core/combined_protein_CDS.fasta",
#         panaroo_done=rules.panaroo_core_genome.output
#     output:
#         db=DIAMOND_DB_FILE,
#     conda:
#         "../envs/diamond.yaml"
#     params:
#         db_prefix = DIAMOND_DB_PREFIX
#     log:
#         "logs/prot_phylogeny/build.log"
#     threads: 16
#     shell:
#         "diamond makedb --in {input.protein_pool} -d {params.db_prefix} > {log} 2>&1"


rule diamond_build:
    input:
        protein_pool=rules.prokka_prediction.output.proteins,
    output:
        db=config["diamond_db"] + "{sample}.dmnd",
    conda:
        "../envs/diamond.yaml"
    params:
        db_prefix = config["diamond_db"] + "{sample}"
    log:
        "logs/prot_phylogeny/build_{sample}.log"
    threads: 16
    shell:
        "diamond makedb --in {input.protein_pool} -d {params.db_prefix} > {log} 2>&1"

rule diamond_find_align_orthologs:
    input:
        db_instance = rules.diamond_build.output.db,
        queries=QUERY_FILE,
    output:
        alignment= "results/protein_queries/{sample}_alignment.tsv"
    log:
        "logs/prot_phylogeny/diamond_blast_{sample}.log"
    conda:
        "../envs/diamond.yaml" # Environment containing blast/python
    threads: 32
    params:
        format = "--outfmt 6 qseqid sseqid pident",
        db_prefix = rules.diamond_build.params.db_prefix
    shell:
        "diamond blastp --threads {threads} -d {params.db_prefix} -q {input.queries} -o {output.alignment} {params.format} > {log} 2>&1"

rule ortholog_search_all:
    input:
        expand(rules.diamond_find_align_orthologs.output.alignment, sample=SAMPLES.index),

# this has to be a checkpoint as we cant be sure about the result of ortholog search - a query protein can be included or not.
checkpoint extract_homologs:
    input:
        tsvs = expand(rules.diamond_find_align_orthologs.output.alignment, sample=SAMPLES.index),
        protein_pools=expand(rules.diamond_build.input.protein_pool, sample=SAMPLES.index),
    output:
        splitted_dir=directory("results/protein_queries/fasta"),
    params:
        threshold = 50
    log:
        "logs/extract_homologs.log"
    conda:
        "../envs/diamond.yaml"
    script:
        "../scripts/extract_homologs.py"

def plot_trees_for_succ_queries(wildcards):
    ck_output = checkpoints.extract_homologs.get(**wildcards).output.splitted_dir
    queries = [
        os.path.splitext(f)[0] 
        for f in os.listdir(ck_output) 
        if f.endswith(".fasta")
    ]
    return expand("results/protein_queries/trees/{query}.png", query=queries)


rule mafft_msa:
    input:
        homos = "results/protein_queries/fasta/{query}.fasta",
        #conenction=plot_trees_for_succ_queries
    output:
        msa =  "results/protein_queries/aln/{query}.aln"
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
        multiext("results/protein_queries/trees/{query}", ".bionj", ".log", ".mldist", ".model.gz", ".treefile", ".iqtree", ".ckp.gz")
    log:
        "logs/prot_phylogeny/iqtree_{query}.log"
    threads: 4
    params:
        prefix = "results/protein_queries/trees/{query}"
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
        "results/protein_queries/trees/{query}.treefile"
    output:
        "results/protein_queries/trees/{query}.png"
    log:
        "logs/prot_phylogeny/plot_tree.log"
    conda:
        "../envs/phylo.yaml"
    script:
        "../scripts/visualize_tree.py"

rule run_prot_module:
    input:
        plot_trees_for_succ_queries,