import glob

# # Check if the user defined a directory of query proteins in the config
# PROTEIN_DIR = config.get("optional_protein_dir", "")

# if PROTEIN_DIR and os.path.exists(PROTEIN_DIR):
#     # Capture all protein names (e.g., "pi")
#     QUERY_PROTEINS = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(f"{PROTEIN_DIR}/*.fasta")]
# else:
#     QUERY_PROTEINS = []

rule find_orthologs:
    input:
        #query = os.path.join(config["optional_protein_dir"], "{protein}.fasta"),
        # We need the predicted proteins from all your samples
        sample_proteins = expand("results/annotation/{sample}/{sample}.faa", sample=SAMPLES_LONG.index)
    output:
        combined_fasta = "results/targeted_phylogeny/{protein}/unaligned_hits.faa"
    # log:
    #     "logs/targeted_phylogeny/blast_{protein}.log"
    # conda:
    #     "../envs/blast.yaml" # Environment containing blast/python
    # run:
    #     import subprocess
    #     from bioservices import BLAST # or use standard local blastp commands via python
        
    #     # Step 1: Create a local BLAST database of all sample proteins combined
    #     # Step 2: BLAST the single query against it
    #     # Step 3: Parse the top hit per sample and write them plus the query to combined_fasta
    #     # (Alternatively, you can loop a local blastp command over your sample files)
        
    #     with open(output.combined_fasta, "w") as out_f:
    #         # First, write the original query protein pi into the file
    #         with open(input.query) as q_f:
    #             out_f.write(q_f.read() + "\n")
                
    #         # Loop through each sample, blast, and grab the best hit
    #         for sample_faa in input.sample_proteins:
    #             # Local command line call to blastp to grab the top hit
    #             # Append that sequence text straight to output.combined_fasta
    #             pass

rule align_protein:
    input:
        fasta = "results/targeted_phylogeny/{protein}/unaligned_hits.faa"
    output:
        alignment = "results/targeted_phylogeny/{protein}/alignment.aln"
    # log:
    #     "logs/targeted_phylogeny/mafft_{protein}.log"
    # threads: 2
    # conda:
    #     "../envs/mafft.yaml"
    # shell:
    #     "mafft --thread {threads} --auto {input.fasta} > {output.alignment} 2> {log}"

rule tree_per_protein:
    input:
        alignment = "results/targeted_phylogeny/{protein}/alignment.aln"
    output:
        tree = "results/targeted_phylogeny/{protein}/gene_tree.treefile"
    # log:
    #     "logs/targeted_phylogeny/iqtree_{protein}.log"
    # threads: 4
    # params:
    #     prefix = "results/targeted_phylogeny/{protein}/gene_tree"
    # conda:
    #     "../envs/iqtree.yaml"
    # shell:
    #     """
    #     iqtree \
    #         -s {input.alignment} \
    #         -pre {params.prefix} \
    #         -m MFP \
    #         -nt {threads} \
    #         -bb 1000 \
    #         > {log} 2>&1
    #     """