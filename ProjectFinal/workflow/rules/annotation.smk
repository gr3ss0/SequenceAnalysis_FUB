rule prokka_prediction:
	input:
		assembly = "results/assembly/{sample}/secondary_assembly.fasta"
	output:
		# Prokka creates a GFF with FASTA appended automatically
		out_dir = directory("results/annotation/{sample}"),
		genes = "results/annotation/{sample}/{sample}.gff",
		proteins = "results/annotation/{sample}/{sample}.faa",
		stats = "results/annotation/{sample}/{sample}.txt"
	threads: 8
	conda:
		"../envs/annotation.yaml"
	log:
		"logs/annotation/{sample}.log"
	shell:
		"""
		prokka --cpus {threads} \
			   --outdir {output.out_dir} \
			   --prefix {wildcards.sample} \
			   --locustag {wildcards.sample} \
			   --force \
			   {input.assembly} > {log} 2>&1
		"""

rule prokka_predict_external:
	input:
		assembly = config["external_genome"]
	output:
		out_dir = directory("results/annotation_ext/external_genome"),
		genes = "results/annotation_ext/external_genome/external_genome.gff",
		proteins = "results/annotation_ext/external_genome/external_genome.faa",
		stats = "results/annotation_ext/external_genome/external_genome.txt"
	threads: 8
	conda:
		"../envs/annotation.yaml"
	log:
		"logs/annotation/external_genome.log"
	shell:
		"""
		prokka --cpus {threads} \
			   --outdir {output.out_dir} \
			   --prefix external_genome \
			   --locustag EXTERNAL \
			   --force \
			   {input.assembly} > {log} 2>&1
		"""
	
rule panaroo_core_genome:
	input:
		annotations = get_annotation_inputs
	output:
		core_aln = "results/annotation/core/core_gene_alignment.aln",
		results=directory("results/annotation/core")
	threads: 64
	params:
		mode="--alignment core",
		treshold="--core_threshold 0.95",
		sequence_identity="--threshold 0.98"
	log:
		"logs/annotation/core_genome.log"
	conda:
		"../envs/annotation.yaml"
	shell:
		"panaroo -i {input.annotations} -o {output.results} --threads {threads} {params.mode} {params.treshold} --clean-mode sensitive > {log} 2>&1"

rule iqtree_phylogeny:
	input:
		core_aln = rules.panaroo_core_genome.output.core_aln,
	output:
		multiext("results/phylo_tree/iqtree_out", ".bionj", ".log", ".mldist", ".model.gz", ".treefile", ".iqtree", ".ckp.gz")
	log:
		"logs/phylo_tree.log"
	threads: 4
	conda:
		"../envs/phylo.yaml"
	shell:
		"""
		iqtree -s {input} -nt {threads} -pre results/phylo_tree/iqtree_out > {log} 2>&1
		"""

rule core_genome_tree_plot:
	input:
		"results/phylo_tree/iqtree_out.treefile"
	output:
		"results/phylo_tree/core_genome_tree.png"
	log:
        "logs/prot_phylogeny/plot_tree.log"
	conda:
		"../envs/phylo.yaml"
	script:
		"../scripts/visualize_tree.py"

