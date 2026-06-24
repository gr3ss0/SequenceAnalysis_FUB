rule down_sample:
	input:
		unpack(get_map_input),
	output:
		r1 = "results/fastq/downsampled/{sample}_1.fastq.gz",
		r2 = "results/fastq/downsampled/{sample}_2.fastq.gz",
		hist = "results/fastq/downsampled/{sample}_histogram.txt",
	threads: 4
	log:
		"logs/downsample/{sample}.log"
	conda:
		"../envs/downsample.yaml"
	shell:
		"bbnorm.sh threads={threads} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} hist={output.hist} > {log} 2>&1"

rule assembly:
	input:
		r1 = "results/fastq/downsampled/{sample}_1.fastq.gz",
		r2 = "results/fastq/downsampled/{sample}_2.fastq.gz"
	output:
		"results/de-novo/assembly/{sample}/contigs.fasta"
	threads: 4
	log:
		"logs/assembly/{sample}.log"
	conda:
		"../envs/assembly.yaml"
	shell:
		"spades.py -t {threads} --phred-offset 33 -1 {input.r1} -2 {input.r2} -o results/de-novo/assembly/{wildcards.sample} > {log} 2>&1"

rule index_assembly:
	input:
		"results/de-novo/assembly/{sample}/contigs.fasta"
	output:
		multiext("results/de-novo/assembly/{sample}/contigs", 
				 ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
	params:
		output_file="results/de-novo/assembly/{sample}/contigs"
	threads: 4
	log:
		"logs/index_assembly/{sample}.log"
	conda:
		"../envs/mapping.yaml"
	shell:
		"bowtie2-build --threads {threads} -f {input} {params.output_file} > {log} 2>&1"


rule map_assembly:
	input:
		index_files = rules.index_assembly.output,
		r1 = "results/fastq/downsampled/{sample}_1.fastq.gz",
		r2 = "results/fastq/downsampled/{sample}_2.fastq.gz"
	output:
		bam = "results/de-novo/mapping/{sample}.bam"
	log:
		"logs/mapping/{sample}.log"
	threads: 4
	params:
		min_frag=config["bowtie2_params"]["min_fragment_ln"],
		seed=config["bowtie2_params"]["seed"],
		index_prefix="results/de-novo/assembly/{sample}/contigs",
	conda:
		"../envs/mapping.yaml"
	shell:
		"""
        bowtie2 --threads {threads} --minins {params.min_frag} --seed {params.seed} \
            -x {params.index_prefix} \
            -1 {input.r1} -2 {input.r2} 2>{log} \
            | samtools sort -@{threads} -o {output.bam} - >>{log} 2>&1
        """

rule index_bam:
    input:
        "results/de-novo/mapping/{sample}.bam"
    output:
        "results/de-novo/mapping/{sample}.bam.bai"
    conda:
        "../envs/mapping.yaml" # Or wherever samtools is installed
    shell:
        "samtools index {input}"


rule polish:
	input:
		assembly = "results/de-novo/assembly/{sample}/contigs.fasta",
		mapped = rules.map_assembly.output,
		mapped_idx = rules.index_bam.output,
	output:
		"results/de-novo/assembly/{sample}/polished_contigs.fasta"
	threads: 4
	log:
		"logs/polish/{sample}.log"
	conda:
		"../envs/assembly.yaml"
	params:
		output_prefix = "results/de-novo/assembly/{sample}/polished_contigs",
	shell:
		"pilon --genome {input.assembly} --frags {input.mapped} --output {params.output_prefix} > {log} 2>&1"


#https://snakemake-wrappers.readthedocs.io/en/v2.10.0/wrappers/ragtag/scaffold.html
rule scaffold:
	input:
		query="results/de-novo/assembly/{sample}/polished_contigs.fasta",
		ref=config["ref"],
		ref_index="results/index/reference.mmi",
	output:
		fasta="results/de-novo/assembly/{sample}/scaffolds.fasta",
		agp="results/de-novo/assembly/{sample}/ragtag.scaffold.agp",
		stats="results/de-novo/assembly/{sample}/ragtag.scaffold.stats",
	params:
		extra="",
	threads: 16
	log:
		"logs/ragtag/{sample}_scaffold.log",
	wrapper:
		"v2.10.0/bio/ragtag/scaffold"


rule concat_assembly:
	input:
		expand("results/de-novo/assembly/{sample}/scaffolds.fasta", sample=SAMPLES.index)
	output:
		"results/de-novo/assembly/combined.fa"
	log:
		"logs/concat_assembly.log"
	threads: 1
	shell:
		"cat {input} > {output}"



rule msa_de_novo:
	input:
		rules.concat_assembly.output
	output:
		"results/de-novo/phylo_tree/msa.fa"
	log:
		"logs/msa.log"
	threads: 20
	resources:
		mem_mb=8000  # Gives the job 8GB of RAM to prevent crashes
	conda:
		"../envs/phylo.yaml"
	shell:
		"mafft --thread {threads} --nomemsave --auto {input} > {output} 2> {log}"


rule phylo_tree_de_novo:
	input:
		rules.msa_de_novo.output
	output:
		multiext("results/de-novo/phylo_tree/iqtree_out", ".bionj", ".log", ".mldist", ".model.gz", ".treefile", ".iqtree", ".ckp.gz")
		
	log:
		"logs/phylo_tree.log"
	threads: 4
	conda:
		"../envs/phylo.yaml"
	shell:
		"""
		iqtree -s {input} -nt {threads} -pre results/de-novo/phylo_tree/iqtree_out > {log} 2>&1
		"""

rule visualize_tree_de_novo:
	input:
		"results/de-novo/phylo_tree/iqtree_out.treefile"
	output:
		"results/de-novo/phylo_tree/adenovirus_tree.png"
	conda:
		"../envs/phylo.yaml"
	script:
		"scripts/visualize_tree.py"