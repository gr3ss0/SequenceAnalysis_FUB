rule phylo_tree:
	input:
		"results/phylo_tree/msa.fa"
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

rule visualize_tree:
    input:
        "results/phylo_tree/iqtree_out.treefile"
    output:
        "results/phylo_tree/adenovirus_tree.png"
    conda:
        "../envs/phylo.yaml"
    script:
        "../scripts/visualize_tree.py"