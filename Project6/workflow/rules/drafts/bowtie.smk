rule index_ref:
	input:
		ref = config["ref"]
	output:
		multiext(config["bowtie2_params"]["ref_index_dir"] + REF_BASENAME, 
                 ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
	log:
		"logs/index_ref/"+REF_BASENAME+".log"
	threads: 4
	params:
		output_file=config["bowtie2_params"]["ref_index_dir"]+REF_BASENAME

	conda:
		"../envs/mapping.yaml"
	shell:
		"bowtie2-build --threads {threads} -f {input.ref} {params.output_file} > {log} 2>&1" #<reference_in> <bt2_base>


def get_map_input(wildcards):
    # Case 1: Trimming was skipped, use raw data from SAMPLES dataframe
    if config["analysis_options"].get("skip_trimming", True):
        return {
            "r1": SAMPLES.at[wildcards.sample, 'fq1'],
            "r2": SAMPLES.at[wildcards.sample, 'fq2']
        }
    
    # Case 2: Trimming was performed, use results from your trimming rule
    # Note: Use the actual filenames produced by your Trimmomatic rule
    else:
        return {
            "r1": f"results/trimmed/{wildcards.sample}_R1_paired.fq.gz",
            "r2": f"results/trimmed/{wildcards.sample}_R2_paired.fq.gz"
        }


rule map:
	input:
		unpack(get_map_input),
		index_files = multiext(config["bowtie2_params"]["ref_index_dir"]+REF_BASENAME, 
		".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
	output:
		"results/sam/{sample}.sam"
	log:
		"logs/map/{sample}.log"
	threads: 4
	params:
		min_frag=config["bowtie2_params"]["min_fragment_ln"],
		seed=config["bowtie2_params"]["seed"],
		index_basename=config["bowtie2_params"]["ref_index_dir"]+REF_BASENAME
	conda:
		"../envs/mapping.yaml"
	shell:
		"bowtie2 --threads {threads} --minins {params.min_frag} --seed {params.seed} -x {params.index_basename} -1 {input.r1} -2 {input.r2} -S {output} 2>{log}"
