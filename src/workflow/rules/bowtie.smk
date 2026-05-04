rule index_ref:
	input:
		ref = config["ref"]
	output:
		temp(multiext(REF_BASENAME, 
                 ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"))
	log:
		"logs/index_ref/"+REF_BASENAME+".log"
	threads: 4
	
	conda:
		"../envs/mapping.yaml"
	shell:
		"bowtie2-build --threads {threads} -f {input.ref} " + REF_BASENAME + " > {log} 2>&1" #<reference_in> <bt2_base>

rule map:
	input:
		r1 = lambda wildcards:SAMPLES.at[wildcards.sample,'fq1'],
		r2 = lambda wildcards:SAMPLES.at[wildcards.sample,'fq2'],
		ref_index = multiext(REF_BASENAME, 
		".1.bt2", ".2.bt2", ".3.bt2",          ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
	output:
		"results/sam/{sample}.sam"
	log:
		"logs/map/{sample}.log"
	threads: 4
	params:
		min_frag=config["bowtie2_params"]["min_fragment_ln"],
		seed=config["bowtie2_params"]["seed"]
	conda:
		"../envs/mapping.yaml"
	shell:
		"bowtie2 --threads {threads} {params.min_frag} {params.seed} -x " + REF_BASENAME + " -1 {input.r1} -2 {input.r2} -S {output} 2>{log}"
