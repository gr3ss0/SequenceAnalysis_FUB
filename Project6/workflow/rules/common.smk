def get_map_input(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return {
            "fq1": SAMPLES.at[wildcards.sample, 'fq1'], 
            "fq2": SAMPLES.at[wildcards.sample, 'fq2']
        }
    else:
        # Updated to point to fastp's compressed outputs
        return {
            "fq1": f"results/trimmed/{wildcards.sample}.1.fastq.gz", 
            "fq2": f"results/trimmed/{wildcards.sample}.2.fastq.gz"
        }