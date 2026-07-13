def get_map_input_short(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return {
            "r1": SAMPLES_SHORT.at[wildcards.sample, 'fq1'], 
            "r2": SAMPLES_SHORT.at[wildcards.sample, 'fq2']
        }
    else:
        # Updated to point to fastp's compressed outputs
        return {
            "r1": f"results/trimmed/short/{wildcards.sample}.1.fastq.gz", 
            "r2": f"results/trimmed/short/{wildcards.sample}.2.fastq.gz"
        }

def get_map_input_long(wildcards):
    if config["analysis_options"].get("skip_trimming", True):
        return SAMPLES_LONG.at[wildcards.sample, 'fq']
        
    else:
        # Updated to point to fastp's compressed outputs
        return f"results/trimmed/long/{wildcards.sample}.fastq.gz"
            