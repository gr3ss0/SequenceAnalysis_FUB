rule variability_per_position:
    input:
        msa="results/phylo_tree/msa.fa"
    output:
        variability="results/variability/variability_per_position.tsv"
    log:
        "logs/variability_per_position.log"
    threads: 1
    conda:
        "../envs/phylo.yaml"
    script:
        "../scripts/calculate_variability.py"

rule variability_per_window:
    input:
        variability = rules.variability_per_position.output.variability
    output:
        variability = f"results/variability/average_variability_{WINDOW_SIZE}.tsv"
    log:
        "logs/variability_per_window.log"
    threads: 1
    conda:
        "../envs/phylo.yaml"
    params:
        window_size = WINDOW_SIZE #this is defined in Snakefile, comes from config
    script:
        "../scripts/average_variability.py"

rule plot_variability:
    input:
        variability = rules.variability_per_window.output.variability
    output:
        variability_plot = f"results/variability/variability_plot_{WINDOW_SIZE}.png"
    log:
        "logs/plot_variability.log"
    threads: 1
    conda:
        "../envs/phylo.yaml"
    params:
        window_size = WINDOW_SIZE
    script:
        "../scripts/plot_variability.py"