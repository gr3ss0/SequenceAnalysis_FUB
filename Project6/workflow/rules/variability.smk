rule plot_variable_genes_heatmap:
    input:
        counts = "results/featureCounts/all.featureCounts"
    output:
        heatmap = "results/plots/top_variable_genes_heatmap.pdf"
    params:
        top_n = config["top_n_variable"]
    log:
        "logs/plots/heatmap.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/plot_heatmap.R"