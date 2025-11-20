rule cluster_genes:
    input:
        input_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        nu_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_Nu.OI.Avg.tsv"
    output:
        output_clusters="../outputs/gi_scores/{screen}/clusters/gene_clusters_{score}.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_gene_clustering.log"
    conda:
        "../envs/clustering-env.yaml"
    params:
        soft_threshold_power=lambda wildcards: config["SOFT_THRESHOLD_POWER"][config["PHENOTYPES_TO_CLUSTER"].index(wildcards.score)]
    script:
        "../scripts/cluster_genes.R"


rule cluster_all_phenotypes:
    input:
        lambda wildcards: _expand_gene_clusters(wildcards)
