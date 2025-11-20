rule plot_figure_1d:
    input:
        # use the consolidated construct scores file and subset inside the plotting script
        input_gamma="../outputs/gi_scores/screen2023/construct_scores/all_gis_Gamma.OI.Avg.tsv"
    output:
        output_figure_1d="../outputs/manuscript_figures/figure_1d.svg"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_1d.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_1d.R"


rule plot_figure_1e:
    input:
        # use the consolidated construct scores file and subset inside the plotting script
        input_tau="../outputs/gi_scores/screen2023/construct_scores/all_gis_Tau.OI.Avg.tsv"
    output:
        output_figure_1e="../outputs/manuscript_figures/figure_1e.svg"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_1e.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_1e.R"


rule plot_figure_2a:
    input:
        # updated paths: differential outputs are now placed with Nu.* naming under construct/gene/discriminant dirs
        input_nu_r1="../outputs/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R1.tsv",
        input_nu_r2="../outputs/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R2.tsv",
        input_nu_gene_level="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv"
    output:
        output_figure_2a="../outputs/manuscript_figures/figure_2a.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_2a.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_2a.R"


rule plot_figure_2b:
    input:
        input_nu_r1="../outputs/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R1.tsv",
        input_nu_r2="../outputs/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R2.tsv",
        input_nu_gene_level="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv"
    output:
        output_figure_2b="../outputs/manuscript_figures/figure_2b.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_2b.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_2b.R"


rule plot_figure_2d:
    input:
        input_nu="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_clusters="../outputs/gi_scores/screen2023/clusters/gene_clusters_Nu.OI.Avg.tsv"
    output:
        output_figure_2d="../outputs/manuscript_figures/figure_2d.pdf"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_2d.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_2d.R"


rule plot_figure_3a:
    input:
        input_gamma="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Gamma.OI.Avg.tsv",
        input_tau="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Tau.OI.Avg.tsv",
        input_nu="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv"
    output:
        output_figure_3a="../outputs/manuscript_figures/figure_3a.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_3a.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_3a.R"


rule plot_figure_3b_negative:
    input:
        input_gamma="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Gamma.OI.Avg.tsv",
        input_tau="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Tau.OI.Avg.tsv",
        input_nu="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_clusters="../outputs/gi_scores/screen2023/clusters/gene_clusters_Nu.OI.Avg.tsv"
    output:
        output_figure_3b_gamma_neg="../outputs/manuscript_figures/figure_3b_gamma_negative.png",
        output_figure_3b_tau_neg="../outputs/manuscript_figures/figure_3b_tau_negative.png",
        output_figure_3b_nu_neg="../outputs/manuscript_figures/figure_3b_nu_negative.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_3b_negative.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_3b_negative.R"


rule plot_figure_3b_positive:
    input:
        input_gamma="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Gamma.OI.Avg.tsv",
        input_tau="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Tau.OI.Avg.tsv",
        input_nu="../outputs/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_clusters="../outputs/gi_scores/screen2023/clusters/gene_clusters_Nu.OI.Avg.tsv"
    output:
        output_figure_3b_gamma_pos="../outputs/manuscript_figures/figure_3b_gamma_positive.png",
        output_figure_3b_tau_pos="../outputs/manuscript_figures/figure_3b_tau_positive.png",
        output_figure_3b_nu_pos="../outputs/manuscript_figures/figure_3b_nu_positive.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_3b_positive.log"
    conda:
        "../envs/manuscript-env.yaml"
    script:
        "../scripts/manuscript_figures/plot_figure_3b_positive.R"
