library(readr)
library(dplyr)
library(WGCNA)
library(impute)
library(tidyr)

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
    source("scripts/dual_logging.R")
    .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
    on.exit({ 
        .dual_cleanup() 
    }, add = TRUE)
}

message(sprintf("[%s] cluster_genes.R starting", Sys.time()))

scores <- read_tsv(snakemake@input[["input_scores"]])
nu_scores <- read_tsv(snakemake@input[["nu_scores"]])
soft_threshold_power <- snakemake@params[["soft_threshold_power"]]

message(sprintf("[%s] Input gene scores: %s (%d rows)", Sys.time(),
                snakemake@input[["input_scores"]], nrow(scores)))
message(sprintf("[%s] Soft-thresholding power: %s",
                Sys.time(), soft_threshold_power))

# Remove interactions involving non-targeting guides
noncontrol_gene_gis <- scores[!grepl("NTPG_", scores$PseudogeneCombinationName),]
message(sprintf("[%s] Filtered non-targeting interactions: %d -> %d rows",
                Sys.time(), nrow(scores), nrow(noncontrol_gene_gis)))
message(sprintf("[%s] Unique genes to cluster: %d",
                Sys.time(), length(unique(c(noncontrol_gene_gis$Pseudogene1, noncontrol_gene_gis$Pseudogene2)))))

# Remove any constructs without nu scores
noncontrol_gene_gis <- noncontrol_gene_gis %>%
    filter(PseudogeneCombinationID %in% nu_scores$PseudogeneCombinationID)

message(sprintf("[%s] Filtered to be consistent with manuscript. Removed interactions without nu scores. Unique genes to cluster: %d",
                Sys.time(), length(unique(c(noncontrol_gene_gis$Pseudogene1, noncontrol_gene_gis$Pseudogene2)))))

# Make a grid with gene names
gene_grid <- expand.grid(unique(noncontrol_gene_gis$Pseudogene1), unique(noncontrol_gene_gis$Pseudogene1))
message(sprintf("[%s] Created gene grid: %d x %d -> %d rows",
                Sys.time(), length(unique(noncontrol_gene_gis$Pseudogene1)), 
                length(unique(noncontrol_gene_gis$Pseudogene1)), nrow(gene_grid)))

# Merge in interaction scores
gene_grid <- left_join(gene_grid, 
                       noncontrol_gene_gis, 
                       by = c("Var1" = "Pseudogene2", "Var2" = "Pseudogene1"))
gene_grid <- left_join(gene_grid, 
                       noncontrol_gene_gis, 
                       by = c("Var1" = "Pseudogene1", "Var2" = "Pseudogene2"))
# Combine the columns
gene_grid$Gene.GI <- ifelse(is.na(gene_grid$InteractionScore.x), gene_grid$InteractionScore.y, gene_grid$InteractionScore.x)
gene_grid$N <- ifelse(is.na(gene_grid$N.x), gene_grid$N.y, gene_grid$N.x)

na_before <- sum(is.na(gene_grid$Gene.GI))
message(sprintf("[%s] After merging scores to grid: %d NA interaction scores (will be imputed)",
                Sys.time(), na_before))

# Remove extraneous columns
gene_grid <- gene_grid[,c("Var1", "Var2", "Gene.GI", "N")]
gene_grid$Gene.GI <- as.numeric(gene_grid$Gene.GI)
gene_grid$N <- as.numeric(gene_grid$N)

# Go from long to wide
gene_mtx <- pivot_wider(gene_grid[,1:3], names_from = Var2, values_from = Gene.GI)
gene_mtx <- as.data.frame(gene_mtx)
### Rename rows
rownames(gene_mtx) <- gene_mtx$Var1
### Remove redundant column
gene_mtx <- gene_mtx[,2:ncol(gene_mtx)]

na_cells <- sum(is.na(as.matrix(gene_mtx)))
message(sprintf("[%s] Converted to wide matrix: %d genes x %d genes, %d missing values",
                Sys.time(), nrow(gene_mtx), ncol(gene_mtx), na_cells))

# Impute missing values
gene_mtx <- impute.knn(as.matrix(gene_mtx))$data

message(sprintf("[%s] Completed imputation: matrix dimensions %d x %d",
                Sys.time(), nrow(gene_mtx), ncol(gene_mtx)))

# Perform clustering
message(sprintf("[%s] Starting blockwiseModules (soft_threshold_power=%s)",
                Sys.time(), soft_threshold_power))
wg <- blockwiseModules(gene_mtx, 
                        power = soft_threshold_power, 
                        TOMType = "signed",
                        minModuleSize = 3,
                        numericLabels = TRUE)
message(sprintf("[%s] blockwiseModules completed", Sys.time()))

n_modules <- length(unique(wg$colors[!is.na(wg$colors)]))
message(sprintf("[%s] Number of assigned modules (including 0): %d", Sys.time(), n_modules))

# Prepare output table
assigned_annotations <- data.frame(gene = rownames(gene_mtx),
                                   cluster = wg$colors)

outpath <- snakemake@output[["output_clusters"]]
message(sprintf("[%s] Writing cluster assignments to: %s", Sys.time(), outpath))
write_tsv(assigned_annotations, outpath)
message(sprintf("[%s] cluster_genes.R completed", Sys.time()))
