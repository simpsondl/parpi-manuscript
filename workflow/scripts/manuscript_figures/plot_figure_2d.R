library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(impute)
library(colorBlindness)
library(ComplexHeatmap)
library(seriation)
library(circlize)

# Load data
nu_gene_gis <- read_tsv(snakemake@input[["input_nu"]])
clusters <- read_tsv(snakemake@input[["input_clusters"]])

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol_nu_gene_gis <- nu_gene_gis[!grepl("NTPG_", nu_gene_gis$PseudogeneCombinationName),]

# Make a grid with gene names
gene_grid <- expand.grid(unique(noncontrol_nu_gene_gis$Pseudogene1), unique(noncontrol_nu_gene_gis$Pseudogene1))
# Merge in interaction scores
gene_grid <- left_join(gene_grid, 
                       noncontrol_nu_gene_gis, 
                       by = c("Var1" = "Pseudogene2", "Var2" = "Pseudogene1"))
gene_grid <- left_join(gene_grid, 
                       noncontrol_nu_gene_gis, 
                       by = c("Var1" = "Pseudogene1", "Var2" = "Pseudogene2"))
# Combine the columns

gene_grid$Gene.GI <- ifelse(is.na(gene_grid$InteractionScore.x), gene_grid$InteractionScore.y, gene_grid$InteractionScore.x)
gene_grid$N <- ifelse(is.na(gene_grid$N.x), gene_grid$N.y, gene_grid$N.x)


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

# Impute missing values
gene_mtx <- impute.knn(as.matrix(gene_mtx))$data

# Get genes in clusters
to_keep <- unique(c(clusters$gene[clusters$cluster > 0]))

tmp_mtx <- gene_mtx[to_keep,to_keep]
d <- as.dist(1 - cor(as.matrix(tmp_mtx)))
# Rearrange using OLO algorithm from seriation
# Rearrange using OLO algorithm from seriation
o1 <- seriate(d, method = "OLO_average")
maxmag <- 3.5
# Define color spectrum ranges
col_fun <- colorRamp2(c(-1 * maxmag, 0, maxmag), c("#33716B", "white", "#D81B60"))

clust <- data.frame(gene = rownames(tmp_mtx))
clust <- left_join(clust, clusters[c("gene", "cluster")])

clust$cluster[clust$cluster == 0] <- NA
clust$gene <- factor(clust$gene, levels = to_keep)
clust <- clust[order(clust$gene),]
clust$cluster <- factor(as.character(clust$cluster), 
                        levels = as.character(1:12))

#row
c_both <- rowAnnotation(NuClust = clust$cluster,
                        show_annotation_name = FALSE, 
                        col = list(NuClust = setNames(c("#E5E5E5", "#F06FAA", "#4D2D89",
                                                        "#9673B3", "#376DB5", "#70BF44", 
                                                        "#BA2C32", "#96D2B0", "#924C21",
                                                        "#DA6F27", "#009292", "#7DB2E0"), as.character(1:12))), 
                        na_col = "white")


c_both_col <- columnAnnotation(NuClust = clust$cluster,
                               show_annotation_name = FALSE, 
                               show_legend = FALSE,
                               col = list(NuClust = setNames(c("#E5E5E5", "#F06FAA", "#4D2D89",
                                                               "#9673B3", "#376DB5", "#70BF44", 
                                                               "#BA2C32", "#96D2B0", "#924C21",
                                                               "#DA6F27", "#009292", "#7DB2E0"), as.character(1:12))), 
                               na_col = "white")

hm <- Heatmap(t(tmp_mtx), 
              rect_gp = gpar(type = "none"), 
              cluster_rows = as.dendrogram(o1[[1]]), 
              cluster_columns = as.dendrogram(o1[[1]]), 
              show_row_dend = TRUE,
              show_column_dend = FALSE,
              row_dend_side = "left",
              heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i == j){
                  grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                } else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
              },
              col = col_fun,
              name = "GI",
              left_annotation = c_both, 
              bottom_annotation = c_both_col,
              show_column_names = TRUE,
              show_row_names = TRUE,
              column_names_side = "bottom",
              row_names_gp = gpar(fontsize = 4),
              column_names_gp = gpar(fontsize = 4),
              column_dend_height = unit(2, "cm"))

pdf(snakemake@output[["output_figure_2d"]], width = 10, height = 10)
draw(hm)
dev.off()
