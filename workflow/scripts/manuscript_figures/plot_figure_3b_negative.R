library(readr)
library(circlize)
library(dplyr)

# Import data
all_clusts <- read_tsv(snakemake@input[["input_clusters"]])
gamma_gene_gi <- read_tsv(snakemake@input[["input_gamma"]])
tau_gene_gi <- read_tsv(snakemake@input[["input_tau"]])
nu_gene_gi <- read_tsv(snakemake@input[["input_nu"]])

# Filter to negative nu
gamma_gene_filt <- gamma_gene_gi[gamma_gene_gi$PseudogeneCombinationID %in%
                                   nu_gene_gi$PseudogeneCombinationID[nu_gene_gi$Sig & 
                                                                        nu_gene_gi$InteractionScore < 0],] %>%
  filter(Category == "X+Y")


tau_gene_filt <- tau_gene_gi[tau_gene_gi$PseudogeneCombinationID %in%
                               nu_gene_gi$PseudogeneCombinationID[nu_gene_gi$Sig & 
                                                                    nu_gene_gi$InteractionScore < 0],] %>%
  filter(Category == "X+Y")


nu_gene_filt <- nu_gene_gi[nu_gene_gi$Sig & nu_gene_gi$InteractionScore < 0,] %>%
  filter(Category == "X+Y")


clust_names <- c("PARP1 interactors", 
                 "FA pathway",
                 "CST pathway",
                 "911 complex",
                 "CIA complex",
                 "RAD54L",
                 "EGFR signaling",
                 "BCDX2",
                 "COP9 signalosome",
                 "BRCA1-A complex",
                 "Protein phosphatase 2A",
                 "MRN complex")


gis <- expand.grid(all_clusts$gene[all_clusts$Cluster > 0], all_clusts$gene[all_clusts$Cluster > 0])
gis$First <- apply(gis[,1:2], 1, min)
gis$Second <- apply(gis[,1:2], 1, max)
gis$name <- paste(gis$First, gis$Second, sep = ":")
gis$Same <- 0

for(i in 1:nrow(gis)){
  j <- all_clusts$Cluster[all_clusts$gene == gis$First[i]]
  k <- all_clusts$Cluster[all_clusts$gene == gis$Second[i]]
  gis$Same[i] <- j == k
}


to_keep <- nu_gene_filt$PseudogeneCombinationName[nu_gene_filt$Sig & nu_gene_filt$InteractionScore < 0 & nu_gene_filt$Category == "X+Y"]
gis <- gis[gis$name %in% to_keep & gis$Same == 0,]
gis <- gis[!duplicated(gis$name),]
gis <- gis[!(gis$First %in% c("PARP1", "PARP2")) & !(gis$Second %in% c("PARP1", "PARP2")),]
gis$gamma <- ifelse(gis$name %in% gamma_gene_gi$PseudogeneCombinationName[gamma_gene_gi$Sig], 1, 0)
gis$tau <- ifelse(gis$name %in% tau_gene_gi$PseudogeneCombinationName[tau_gene_gi$Sig], 1, 0)
gis$nu <- ifelse(gis$name %in% nu_gene_gi$PseudogeneCombinationName[nu_gene_gi$Sig], 1, 0)


gis2 <- left_join(gis[,3:9], gamma_gene_gi[,c(3,6)], by = c("name" = "PseudogeneCombinationName"))
gis2 <- left_join(gis2, tau_gene_gi[,c(3,6)], c("name" = "PseudogeneCombinationName"))
gis2 <- left_join(gis2, nu_gene_gi[,c(3,6)], c("name" = "PseudogeneCombinationName"))


colnames(gis2)[1:2] <- c("from", "to")
grouping <- structure(all_clusts$Cluster[all_clusts$gene %in% unique(c(gis2$from, gis2$to))], names = all_clusts$gene[all_clusts$gene %in% unique(c(gis2$from, gis2$to))])
gene_info <- all_clusts[all_clusts$gene %in% names(grouping),]
gene_info <- gene_info[order(gene_info$Cluster),]
# gene.info$color <- "gray20"
# gene.info$color[gene.info$Cluster %% 2 == 0] <- "gray80"
gene_info$color <- sapply(gene_info$Cluster, 
                          function(x) c("#E5E5E5", "#F06FAA", "#4D2D89",
                                        "#9673B3", "#376DB5", "#70BF44", 
                                        "#BA2C32", "#96D2B0", "#924C21",
                                        "#DA6F27", "#009292", "#7DB2E0")[as.numeric(x)])


gis2$gammacolor <- "#FFFFFF00"
gis2$gammacolor[gis2$gamma == 1 & gis2$InteractionScore.x > 0] <- "#ff7f0030"
gis2$gammacolor[gis2$gamma == 1 & gis2$InteractionScore.x < 0] <- "#1e90ff30"
gis2$tau_color <- "#FFFFFF00"
gis2$nu_color <- "#FFFFFF00"
gis2$tau_color[gis2$tau == 1 & gis2$InteractionScore.y > 0] <- "#ff7f0030"
gis2$tau_color[gis2$tau == 1 & gis2$InteractionScore.y < 0] <- "#1e90ff30"
gis2$nu_color[gis2$nu == 1 & gis2$InteractionScore < 0] <- "#33716B30"


coloring <- structure(gene_info$color, names = gene_info$gene)


pdf(snakemake@output[["output_figure_3b_gamma_neg"]], width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$gammacolor, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.05), 
             preAllocateTracks = list(track.height = .05))
dev.off()


pdf(snakemake@output[["output_figure_3b_tau_neg"]], width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$tau_color, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.05), 
             preAllocateTracks = list(track.height = .05))
dev.off()


pdf(snakemake@output[["output_figure_3b_nu_neg"]], width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$nu_color, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.05), 
             preAllocateTracks = list(track.height = .05))
dev.off()
