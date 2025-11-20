library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(impute)

# Load data
nu_r1 <- read_tsv(snakemake@input[["input_nu_r1"]])
nu_r2 <- read_tsv(snakemake@input[["input_nu_r2"]])
nu_gene <- read_tsv(snakemake@input[["input_nu_gene_level"]])

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol_nu_r1 <- nu_r1[nu_r1$Category == "X+Y", ]
noncontrol_nu_r2 <- nu_r2[nu_r2$Category == "X+Y", ]

# Make a grid with gene names
gene_grid_r1 <- expand.grid(unique(noncontrol_nu_r1$Pseudogene1), unique(noncontrol_nu_r1$Pseudogene1))
gene_grid_r2 <- expand.grid(unique(noncontrol_nu_r2$Pseudogene1), unique(noncontrol_nu_r2$Pseudogene1))
# Merge in interaction scores
gene_grid_r1 <- left_join(gene_grid_r1, 
                          noncontrol_nu_r1, 
                          by = c("Var1"= "Pseudogene2", "Var2" = "Pseudogene1"))
gene_grid_r1 <- left_join(gene_grid_r1, 
                          noncontrol_nu_r1, 
                          by = c("Var1"= "Pseudogene1", "Var2" = "Pseudogene2"))
gene_grid_r2 <- left_join(gene_grid_r2, 
                          noncontrol_nu_r2, 
                          by = c("Var1" = "Pseudogene2", "Var2" = "Pseudogene1"))
gene_grid_r2 <- left_join(gene_grid_r2, 
                          noncontrol_nu_r2, 
                          by = c("Var1" = "Pseudogene1", "Var2" = "Pseudogene2"))
# Combine the columns
gene_grid_r1$Gene.GI <- ifelse(is.na(gene_grid_r1$InteractionScore.x), 
                               gene_grid_r1$InteractionScore.y, gene_grid_r1$InteractionScore.x)
gene_grid_r1$N <- ifelse(is.na(gene_grid_r1$N.x), 
                         gene_grid_r1$N.y, gene_grid_r1$N.x)
gene_grid_r2$Gene.GI <- ifelse(is.na(gene_grid_r2$InteractionScore.x), 
                               gene_grid_r2$InteractionScore.y, gene_grid_r2$InteractionScore.x)
gene_grid_r2$N <- ifelse(is.na(gene_grid_r2$N.x), 
                         gene_grid_r2$N.y, gene_grid_r2$N.x)


# Remove extraneous columns
gene_grid_r1 <- gene_grid_r1[,c("Var1", "Var2", "Gene.GI", "N")]
gene_grid_r1$Gene.GI <- as.numeric(gene_grid_r1$Gene.GI)
gene_grid_r1$N <- as.numeric(gene_grid_r1$N)
gene_grid_r2 <- gene_grid_r2[,c("Var1", "Var2", "Gene.GI", "N")]
gene_grid_r2$Gene.GI <- as.numeric(gene_grid_r2$Gene.GI)
gene_grid_r2$N <- as.numeric(gene_grid_r2$N)

# Go from long to wide
gene_mtx_r1 <- pivot_wider(gene_grid_r1[,1:3], names_from = Var2, values_from = Gene.GI)
gene_mtx_r1 <- as.data.frame(gene_mtx_r1)
gene_mtx_r2 <- pivot_wider(gene_grid_r2[,1:3], names_from = Var2, values_from = Gene.GI)
gene_mtx_r2 <- as.data.frame(gene_mtx_r2)
# Rename rows
rownames(gene_mtx_r1) <- gene_mtx_r1$Var1
rownames(gene_mtx_r2) <- gene_mtx_r2$Var1
# Remove redundant column
gene_mtx_r1 <- gene_mtx_r1[,2:ncol(gene_mtx_r1)]
gene_mtx_r2 <- gene_mtx_r2[,2:ncol(gene_mtx_r2)]

# Impute missing values
gene_mtx_r1 <- impute.knn(as.matrix(gene_mtx_r1))$data
gene_mtx_r2 <- impute.knn(as.matrix(gene_mtx_r2))$data
# Calculate correlations
gene_cor_r1 <- cor(as.matrix(gene_mtx_r1))
gene_cor_r2 <- cor(as.matrix(gene_mtx_r2))

# Wide to long
cors_r1 <- melt(gene_cor_r1)
cors_r2 <- melt(gene_cor_r2)
cors <- inner_join(cors_r1, cors_r2,
                   by = c("Var1", "Var2"),
                   suffix = c(".R1",".R2"))
cors2 <- cors[cors$Var1 != cors$Var2,]

# Remove duplicates
cors2$First <- apply(cors2[,c(1,2)], 1, min)
cors2$Second <- apply(cors2[,c(1,2)], 1, max)
cors2$PseudogeneCombinationName <- paste(cors2$First, cors2$Second, sep = ":")
#cors3 <- inner_join(cors2, id.map)
cors4 <- unique(cors2[,c("PseudogeneCombinationName", "value.R1", "value.R2")])

# Add interaction significance
cors4$Significance <- "NS"
cors4$Significance[cors4$PseudogeneCombinationName %in% 
                     nu_gene$PseudogeneCombinationName[nu_gene$Discriminant >
                                                       quantile(nu_gene$Discriminant[nu_gene$Category %in% 
                                                                                       c("X+NT", "NT+NT")], 1)[[1]]]] <- 
  "Significant"

# Rearrange for plotting
cors_tmp <- rbind(cors4[cors4$Significance == "NS",],
                  cors4[cors4$Significance == "Significant",])


cors_tmp$Significance <- factor(cors_tmp$Significance, levels = c("NS", "Significant"))

# Plot - FIGURE 2B
p2b <- ggplot(cors_tmp, aes(value.R1, value.R2, col = Significance)) +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .7, size = .25, pch = 16) +
  geom_point(data = cors_tmp[cors_tmp$Significance == "Significant",], alpha = .7, size = .3, pch = 16) +
  geom_point(data = cors_tmp[cors_tmp$Significance == "FA pathway",], alpha = 1, size = .5, pch = 16) +
  theme_bw() + 
  removeGrid() +
  theme(legend.position = "none",
        text = element_text(size = 10)) +
  xlab("Replicate 1") + ylab("Replicate 2") +
  scale_color_manual(values = c("#999999", "#BB5566")) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-.6, .6, .3)) + 
  scale_y_continuous(breaks = seq(-.6, .6, .3)) +
  annotate(geom = "text", x = .3, y = -.25, label = paste0("r = ", 
                                                           round(cor(cors4$value.R1, cors4$value.R2), 3))) +
  annotate(geom = "text", x = .3, y = -.3, 
            label = paste0("r = ", 
                           round(cor(cors4$value.R1[cors4$Significance == "Significant"], 
                           cors4$value.R2[cors4$Significance == "Significant"]), 3)))


# Save plots
ggsave(snakemake@output[["output_figure_2b"]], p2b, 
       device = "png", width = 3, height = 3, dpi = 300)
