library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)

nu_gene_gi_r1 <- read_tsv(snakemake@input[["input_nu_r1"]])
nu_gene_gi_r2 <- read_tsv(snakemake@input[["input_nu_r2"]])
nu_sig <- read_tsv(snakemake@input[["input_nu_gene_level"]])

# Merge replicate datasets
gc_nu <- full_join(nu_gene_gi_r1[, c("PseudogeneCombinationID", "Category", "InteractionScore")], 
                   nu_gene_gi_r2[, c("PseudogeneCombinationID", "Category", "InteractionScore")],
                   by = c("PseudogeneCombinationID", "Category"),
                   suffix = c(".R1", ".R2"))


gc_nu$Sig <- gc_nu$PseudogeneCombinationID %in% nu_sig$PseudogeneCombinationID[nu_sig$Hit]

gc_nu$CatSig <- gc_nu$Category
gc_nu$CatSig[gc_nu$Sig & gc_nu$Category != "X+X"] <- gc_nu$Sig[gc_nu$Sig & gc_nu$Category != "X+X"]

# Rearrange data for plotting
tmp_df <- rbind(gc_nu[gc_nu$Category == "X+Y", ],
                gc_nu[gc_nu$Category == "X+NT", ],
                gc_nu[gc_nu$Category == "NT+NT", ],
                gc_nu[gc_nu$CatSig == TRUE, ])

tmp_df$CatSig <- factor(tmp_df$CatSig, levels = c("X+Y", "X+NT", "NT+NT", "TRUE"))


# Plot - FIGURE 2A
p2a <- ggplot(tmp_df, aes(InteractionScore.R1, InteractionScore.R2, color = CatSig)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size = 1, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#999999", "#abddde", "#046c9a", "#BB5566")) +
  annotate("text", hjust = 0, x = -6, y = 10, 
           label = paste("r[a] ==", round(cor(tmp_df$InteractionScore.R1, 
                                              tmp_df$InteractionScore.R2), 3)), parse = TRUE) +
  annotate("text", hjust = 0, x = -6, y = 9.2, 
           label = paste("r[s] ==", round(cor(tmp_df$InteractionScore.R1[tmp_df$Sig], 
                                              tmp_df$InteractionScore.R2[tmp_df$Sig]), 3)), parse = TRUE) +
  annotate("text", hjust = 0, x = -6, y = 8.4, 
           label = paste("r[s > 0] ==", round(cor(tmp_df$InteractionScore.R1[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 > 0], 
                                                  tmp_df$InteractionScore.R2[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 > 0]), 3)), 
                                                                               parse = TRUE) +
  annotate("text", hjust = 0, x = -6, y = 7.6, 
           label = paste("r[s < 0] ==", round(cor(tmp_df$InteractionScore.R1[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 < 0], 
                                                  tmp_df$InteractionScore.R2[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 < 0]), 3)), 
                                                                               parse = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Nu Gene-Gene Interaction Score Replicate 1") +
  ylab("Nu Gene-Gene Interaction Score Replicate 2") +
  labs(color = "") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-8, 12, 4)) +
  scale_y_continuous(breaks = seq(-8, 12, 4))

# Add marginal histograms - FIGURE 2A
p2ah <- ggMarginal(p2a, groupColour = TRUE, groupFill = FALSE)


# Save plots
ggsave(snakemake@output[["output_figure_2a"]], 
       p2ah, device = "png", width = 5, height = 5, dpi = 300)