library(readr)
library(ggplot2)
library(dplyr)

# Load data
gamma_gene_gi <- read_tsv(snakemake@input[["input_gamma"]])
tau_gene_gi <- read_tsv(snakemake@input[["input_tau"]])
nu_gene_gi <- read_tsv(snakemake@input[["input_nu"]])


gi_cmp <- inner_join(gamma_gene_gi[,c(1:6,12)], tau_gene_gi[,c(1:6,12)], by = colnames(gamma_gene_gi)[1:5], 
                    suffix = c(".Gamma", ".Tau"))
gi_cmp <- inner_join(gi_cmp, nu_gene_gi[,c(1:6, 12)], by = colnames(gi_cmp)[1:5],)
colnames(gi_cmp)[10:11] <- c("InteractionScore.Nu", "Sig.Nu")

gi_cmp2 <- gi_cmp %>% filter(Category == "X+Y")

gi_cmp2$Highlight <- NA
gi_cmp2$Highlight[gi_cmp2$InteractionScore.Nu > 0 & gi_cmp2$Sig.Nu] <- "Positive"
gi_cmp2$Highlight[gi_cmp2$InteractionScore.Nu < 0 & gi_cmp2$Sig.Nu] <- "Negative"
gi_cmp2$Highlight[is.na(gi_cmp2$Highlight)] <- "None"
gi_cmp2$Highlight <- factor(gi_cmp2$Highlight, levels = c("None",
                                                          "Positive", "Negative"))


gi_cmp2 <- gi_cmp2[order(gi_cmp2$Highlight),]

# Plot - FIGURE 3A
p3a <- ggplot(gi_cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) +
  geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_point(alpha = .5, size = .8) +
  geom_point(data = gi_cmp2[gi_cmp2$Highlight != "None",], 
             aes(color = Highlight), alpha = .5, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  coord_fixed() +
  scale_x_continuous(limits = c(-12,8), breaks = seq(-12, 8, 4)) +
  scale_y_continuous(limits = c(-11,16.5), breaks = seq(-8, 16, 4)) +
  xlab("") + ylab("") +
  scale_color_manual(values = c("#D81B60", "#33716B")) +
  xlab("Gamma IS") + ylab("Tau IS") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

print("Made figure 3a plot")


# plt_lab <- ggplot(gi.samp[gi.samp$Category == "X+Y",], aes(InteractionScore.Gamma, InteractionScore.Tau)) +
#   geom_abline(alpha = .7) +
#   geom_hline(yintercept = 0, alpha = 0.7) +
#   geom_vline(xintercept = 0, alpha = 0.7) +
#   geom_point(alpha = .5, size = .8) +
#   geom_point(data = gi.cmp2[gi.cmp2$Highlight != "None",], aes(color = Highlight), alpha = .5, size = 1.2) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.grid = element_blank(),
#         axis.ticks.length = unit(.2, "cm"),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line()) + coord_fixed() +
#   scale_x_continuous(limits = c(-12,8), breaks = seq(-12, 8, 1)) +
#   scale_y_continuous(limits = c(-11,16.5), breaks = seq(-8, 16, 1)) +
#   scale_color_manual(values = c("#D81B60", "#33716B")) +
#   xlab("Gamma IS") + ylab("Tau IS") +
#   guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))


ggsave(snakemake@output[["output_figure_3a"]], p3a,
       device = "png", height = 5, width = 3)
