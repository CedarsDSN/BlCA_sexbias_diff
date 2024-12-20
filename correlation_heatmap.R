rm(list=ls())
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(cowplot)
library(grid)

setwd("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/")

# Load the data
group_A <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/female_ciber.csv", row.names = 1,check.names=F)
group_C <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/female_all_gene_tpm_tumor_only.csv", row.names = 1,check.names=F)


group_A<-group_A[-c(2,5,10,21,22),]

group_B <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/male_ciber.csv", row.names = 1,check.names=F)
group_D <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/male_all_gene_tpm_tumor_only.csv", row.names = 1,check.names=F)

group_B<-group_B[-c(2,5,10,21,22),]

genes_male<-read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/survival_analysis/male_hub/male_survival_sig_genes.csv",header=T)
genes_female<-read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/survival_analysis/female_hub/female_survival_sig_genes.csv",header=T)

genes<-c(genes_male[,1],genes_female[,1])

genes<-c("PDGFRA","DAXX","PPARG","IKBKB")

group_C<-group_C[genes,]
group_D<-group_D[genes,]

group_A<-group_A[,colnames(group_C)]
group_B<-group_B[,colnames(group_D)]

dim(group_A)
dim(group_C)
dim(group_B)
dim(group_D)


# Function to calculate correlation and p-value using Spearman method
calculate_correlations <- function(immune_cells, genes) {
  cor_results <- data.frame()
  for (cell in rownames(immune_cells)) {
    for (gene in rownames(genes)) {
      cor_test <- cor.test(as.numeric(immune_cells[cell, ]), as.numeric(genes[gene, ]), method = "spearman")
      cor_value <- cor_test$estimate
      p_value <- cor_test$p.value
      cor_results <- rbind(cor_results, data.frame(Cell = cell, Gene = gene, Correlation = cor_value, PValue = p_value))
    }
  }
  return(cor_results)
}

# Calculate correlations for groups A (immune cells) and C (gene expressions)
cor_AC <- calculate_correlations(group_A, group_C)
cor_AC$Group <- "AC"

# Calculate correlations for groups B (immune cells) and D (gene expressions)
cor_BD <- calculate_correlations(group_B, group_D)
cor_BD$Group <- "BD"

# Function to perform Fisher's z-test and include additional columns
fisher_z_test <- function(cor_AC, cor_BD, n_AC, n_BD) {
  cor_results <- data.frame()
  for (i in 1:nrow(cor_AC)) {
    z_AC <- 0.5 * log((1 + cor_AC$Correlation[i]) / (1 - cor_AC$Correlation[i]))
    z_BD <- 0.5 * log((1 + cor_BD$Correlation[i]) / (1 - cor_BD$Correlation[i]))
    SE <- sqrt(1 / (n_AC - 3) + 1 / (n_BD - 3))
    z_diff <- (z_AC - z_BD) / SE
    p_value <- 2 * (1 - pnorm(abs(z_diff)))
    cor_results <- rbind(cor_results, data.frame(
      Cell = cor_AC$Cell[i], 
      Gene = cor_AC$Gene[i], 
      Correlation_AC = cor_AC$Correlation[i], 
      PValue_AC = cor_AC$PValue[i], 
      Correlation_BD = cor_BD$Correlation[i], 
      PValue_BD = cor_BD$PValue[i], 
      z_diff = z_diff, 
      p_value = p_value
    ))
  }
  return(cor_results)
}

# Perform Fisher's z-test
n_AC <- ncol(group_A)
n_BD <- ncol(group_B)
fisher_results <- fisher_z_test(cor_AC, cor_BD, n_AC, n_BD)

# Save correlation results to CSV with group identifiers and p-values
correlation_results <- rbind(cor_AC, cor_BD)
write.csv(correlation_results, "correlation_results.csv", row.names = FALSE)

# Save Fisher's z-test results to CSV with additional columns
write.csv(fisher_results, "fisher_results.csv", row.names = FALSE)

# Output significant results
significant_results <- fisher_results %>% filter(p_value < 0.05)
print(significant_results)

library(pheatmap)
library(reshape2)

# Load correlation results
#correlation_results <- read.csv("correlation_results.csv")

# Prepare the data for heatmap for group A&C
cor_AC <- correlation_results %>% filter(Group == "AC")
cor_AC_matrix <- dcast(cor_AC, Cell ~ Gene, value.var = "Correlation")
rownames(cor_AC_matrix) <- cor_AC_matrix$Cell
cor_AC_matrix <- cor_AC_matrix[,-1]

# Prepare the data for heatmap for group B&D
cor_BD <- correlation_results %>% filter(Group == "BD")
cor_BD_matrix <- dcast(cor_BD, Cell ~ Gene, value.var = "Correlation")
rownames(cor_BD_matrix) <- cor_BD_matrix$Cell
cor_BD_matrix <- cor_BD_matrix[,-1]



# Create significance matrix for uniquely significant in group A&C
sig_unique_AC <- fisher_results %>%
#  mutate(Significant = ifelse((PValue_AC < 0.05 & PValue_BD >= 0.05), #####unique significant or show all
  mutate(Significant = ifelse((PValue_AC < 0.05),                            
                              ifelse(PValue_AC < 0.01, "*", ""), "")) %>%
  dcast(Cell ~ Gene, value.var = "Significant")
rownames(sig_unique_AC) <- sig_unique_AC$Cell
sig_unique_AC <- sig_unique_AC[,-1]

# Create significance matrix for uniquely significant in group B&D
sig_unique_BD <- fisher_results %>%
#  mutate(Significant = ifelse((PValue_BD < 0.05 & PValue_AC >= 0.05), #####unique significant or show all
  mutate(Significant = ifelse((PValue_BD < 0.05),
                              ifelse(PValue_BD < 0.01, "*", ""), "")) %>%
  dcast(Cell ~ Gene, value.var = "Significant")
rownames(sig_unique_BD) <- sig_unique_BD$Cell
sig_unique_BD <- sig_unique_BD[,-1]

# Function to add asterisks only for significant correlations
add_asterisks <- function(mat, sig_mat) {
  mat_annotated <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      if (sig_mat[i, j] != "") {
        mat_annotated[i, j] <- sig_mat[i, j]
      }
    }
  }
  rownames(mat_annotated) <- rownames(mat)
  colnames(mat_annotated) <- colnames(mat)
  return(mat_annotated)
}

# Annotate matrices with asterisks only for significant correlations
cor_AC_annotated <- add_asterisks(cor_AC_matrix, sig_unique_AC)
cor_BD_annotated <- add_asterisks(cor_BD_matrix, sig_unique_BD)

heatmap_AC <- pheatmap(cor_AC_matrix, 
                       #main = "Female", 
                       cluster_rows = FALSE, cluster_cols = FALSE, 
                       show_rownames = TRUE, show_colnames = TRUE,
                       breaks = seq(-0.5, 0.5, length.out = 100),
                       display_numbers = cor_AC_annotated, number_color = "black", fontsize_number = 30,
                       fontsize_row = 18,
                       fontsize_col = 18,
                       legend = FALSE)

# Generate the heatmap for group B&D with the legend
heatmap_BD <- pheatmap(cor_BD_matrix, 
                       #main = "Male", 
                       cluster_rows = FALSE, cluster_cols = FALSE, 
                       show_rownames = TRUE, show_colnames = TRUE,
                       breaks = seq(-0.5, 0.5, length.out = 100),
                       display_numbers = cor_BD_annotated, number_color = "black", fontsize_number = 30,
                       fontsize_row = 18,
                       fontsize_col = 18,
                       #fontsize = 20,
                       legend = FALSE)

pdf("female_temp.pdf")
print(heatmap_BD)
dev.off()

# Extract the legend
legend <- cowplot::get_legend(heatmap_BD$gtable)

# Create custom legend for the significance annotation
custom_legend <- grobTree(
  textGrob("* P < 0.01", x=0.92, y=7.5, hjust=0, gp=gpar(fontsize=18))
)

# Arrange the heatmaps and legends
combined <- grid.arrange(
  arrangeGrob(heatmap_AC$gtable, heatmap_BD$gtable, ncol = 2),
#  legend,
#  custom_legend,
  ncol = 1,
  heights = c(10, 0, 0)
)

# Add the main title
final_plot <- arrangeGrob(
  textGrob("Spearman correlation of the immune related key genes and the immune cells", gp = gpar(fontsize = 25, fontface = "plain")),
  combined,
  ncol = 1,
  heights = c(1, 10)
)

# Draw the final plot
grid.newpage()
grid.draw(final_plot)

# Save the combined plot to a file if needed
ggsave("combined_heatmaps.pdf", plot = final_plot, width = 15, height = 10)

# Generate the heatmap for group A&C without clustering and with custom scale
#pheatmap(cor_AC_matrix, main = "Spearman correlation of the key genes and the immune cells (Female)", 
#         cluster_rows = FALSE, cluster_cols = FALSE, 
#         show_rownames = TRUE, show_colnames = TRUE,
#         breaks = seq(-0.5, 0.5, length.out = 100),
#         display_numbers = cor_AC_annotated, number_color = "black", fontsize_number = 30)

# Generate the heatmap for group B&D without clustering and with custom scale
#pheatmap(cor_BD_matrix, main = "Spearman correlation of the key genes and the immune cells (Male)", 
#         cluster_rows = FALSE, cluster_cols = FALSE, 
#         show_rownames = TRUE, show_colnames = TRUE,
#         breaks = seq(-0.5, 0.5, length.out = 100),
#         display_numbers = cor_BD_annotated, number_color = "black", fontsize_number = 30)
