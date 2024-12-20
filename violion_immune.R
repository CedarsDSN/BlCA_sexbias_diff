rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggsignif)

female_ciber <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/female_ciber.csv", row.names = 1,check.names=F)
male_ciber <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/male_ciber.csv", row.names = 1,check.names = F)

female <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/female_all_gene_tpm_tumor_only.csv", row.names = 1,check.names=F)
male <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/immune_filtration/cibersort/male_all_gene_tpm_tumor_only.csv", row.names = 1,check.names=F)

female_ciber<-female_ciber[,colnames(female)]
male_ciber<-male_ciber[,colnames(male)]

female_ciber$Immune_Cell_Type <- rownames(female_ciber)
male_ciber$Immune_Cell_Type <- rownames(male_ciber)

# Melt the data to long format
female_ciber_long <- melt(female_ciber, id.vars = "Immune_Cell_Type", variable.name = "Sample", value.name = "Proportion")
female_ciber_long$Gender <- "Female Tumor"
male_ciber_long <- melt(male_ciber, id.vars = "Immune_Cell_Type", variable.name = "Sample", value.name = "Proportion")
male_ciber_long$Gender <- "Male Tumor"

# Combine the data
combined_ciber_long <- rbind(female_ciber_long, male_ciber_long)

# Perform the Wilcox test for each immune cell type
immune_cell_types <- unique(combined_ciber_long$Immune_Cell_Type)
wilcox_results <- data.frame(Immune_Cell_Type = character(), P_Value = numeric(), stringsAsFactors = FALSE)

for (cell_type in immune_cell_types) {
  female_values <- combined_ciber_long %>% filter(Immune_Cell_Type == cell_type & Gender == "Female Tumor") %>% pull(Proportion)
  male_values <- combined_ciber_long %>% filter(Immune_Cell_Type == cell_type & Gender == "Male Tumor") %>% pull(Proportion)
  test_result <- wilcox.test(female_values, male_values)
  wilcox_results <- rbind(wilcox_results, data.frame(Immune_Cell_Type = cell_type, P_Value = test_result$p.value))
}

# Print the Wilcox test results
print(wilcox_results)

# Filter for significant p-values
significant_results <- wilcox_results %>% filter(P_Value < 0.05)

# Calculate y positions for the p-values
y_positions <- combined_ciber_long %>% group_by(Immune_Cell_Type) %>% summarise(y_pos = max(Proportion) + 0.05)

# Merge y positions with significant results
significant_results <- merge(significant_results, y_positions, by = "Immune_Cell_Type")

# Generate the violin plot
p <- ggplot(combined_ciber_long, aes(x = Immune_Cell_Type, y = Proportion, fill = Gender)) +
  geom_violin(scale = "width", adjust = 1, width = 0.7) +
  scale_fill_manual(values = c("Female Tumor" = "red", "Male Tumor" = "blue")) +
  labs(title = NULL, x = NULL, y = "Proportion")+
#  labs(title = "Violin Plot of Immune Cell Proportions", x = "Immune Cell Type", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.position = "none")

# Add asterisks and p-values to the plot for significant results
p <- p + geom_text(data = significant_results, aes(x = Immune_Cell_Type, y = y_pos, label = "*"), inherit.aes = FALSE, vjust = -0.5, size = 9) +
  geom_text(data = significant_results, aes(x = Immune_Cell_Type, y = y_pos + 0.05, label = sprintf("p = %.3f", P_Value)), inherit.aes = FALSE, vjust = -0.5, size = 9)

# Print the plot
print(p)
