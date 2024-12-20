rm(list=ls())
library(VennDiagram)

# Load the data directly
female_data <- read.csv("~/Desktop//female_RNAseq_WXS_mf002_p005_info.csv")
male_data <- read.csv("~/Desktop//male_RNAseq_WXS_mf002_p005_info.csv")

# Extract upregulated and downregulated genes
female_upregulated <- female_data[female_data$log2FoldChange > 0, "Gene"]
female_downregulated <- female_data[female_data$log2FoldChange < 0, "Gene"]
male_upregulated <- male_data[male_data$log2FoldChange > 0, "Gene"]
male_downregulated <- male_data[male_data$log2FoldChange < 0, "Gene"]

venn.plot <- venn.diagram(
  x = list(
    "Female Upregulated" = female_upregulated,
    "Female Downregulated" = female_downregulated,
    "Male Upregulated" = male_upregulated,
    "Male Downregulated" = male_downregulated
  ),
  category.names = c("Female Upregulated", "Female Downregulated", "Male Upregulated", "Male Downregulated"),
  filename = NULL,
  output = TRUE,
  col = "transparent",  # Remove ellipse outlines
  fill = c("red", "blue", "green", "yellow"),
  alpha = 0.5,
  label.col = "black",  # Set font color to black
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = "black",  # Set category font color to black
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  margin = 0.1
)

# Plot the Venn diagram
grid.draw(venn.plot)
