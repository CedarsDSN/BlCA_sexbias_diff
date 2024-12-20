# Load necessary packages
#install.packages("survival")
#install.packages("survminer")
rm(list = ls())

library(survival)
library(survminer)
library(survivalROC)
library(ggplot2)
library(pROC)
setwd("~/Desktop/Bladder_TCGA_GTEx/Gender_diff/protein_coding/Survival_analysis/")
# Load the data
survival_data <- read.table("~/Desktop/Bladder_TCGA_GTEx/Gender_diff/protein_coding/Survival_analysis/TCGA-BLCA.survival.tsv", header=T)
expr_female <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Gender_diff/protein_coding/Survival_analysis/female_all_gene_tpm_tumor_only.csv", header=T, check.names=F, row.names = 1)
expr_male <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Gender_diff/protein_coding/Survival_analysis/male_all_gene_tpm_tumor_only.csv", header=T, check.names=F, row.names = 1)

expr_male<-expr_male[rownames(expr_female),]
expr_all<-cbind(expr_male,expr_female)


genes1 <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/PPI_network/compare/male_genes_info_hub.csv",header=T)
genes2 <- read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/PPI_network/compare/female_genes_info_hub.csv",header = T)

genes<-read.csv("~/Desktop/Bladder_TCGA_GTEx/Biomarkers/survival_analysis/male_hub/male_survival_sig_genes.csv",header=T)

#genes <- read.csv("~/Desktop/gene_list.txt")
genes<- as.vector(genes[,1])

length(genes)

expr_female <- expr_female[genes, ]
expr_male <- expr_male[genes, ]
expr_all <- expr_all[genes,]
dim(expr_female)
dim(expr_male)
dim(expr_all)

# Process survival data
survival_data[, 1] <- sub(".$", "", survival_data[, 1])
survival_female <- survival_data[match(colnames(expr_female), survival_data[, 1], nomatch = 0), ]
survival_male <- survival_data[match(colnames(expr_male), survival_data[, 1], nomatch = 0), ]
survival_all <- survival_data[match(colnames(expr_all), survival_data[, 1], nomatch = 0), ]
# Subset expression data
expr_female <- expr_female[, survival_female[, 1]]
expr_male <- expr_male[, survival_male[, 1]]
expr_all <- expr_all[, survival_all[, 1]]
# Function to perform survival analysis and return results

#  male_results <- surfun(gene, expr_male, survival_male, "male")

surfun <- function(gene_name, gender_matrix, survival_matrix, gender) {
  common_samples <- intersect(colnames(gender_matrix), survival_matrix[, 1])
  gender_matrix <- gender_matrix[, common_samples, drop = FALSE]
  survival_matrix <- survival_matrix[survival_matrix[, 1] %in% common_samples, ]
  
  expr_values <- gender_matrix[gene_name, ]
  survival_matrix$expr_values <- as.numeric(expr_values)
  survival_matrix <- survival_matrix[, c(1, 4, 2, 5)]
  colnames(survival_matrix) <- c("id", "futime", "fustat", gene_name)
  rt <- survival_matrix[, c("futime", "fustat", gene_name)]
  rt$futime <- rt$futime / 365 # Convert survival time to months
  group <- ifelse(rt[, 3] > median(rt[, 3]), "High", "Low")
  rt$group <- group
  
  # Log-rank test
  diff <- survdiff(Surv(futime, fustat) ~ group, data = rt)
  logrank_pValue <- 1 - pchisq(diff$chisq, df = 1)
  logrank_pValue_text <- paste0("Log-rank p=", sprintf("%.3f", logrank_pValue)) # Define logrank_pValue_text
  
  # Cox proportional hazards model
  cox_model <- coxph(Surv(futime, fustat) ~ group, data = rt)
  hr <- summary(cox_model)$coefficients[1, 2]
  hr_ci_lower <- summary(cox_model)$conf.int[1, 3]
  hr_ci_upper <- summary(cox_model)$conf.int[1, 4]
  hr_text <- paste0("HR=", sprintf("%.2f", hr), " (", sprintf("%.2f", hr_ci_lower), "-", sprintf("%.2f", hr_ci_upper), ")")
  
  fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
  
  surPlot <- ggsurvplot(fit, 
                        data = rt,
                        conf.int = FALSE,
                        pval = FALSE, # Disable default p-value
                        legend.labs = c("High", "Low"),
                        title=gene_name,
                        legend.title = "Expression",
                        legend=c(0.15,0.15),
                        xlab = "Time (years)",
                        break.time.by = 1, # Change break time to 6 months
                        risk.table.title = "",
                        palette = c("red", "black"),
                        risk.table = FALSE,
                        risk.table.height = .25,
                        font.x=c(20,"plain"),
                        font.y=c(20,"plain"),
                        font.title=c(25,"bold","red"),
                        font.legend=c(16),
                        ggtheme = theme_minimal() + 
                          theme(plot.title=element_text(hjust=0.5),
                            panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(), 
                                panel.border = element_rect(color = "black", fill = NA, size = 1)))
  
  # Add HR and log-rank p-value to the plot
  surPlot$plot <- surPlot$plot + 
    annotate("text", x = max(rt$futime) * 0.5, y = 0.8, label = hr_text, size = 7.5, hjust = 0) +
    annotate("text", x = max(rt$futime) * 0.5, y = 0.7, label = logrank_pValue_text, size = 7.5, hjust = 0)
  
  pdf(file = paste0(gender, "_", gene_name, "_survival_plot.pdf"), onefile = FALSE, width = 6, height = 5)
  print(surPlot)
  dev.off()
  
  # Generate and save ROC plot
  roc <- survivalROC(Stime = rt$futime, status = rt$fustat, marker = rt[, gene_name], predict.time = 5, method = "KM") # predict.time is now 5 years
  
  # Check if ROC AUC calculation is correct
  if (!is.na(roc$AUC) && roc$AUC >= 0.5 && roc$AUC <= 1) {
    auc <- roc$AUC
  } else {
    auc <- NA
  }
  
  par(pty="s")
  roc(rt$fustat,as.numeric(rt[, gene_name]), plot=TRUE, legacy.axes=TRUE, percent=FALSE, xlab="False Positive Percentage", ylab="True Postive Percentage",
      col=palette_14[1],lwd=4, print.auc=TRUE,print.auc.x=45)
  
  
  roc(rt$fustat,as.numeric(rt[, 4]), plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage",
      col=palette_14[1],lwd=4, print.auc=TRUE,print.auc.y=45)

  
  for (i in 2:14){
    roc(rt$fustat,as.numeric(rt[, 3+i]), plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage",
        col=palette_14[i],lwd=4, print.auc=TRUE,print.auc.y=47-i-i,add=TRUE)
  }
  
  legend("bottomright",legend=colnames(rt)[4:17],col=palette_14,lwd=4)
  
  
  roc(rt$fustat,rt[,5], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage",
      col=palette_14[1],lwd=4, print.auc=TRUE,print.auc.x=45)
  
  
  roc_plot <- ggplot(data.frame(FP = roc$FP, TP = roc$TP), aes(x = FP, y = TP)) +
    geom_line(color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = paste0(gene_name, " ROC Curve"), x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA, size = 1))
  
#  pdf(file = paste0(gender, "_", gene_name, "_ROC_plot.pdf"), onefile = FALSE, width = 6, height = 5)
  print(roc_plot)
  dev.off()
  
  return(c(logrank_pValue, hr, auc))
}

# Initialize a data frame to store the results
results <- data.frame(
  gene_name = character(),
  male_pvalue = numeric(),
  female_pvalue = numeric(),
  all_pvalue = numeric(),
  male_hr = numeric(),
  female_hr = numeric(),
  all_hr = numeric(),
  male_auc = numeric(),
  female_auc = numeric(),
  all_auc = numeric(),
  stringsAsFactors = FALSE
)

# Run survival analysis for selected genes

genes_to_analyze<-genes


for (gene in genes_to_analyze) {
  female_results <- surfun(gene, expr_female, survival_female, "female")
  male_results <- surfun(gene, expr_male, survival_male, "male")
  all_results <- surfun(gene, expr_all, survival_all, "all")
  results <- rbind(results, data.frame(
    gene_name = gene,
    male_pvalue = male_results[1],
    female_pvalue = female_results[1],
    all_pvaule = all_results[1],
    male_hr = male_results[2],
    female_hr = female_results[2],
    all_hr = all_results[2],
    male_auc = male_results[3],
    female_auc = female_results[3],
    all_auc = all_results[3]
  ))
}

# Save the results to a CSV file
write.csv(results, file = "all_survival_analysis_results.csv", row.names = FALSE)
