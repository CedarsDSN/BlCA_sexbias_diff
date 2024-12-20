rm(list=ls())
library(ggplot2)
library(cowplot)


kegg<-read.csv("~/Desktop/male_hub_selected.csv",check.names = F,header = T,row.names = NULL,quote="",stringsAsFactors = FALSE)
kegg<-read.csv("~/Desktop/female_hub_genes_GO_selected.csv",check.names = F,header = T,row.names = NULL,quote="",stringsAsFactors = FALSE)




kegg[,2]<-paste(kegg[,2],kegg[,3],sep="~")
kegg<-kegg[,c(2,5,8,4,6,7)]
colnames(kegg)<-c("Term","Gene_Number","Gene","Fold_Enrichment_Score","pvalue","adjp")

kegg$adjp <- -log10(kegg$adjp)

kegg$Gene_Number<-as.numeric(kegg$Gene_Number)

kegg$Fold_Enrichment_Score<-as.numeric(kegg$Fold_Enrichment_Score)

kegg$Term <- reorder(kegg$Term, kegg$adjp)

# Make the plot
plot<-ggplot(kegg,
       aes(x=Fold_Enrichment_Score,
           y=Term,
           size=Gene_Number,
           colour=adjp))+
  geom_point()+
  scale_colour_gradient(low="red",
                        high="blue",
                        name="-log10(Qvalue)")+
  scale_size(name="Gene number") +
  xlab("Fold Enrichment Score")+
  ylab(NULL)+
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
#    axis.title.y = element_text(size = 2), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_blank()
  )


###+ggtitle("Female hub genes GO terms")


ggsave("Female_hub_genes_GO_terms.pdf", plot = plot, width = 16, height = 9)

ggsave("Male_hub_selected.png",plot = plot, device = "png",units="in", width=, height=9, dpi=300)
