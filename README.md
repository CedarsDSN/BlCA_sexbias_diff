## Analysis
### Workflow Overview

graph LR;

subgraph TCGA bulkRNA-seq
    A1(TCGA bulkRNA-seq) --> B1(Identify DE Genes <br> (DESeq2, M vs F));
    B1 --> C1(Pathway analysis);
    B1 --> D1(Survival analysis <br> (M vs F));
end

subgraph TCGA miRNA-seq
    A2(TCGA miRNA-seq) --> B2(Identify DE miRNAs <br> (DESeq2, M vs F));
    B2 --> C2(Survival analysis <br> (M vs F));
    B2 --> D2(Predict DE miRNA target genes <br> (miRNAtap));
    D2 --> E2(Filter target genes <br> (neg correlation & DE));
    E2 --> F2(Pathway analysis);
    E2 --> G2(Target gene scores);
end

subgraph scRNA-seq Analysis
    A3(BLCA scRNA-seq) --> B3(QC, Integration, <br> Cell type ID <br> (Seurat));
    B3 --> C3(Average expression <br> of target genes <br> in each cell type);
    C3 --> D3(Survival analysis <br> (M vs F));
    D3 --> E3(Expression diff <br> (M vs F));
    D3 --> F3(Hazard ratio diff);
    C3 --> G3(Immune-related <br> expression patterns);
    C3 --> H3(Z-score normalization <br> across cell types);
end

subgraph Validation
    I1(Independent GEO mRNA-seq) --> J1(Validation of <br> key findings);
end
