library(readxl)
library(edgeR)
library(survival)
library(openxlsx)
library(dplyr)
library(DESeq2)

# Validation of 121 overlapping genes in external Gliotrain dataset
RNAseq = read.csv("Data/Gliotrain/RNA_Gliotrain_complete_withoutEMC.csv")
rownames(RNAseq) <- RNAseq$gene
RNAseq$gene <- NULL
RNAseq[RNAseq == "NA"] <- ""
RNAseq <- round(RNAseq)
RNAseq[is.na(RNAseq)] <- 0

metadata <- read.table('Data/Gliotrain/Metadata_Gliotrain.txt',sep='\t', header=TRUE)
metadata <- metadata[!grepl("GT.03", metadata$patient_id_gliotrain),]
metadata = metadata[metadata$IDH_status %in% 'WT',]

sample_ids <- colnames(RNAseq)
metadata <- metadata[metadata$patient_id_gliotrain %in% sample_ids,]

gene_list <- scan('Overlap.txt', what="character")
RNAseq <- RNAseq[rownames(RNAseq) %in% gene_list,]

# Set all counts of 4 or lower to 0
# Replace counts less than or equal to 4 with 0 for all samples
samples <- colnames(RNAseq)  # Extract sample column names
for (sample in samples) {
 RNAseq[[sample]] <- ifelse(RNAseq[[sample]] <= 4, 0, RNAseq[[sample]])
}

tissue.obj = DGEList(counts=as.matrix(RNAseq), genes=rownames(RNAseq))

# norm
norm.factorTissue = calcNormFactors(tissue.obj, method = "TMM")
tissue.cpm = cpm(norm.factorTissue, log=TRUE)
tissue.cpm = as.data.frame(tissue.cpm)

metadata$overall_survival_event[metadata$overall_survival_event == "alive"] <- 0
metadata$overall_survival_event[metadata$overall_survival_event == "dead"] <- 1

survival_corr <- function(dataset,time,events){
  tumor_survival = sapply(dataset, function(x) summary(coxph(Surv(time, events) ~ x )))
  # hazard ratio
  tumor_survival_HR = sapply(tumor_survival["coefficients",], function(x) x[[2]])
  # p value
  tumor_survival_p = sapply(tumor_survival["coefficients",], function(x) x[[5]])
  
  tumor_survival_all = data.frame('hazard ratio' = tumor_survival_HR,
                                  # 'stat' = tumor_survival_stat,
                                  'p-value' = tumor_survival_p)
  tumor_survival_all$p.adj = p.adjust(tumor_survival_all$p.value, method = "BH")
  
  
  return(tumor_survival_all)
}

metadata$overall_survival_event <- as.numeric(as.character(metadata$overall_survival_event))

tissue.cpm = as.data.frame(t(tissue.cpm))


tissue.gene_TCGA = survival_corr(tissue.cpm, metadata$OS_months,
                                 metadata$overall_survival_event)
tissue.gene_TCGA = tissue.gene_TCGA[order(tissue.gene_TCGA$hazard.ratio),]
tissue.gene.01_TCGA = tissue.gene_TCGA[tissue.gene_TCGA$p.value < 0.1,]
tissue.gene.005_TCGA = tissue.gene_TCGA[tissue.gene_TCGA$p.value < 0.05,]
tissue.gene.005_TCGA$gene = rownames(tissue.gene.005_TCGA)
tissue.gene_TCGA$gene = rownames(tissue.gene_TCGA)

tissue.gene.001_TCGA = tissue.gene_TCGA[tissue.gene_TCGA$p.value < 0.01,]
tissue.gene.001_TCGA$gene = rownames(tissue.gene.001_TCGA)

write.csv(as.data.frame(tissue.gene.005_TCGA),
          file="Result_Gliotrain/CoxPH_gliotrain_005.csv")


#####################################################################################################
#####################################################################################################
#####################################################################################################

# Check whether significant genes point into same direction
df1 <- read.csv("Result_tissue/Coxph_genes_005.csv")
df2 <- read.csv("Result_Gliotrain/CoxPH_gliotrain_005.csv")


df1 <- df1[order(df1$X), ]
df2 <- df2[order(df2$X), ]

common_genes <- intersect(df1$X, df2$X)

df1_common <- filter(df1, X %in% common_genes)
df2_common <- filter(df2, X %in% common_genes)

all(df1_common$X == df2_common$X)

# Check log2FoldChange values and remove genes with opposite signs
common_genes_filtered <- common_genes[!(df1_common$hazard.ratio < 1 & df2_common$hazard.ratio < 1) &
                                        !(df1_common$hazard.ratio > 1 & df2_common$hazard.ratio > 1)]

df2_filtered <- df2_common[!(df2_common$X %in% common_genes_filtered),]

write.csv(as.data.frame(df2_filtered),
          file="Result_Gliotrain/CoxPH_gliotrain_005_checked.csv")
