library(readxl)
library(edgeR)
library(survival)
library(openxlsx)
library(dplyr)
library(DESeq2)

ccounts = read.table('Data/TCGA/data_mrna_seq_v2_rsem.txt',sep='\t', header=TRUE)
ccounts <- na.omit(ccounts)
ccounts <- ccounts[apply(ccounts[,-1], 1, function(x) !all(x==0)),]
ccounts <- ccounts[!(is.na(ccounts$Hugo_Symbol) | ccounts$Hugo_Symbol==""), ]
ccounts$Entrez_Gene_Id <- NULL
ccounts$Hugo_Symbol <- gsub("-", "_", ccounts$Hugo_Symbol)

ccounts = aggregate(ccounts[c(2:153)], by=list(ccounts$Hugo_Sym), 
                    FUN=function(x) x=sum(x))

colnames(ccounts) = substr(colnames(ccounts),1,nchar(colnames(ccounts))-3)

metadata = read.table('Data/TCGA/data_clinical_patient.txt',sep='\t', header=TRUE)
metadata$Patient_ID <- gsub("-", ".", metadata$Patient_ID)

metadata = metadata[metadata$therapy %in% c('TMZ Chemoradiation, TMZ Chemo', 'Standard Radiation, TMZ Chemo', 'TMZ Chemo'),]


sample_ids <- colnames(ccounts)
metadata <- metadata[metadata$Patient_ID %in% sample_ids,]

sample_ids2 <- metadata$Patient_ID
sample_ids2 <- append(sample_ids2, "Grou")
counts <- ccounts[, (colnames(ccounts) %in% sample_ids2)]
rownames(counts) <- counts$Grou
counts$Grou <- NULL
counts <- round(counts)

gene_list <- scan('Overlap.txt', what="character")
counts <- counts[rownames(counts) %in% gene_list,]

intcount<-lapply(counts,as.integer)
intcount<-as.data.frame(intcount)
rownames(intcount) <- rownames(counts)
str(intcount)

# Set all counts of 4 or lower to 0
# Replace counts less than or equal to 4 with 0 for all samples
samples <- colnames(intcount)  # Extract sample column names
for (sample in samples) {
 intcount[[sample]] <- ifelse(intcount[[sample]] <= 4, 0, intcount[[sample]])
}

tissue.obj = DGEList(counts=as.matrix(intcount), genes=rownames(intcount))

# norm
tissue.cpm = cpm(tissue.obj, log=TRUE)
tissue.cpm = as.data.frame(tissue.cpm)

metadata$Overall_Survival_Status[metadata$Overall_Survival_Status == "0:LIVING"] <- 0
metadata$Overall_Survival_Status[metadata$Overall_Survival_Status == "1:DECEASED"] <- 1

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

metadata$Overall_Survival_Status <- as.numeric(as.character(metadata$Overall_Survival_Status))

tissue.cpm = as.data.frame(t(tissue.cpm))
tissue.gene_TCGA = survival_corr(tissue.cpm, metadata$OS_MONTHS,
                                 metadata$Overall_Survival_Status)
tissue.gene_TCGA = tissue.gene_TCGA[order(tissue.gene_TCGA$hazard.ratio),]
tissue.gene.01_TCGA = tissue.gene_TCGA[tissue.gene_TCGA$p.value < 0.1,]
tissue.gene.005_TCGA = tissue.gene_TCGA[tissue.gene_TCGA$p.value < 0.05,]
tissue.gene.005_TCGA$gene = rownames(tissue.gene.005_TCGA)

tissue.gene.001_TCGA = tissue.gene_TCGA[tissue.gene_TCGA$p.value < 0.01,]
tissue.gene.001_TCGA$gene = rownames(tissue.gene.001_TCGA)

write.csv(as.data.frame(tissue.gene.005_TCGA),
          file="Result_TCGA/CoxPH_TCGA_005.csv")

#####################################################################################################
#####################################################################################################
#####################################################################################################

df1 <- read.csv("Result_tissue/Coxph_genes_005.csv")
df2 <- read.csv("Result_TCGA/CoxPH_TCGA_005.csv")


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
          file="Result_TCGA/CoxPH_TCGA_005_checked.csv")
