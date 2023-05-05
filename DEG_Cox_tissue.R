library(readxl)
library(edgeR)
library(survival)
library(openxlsx)
library(dplyr)
library(DESeq2)

##### tissue #####
outcomes.t = read_excel("Data/Tissue/Metadata_Tissue.xlsx", 1)
sampleNames.t = outcomes.t$GS.number

##### read data #####
tissue.count = read.table("Data/Tissue/RNASeq_rawCounts-EMC.txt",sep = "\t", header = TRUE)#, row.names = 1
meta.t = read_excel("Data/Tissue/GLIOTRAIN list_Labeling.xlsx", 1)
meta.t = na.omit(meta.t)
meta.t$GsID = paste0("GS.",meta.t$GsID)
colnames(tissue.count)[c(2:57)] = meta.t$GsID

#### Aggregate
tissue.count.aggre = aggregate(tissue.count[c(2:57)], by=list(tissue.count$Gene_Symbol), 
                               FUN=function(x) x=sum(x))
rownames(tissue.count.aggre) = tissue.count.aggre$Group.1
tissue.count.aggre$Group.1 = NULL

tissue.count20 <- tissue.count.aggre[, colnames(tissue.count.aggre) %in% sampleNames.t]


# Define the responders and non-responders group / make list of corresponding sample names
LTS <- outcomes.t[outcomes.t$Survival=="LTS", c("GS.number", "Survival")]
STS <- outcomes.t[outcomes.t$Survival=="STS", c("GS.number", "Survival")]

sampleNames.L = LTS$GS.number
sampleNames.S = STS$GS.number

# Create dataframe with each samplename corresponding to the right group
dft = data.frame()
dft <- rbind(dft, as.data.frame(sampleNames.t))
dft$condition <- outcomes.t$Survival
rownames(dft) <- dft$sampleNames.t
dft$sampleNames.t = NULL

# Remove rows containing only 0 (to compare with threshold)
filtered_dft_naomit <- tissue.count20 %>% filter_all(any_vars(. != 0))


# Set all counts of 4 or lower to 0
# Replace counts less than or equal to 4 with 0 for all samples
samples <- colnames(tissue.count20)  # Extract sample column names
for (sample in samples) {
  tissue.count20[[sample]] <- ifelse(tissue.count20[[sample]] <= 4, 0, tissue.count20[[sample]])
}

# Filter step, filter rows which have less than 5% counts in both groups (responders or nonresponders) --> keep rows with more than 5% counts in 1 of the groups 
# Calculate the percentage of non-zero samples separately for the first two columns and the last two columns
tissue.count20$non_zero_precentage_LTS <- rowMeans(tissue.count20[, sampleNames.L] > 0) * 100
tissue.count20$non_zero_precentage_STS <- rowMeans(tissue.count20[, sampleNames.S] > 0) * 100

# Set the threshold value
threshold <- 10

# Filter the rows where the non-zero percentage is higher than the threshold in at least one of the sets of columns
filtered_dft <- tissue.count20 %>% filter(non_zero_precentage_LTS > threshold | non_zero_precentage_STS > threshold)
filtered_dft$non_zero_precentage_LTS = NULL
filtered_dft$non_zero_precentage_STS = NULL

# Check whether the samplenames are in the same order
filtered_dft <- filtered_dft[, rownames(dft)]
all(rownames(dft) == colnames(filtered_dft))

'MGMT' %in% rownames(filtered_dft)
'TERT' %in% rownames(filtered_dft)

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_dft,
                              colData = dft,
                              design = ~ condition)


# Define the reference level 
dds$condition <- factor(dds$condition, levels = c("LTS","STS"))

# DEG analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","LTS","STS"))

resultsNames(dds)

resOrdered <- res[order(res$pvalue),]

# Only keep significant genes
DE_genes <- subset(resOrdered, pvalue < 0.05)
DE_genes_Padjust <- subset(resOrdered, padj < 0.05)
DE_genes_005 <- subset(resOrdered, pvalue < 0.05)

# Define the up and down regulated genes within the responders group
log2FC <- DE_genes$log2FoldChange
upregulated_genes <- subset(DE_genes, log2FoldChange > 1)
downregulated_genes <- subset(DE_genes, log2FoldChange < -1)

# Plot the MGMT counts per group
plotCounts(dds, gene="MGMT", intgroup="condition")
plotCounts(dds, gene="TERT", intgroup="condition")

# Save analysis
write.csv(as.data.frame(DE_genes_005),
          file="Result_tissue/DE_genes_tissue_005.csv")

#####################################################################################################
#####################################################################################################
#####################################################################################################

##### tissue #####
tissue.count = read.table("Data/Tissue/RNASeq_rawCounts-EMC.txt",sep = "\t", header = TRUE)#, row.names = 1
meta.t = read_excel("Data/Tissue/GLIOTRAIN list_Labeling.xlsx", 1)
meta.t = na.omit(meta.t)
sampleNames.t = colnames(tissue.count)[c(2:57)]
sampleNames.t = gsub(".R","",sampleNames.t)
sampleNames.t = gsub(".","-",sampleNames.t,fixed = TRUE)

outcomes.t = read_excel("Data/Tissue/gliotrain survival data_IN140423.xlsx")
outcomes.t$patient_id_gliotrain = gsub(".","-",outcomes.t$patient_id_gliotrain,fixed = TRUE)
sampleNames.t = gsub("-01-WT-02","",sampleNames.t)
outcomes.t = outcomes.t[outcomes.t$patient_id_gliotrain %in% sampleNames.t,]
# outcomes.t$patient_id_gliotrain == sampleNames.t
meta.t$GsID = paste0("GS.",meta.t$GsID)

colnames(tissue.count)[c(2:57)] = meta.t$GsID
outcomes.t$patient_id_gliotrain = meta.t$GsID

outcomes.t$overall_survival_event = as.factor(outcomes.t$overall_survival_event)
outcomes.t$overall_survival_event = as.numeric(outcomes.t$overall_survival_event)
outcomes.t$OS_months = as.numeric(outcomes.t$OS_months)

# aggregate 
tissue.count.aggre = aggregate(tissue.count[c(2:57)], by=list(tissue.count$Gene_Symbol), 
                               FUN=function(x) x=sum(x))
rownames(tissue.count.aggre) = tissue.count.aggre$Group.1
tissue.count.aggre$Group.1 = NULL

gene_list <- scan('Result_tissue/DE_genes_tissue_005.txt', what="character")
tissue.count.aggre <- tissue.count.aggre[rownames(tissue.count.aggre) %in% gene_list,]


# Set all counts of 4 or lower to 0
# Replace counts less than or equal to 4 with 0 for all samples
samples <- colnames(tissue.count.aggre)  # Extract sample column names
for (sample in samples) {
  tissue.count.aggre[[sample]] <- ifelse(tissue.count.aggre[[sample]] <= 4, 0, tissue.count.aggre[[sample]])
}

tissue.obj = DGEList(counts=as.matrix(tissue.count.aggre), genes=rownames(tissue.count.aggre))

# norm
norm.factorTissue = calcNormFactors(tissue.obj, method = "TMM")
tissue.cpm = cpm(norm.factorTissue, log=TRUE)
tissue.cpm = as.data.frame(tissue.cpm)

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

tissue.cpm = as.data.frame(t(tissue.cpm))
tissue.gene = survival_corr(tissue.cpm, outcomes.t$OS_months,
                            outcomes.t$overall_survival_event)
tissue.gene = tissue.gene[order(tissue.gene$hazard.ratio),]
tissue.gene.01 = tissue.gene[tissue.gene$p.value < 0.1,]
tissue.gene.005 = tissue.gene[tissue.gene$p.value < 0.05,]
tissue.gene.005$gene = rownames(tissue.gene.005)

tissue.gene.001 = tissue.gene[tissue.gene$p.value < 0.01,]
tissue.gene.001$gene = rownames(tissue.gene.001)

write.csv(as.data.frame(tissue.gene.005),
          file="Result_tissue/Coxph_genes_005.csv")

#####################################################################################################
#####################################################################################################
#####################################################################################################

df1 <- read.csv("Result_cell/Correlation_genes_005.csv")
df2 <- read.csv("Result_tissue/Coxph_genes_005.csv")


df1 <- df1[order(df1$gene), ]
df2 <- df2[order(df2$gene), ]

common_genes <- intersect(df1$gene, df2$gene)

df1_common <- filter(df1, gene %in% common_genes)
df2_common <- filter(df2, gene %in% common_genes)

all(df1_common$gene == df2_common$gene)

# Check log2FoldChange values and remove genes with opposite signs
common_genes_filtered <- common_genes[!(df1_common$coefficient > 0 & df2_common$hazard.ratio > 1) &
                                        !(df1_common$coefficient < 0 & df2_common$hazard.ratio < 1)]

df1_filtered <- df1_common[!(df1_common$gene %in% common_genes_filtered),]
df2_filtered <- df2_common[!(df2_common$gene %in% common_genes_filtered),]

write.csv(as.data.frame(df1_filtered),
          file="Result_cell/Correlation_genes_005_overlap.csv")
write.csv(as.data.frame(df2_filtered),
          file="Result_tissue/Coxph_genes_005_overlap.csv")


##################################################
##################################################
##################################################

##### heatmap tumor19 OS #####
outcomes.t.order = outcomes.t[order(outcomes.t$OS_months,decreasing = TRUE),]
outcomes.t.order$MGMT_methylation_status[outcomes.t.order$MGMT_methylation_status == "methylated"] = "Methylated"
outcomes.t.order$MGMT_methylation_status[outcomes.t.order$MGMT_methylation_status == "unmethylated"] = "Unmethylated"

tissue56.sig = as.data.frame(t(tissue.cpm))[tissue.gene.005$gene,]

tissue56.sig = tissue56.sig[outcomes.t.order$patient_id_gliotrain]


annotation_col  = data.frame(OverallSurvival=outcomes.t.order$OS_months,
                             MGMTtumorStatus = outcomes.t.order$MGMT_methylation_status
)
rownames(annotation_col) = colnames(tissue56.sig)

ann_colors = list(OverallSurvival= c("deeppink4","chocolate2","darkseagreen"),
                  MGMTtumorStatus =c(Methylated =  "darkslateblue", Unmethylated = "darkslategray3")
)

pheatmap(tissue56.sig, col = bluered(100), show_rownames=T, cluster_cols=F, cluster_rows=F, scale="row",
         breaks= seq(-2,2,length.out = 101), fontsize = 10,fontsize_col = 4,labels_row = "",
         annotation_col = annotation_col, 
         annotation_colors=ann_colors,
         annotation_names_col = TRUE)

png("Tissue_heatmap.png", width = 7, height = 7, units = 'in', res = 600)
pheatmap(tissue56.sig, col = bluered(100), show_rownames=T, cluster_cols=F, cluster_rows=F, scale="row",
         breaks= seq(-2,2,length.out = 101), fontsize = 10,fontsize_col = 4,labels_row = "",
         annotation_col = annotation_col, 
         annotation_colors=ann_colors,
         annotation_names_col = TRUE)
dev.off()
