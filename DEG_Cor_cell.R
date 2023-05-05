if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")

library(readxl)
library(edgeR)
library(survival)
library(openxlsx)
library(dplyr)
library(DESeq2)

##### tissue #####
outcomes.c = read_excel("Data/Cell/Metadata_cellculture.xlsx", 1)
sampleNames.c = outcomes.c$GS.number
sampleNames.c = gsub("GS", "GS.", sampleNames.c, fixed = TRUE)

# Define the responders and non-responders group / make list of corresponding sample names
responders <- outcomes.c[outcomes.c$`Response group`=="Responder", c("GS.number", "% cell viability_100uM TMZ")]
nonresponders <- outcomes.c[outcomes.c$`Response group`=="Non-Responder", c("GS.number", "% cell viability_100uM TMZ")]

sampleNames.R = responders$GS.number
sampleNames.R = gsub("GS", "GS.", sampleNames.R, fixed = TRUE)
sampleNames.N = nonresponders$GS.number
sampleNames.N = gsub("GS", "GS.", sampleNames.N, fixed = TRUE)

# Create dataframe with each samplename corresponding to the right group
df = data.frame()
df <- rbind(df, as.data.frame(sampleNames.c))
df$condition <- outcomes.c$`Response group`
rownames(df) <- df$sampleNames.c
df$sampleNames.c = NULL

#### read data ####
#sampleNames.c = rownames(df)
cell.count = read.table("Data/Cell/count_data_24_cells.txt",sep = "\t", header = TRUE, row.names = 1)
cell.count19 = cell.count[sampleNames.c]

# Set all counts of 4 or lower to 0
# Replace counts less than or equal to 4 with 0 for all samples
samples <- colnames(cell.count19)  # Extract sample column names
for (sample in samples) {
  cell.count19[[sample]] <- ifelse(cell.count19[[sample]] <= 4, 0, cell.count19[[sample]])
}

# Filter step, filter rows which have less than 5% counts in both groups (responders or nonresponders) --> keep rows with more than 10% counts in 1 of the groups 
# Calculate the percentage of non-zero samples separately for the first two columns and the last two columns
cell.count19$non_zero_precentage_responders <- rowMeans(cell.count19[, sampleNames.R] > 0) * 100
cell.count19$non_zero_precentage_nonresponders <- rowMeans(cell.count19[, sampleNames.N] > 0) * 100

# Set the threshold value
threshold <- 10

# Filter the rows where the non-zero percentage is higher than the threshold in at least one of the sets of columns
filtered_df <- cell.count19 %>% filter(non_zero_precentage_responders > threshold | non_zero_precentage_nonresponders > threshold)
filtered_df$non_zero_precentage_responders = NULL
filtered_df$non_zero_precentage_nonresponders = NULL

# Check whether the samplenames are in the same order
all(rownames(df) == colnames(filtered_df))

'MGMT' %in% rownames(filtered_df)
'TERT' %in% rownames(filtered_df)

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_df,
                              colData = df,
                              design = ~ condition)


# Define the reference level 
dds$condition <- factor(dds$condition, levels = c("Responder","Non-Responder"))

# DEG analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Responder","Non-Responder"))

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

# Save analysis
write.csv(as.data.frame(DE_genes_005),
          file="Result_cell/DE_genes_cell_005.csv")

#####################################################################################################
#####################################################################################################
#####################################################################################################

##### cell #####

outcomes.c = read_excel("Data/Cell/GLIOTRAIN samples including MGMT status_IN.xlsx", 1)
outcomes.c = outcomes.c[order(outcomes.c$`% cell viability_100uM TMZ`),]

sampleNames.c = outcomes.c$GS.number
sampleNames.c = gsub("GS", "GS.", sampleNames.c, fixed = TRUE)


#### read data ####
cell.count = read.table("Data/Cell/count_data_24_cells.txt",sep = "\t", header = TRUE, row.names = 1)
cell.count19 = cell.count[sampleNames.c]

gene_list <- scan('Result_cell/DE_genes_cell_005.txt', what="character")
cell.count19 <- cell.count19[rownames(cell.count19) %in% gene_list,]

# Set all counts of 4 or lower to 0
# Replace counts less than or equal to 4 with 0 for all samples
samples <- colnames(cell.count19)  # Extract sample column names
for (sample in samples) {
  cell.count19[[sample]] <- ifelse(cell.count19[[sample]] <= 4, 0, cell.count19[[sample]])
}

cell19.obj = DGEList(counts=as.matrix(cell.count19), genes=rownames(cell.count19))

# norm
norm.factor19 = calcNormFactors(cell19.obj, method = "TMM")
cell.cpm19 = cpm(norm.factor19, log=TRUE)
cell.cpm19 = as.data.frame(cell.cpm19)

cor_analysis <- function(data,outcomes){
  
  data.t = as.data.frame(t(data))
  
  cor.coe= sapply(data.t, function(x) cor.test(x, outcomes, method ="spearman",exact = FALSE)$estimate)
  cor.p =  sapply(data.t, function(x) cor.test(x, outcomes, method ="spearman",exact = FALSE)$p.value)
  cor.all = cbind("coefficient"=cor.coe, "p"=cor.p)
  cor.all = as.data.frame(cor.all)
  
  # adjust the p-value
  cor.all$p.adj = p.adjust(cor.all$p, method = "BH")
  
  # select genes with p-value lower than 0.05
  cor.all$gene = rownames(cor.all)
  cor.all$gene = gsub(".rho","",cor.all$gene)
  
  # order by coefficient
  cor.all = cor.all[order(cor.all$coefficient),]
  
  return(cor.all)
}

cell.gene19 = cor_analysis(cell.cpm19, outcomes.c$`% cell viability_100uM TMZ`)
cell.gene19.005 = cell.gene19[cell.gene19$p < 0.05,]
cell.gene19.001 = cell.gene19[cell.gene19$p < 0.01,]

write.csv(as.data.frame(cell.gene19.005),
          file="Result_cell/Correlation_genes_005.csv")

####################################################################
####################################################################
####################################################################

library(pheatmap)
library(gplots)

cell19.sig = cell.cpm19[cell.gene19.005$gene,]

outcomes.c19 = outcomes.c

cell.gene19.0005 = rbind(head(cell.gene19,50),tail(cell.gene19,50))
cell19.sig = cell.cpm19[cell.gene19.0005$gene,]

grp = rep("Non-Responders",19)
grp[outcomes.c19$`% cell viability_100uM TMZ`<75] = "Intermediates"
grp[outcomes.c19$`% cell viability_100uM TMZ`<50] = "Responders"
grp = factor(grp,levels = c("Non-Responders","Intermediates","Responders"))

annotation_col  = data.frame(group = grp,
                             Viability100uM= outcomes.c19$`% cell viability_100uM TMZ`,
                             MGMTcultureStatus = outcomes.c19$MGMT_status_culture
)
rownames(annotation_col) = colnames(cell19.sig)

ann_colors = list(group = c('Non-Responders'="brown4",Intermediates="deepskyblue2",Responders="darkseagreen"),
                  Viability100uM= c("darkseagreen","chocolate2","deeppink4"),
                  MGMTcultureStatus =c(Methylated =  "darkslateblue", Unmethylated = "darkslategray3")
                  
)
png("Cell_culture_heatmap.png", width = 7, height = 7, units = 'in', res = 600)
pheatmap(cell19.sig, col = bluered(100), show_rownames=T, cluster_cols=F, cluster_rows=F, scale="row",
         breaks= seq(-2,2,length.out = 101), fontsize = 8,fontsize_row = 4,labels_col = "",#labels_row = "",
         annotation_col = annotation_col, 
         annotation_colors=ann_colors,
         #cex = 0.9,
         cellwidth=15, cellheight=4,
         annotation_names_col = TRUE)
dev.off()

