# script to perform differential gene expression analysis using DESeq2 package
setwd("C:/Users/qinqin.zha/OneDrive - NeoGenomics Laboratories, Inc/Documents/04_Literature_AACRposters/Ohio_BreastCancer_Coorlab_DSP/Mar042023/")


library(dplyr)
# load libraries
library(DESeq2)
library(tidyverse)
library("readxl")

# Step 1: preparing count data ----------------

# read in counts data
cts <- read.csv('DSP-NGSQ3_MatrixCount_10232022.csv', header = TRUE, row.names = 1)
head(cts)
str(cts)
colData1 <- read.csv("MetaData_colomn_CRC_health.csv", header = TRUE, row.names = 1)
str(colData1)
df1 <- cts %>%
     select(row.names(colData1))

str(df1)
# making sure the row names in colData matches to column names in counts_data
all(colnames(df1) %in% rownames(colData1))

# are they in the same order?
all(colnames(df1) == rownames(colData1))
#save to csv file
write.csv(df1, file = "DSP_countCRC_H.csv")
cts1 <- read.csv("AvergaedCounts1.csv", header = TRUE, row.names=1)
cts1
colData2 <- read.csv("MetaData_re_CRC_H1.csv", header = TRUE, row.names = 1)
cts2 <- cts1

all(colnames(cts2) %in% rownames(colData2))

# are they in the same order?
all(colnames(cts2) == rownames(colData2))

# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = round(cts2),
                              colData = colData2,
                              design = ~ Indication)

dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds
dds$Indication <- relevel(dds$Indication, ref = "Healthy")

dds <- DESeq(dds)
res <- results(dds)

res
summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

res <- results(dds, contrast = c("Indication", "CRC", "Healthy"))

vsdata <- vst(dds, blind = FALSE)

plotPCA(vsdata, intgroup = "Indication" )


sigs <- na.omit(res)

#filtere by padj value, from 1673  to 1028 genes
sigs <- sigs[sigs$padj < 0.05,]
sigs

write.csv(sigs, file = "deseq_results_sig.csv")


df <- as.data.frame(sigs)

df.top <- df[ ( df$baseMean > 100) & (abs(df$log2FoldChange) > 0.5) ,]
df.top
str(df.top)

df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE), ]
df.top

rlog_out <- rlog(dds, blind=FALSE)  #get normalized count data from dds object
mat <-assay(rlog_out)[rownames(df.top), rownames(colData2)]
colnames(mat) <- rownames(colData2)

base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale))

colnames(mat.scaled) <- colnames(mat)
mat.scaled  # people use this mat.scaled for heatmap

num_keep <- 25

#1 to num_keep  len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))
rows_keep


l2_val <- as.matrix((df.top[rows_keep,]$log2FoldChange))
colnames(l2_val) <- "logFC"     #adding column name for log2 value matrix

mean <- as.matrix((df.top[rows_keep,]$baseMean))
colnames(mean) <- "AveExpr"   # adding columname for mean matrix


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

library(RColorBrewer)

library(circlize)

#maps value between b/w/r for min and max 12 values

col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red"))

# maps between 0% quantile, and 75% quantile of mean values  ---0, 25, 50, 75, 100 
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))


ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F,
              column_labels = colnames(mat.scaled), name = "Z-score",
              cluster_columns = T)

h2 <- Heatmap(l2_val, row_labels = rownames(df.top)[rows_keep],
              cluster_rows =F, name ="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col)
              {#add text to each grid 
                grid.text(round(l2_val[i, j], 2), x, y)}              )



h3 <- Heatmap(mean, row_labels = rownames(df.top)[rows_keep],
              cluster_rows =F, name ="AveExpr",  col = col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col)
              {#add text to each grid 
                grid.text(round(mean[i, j], 2), x, y)}  ) 



h <- h1+h2+h3
h


png("./heatmap_CRC_H3.png", res =400, width = 6000, height = 5500)

print(h)

dev.off()


# from mat.scaled
df.top
write.csv(df.top, file = "topDeseq_results_sig.csv")

BiocManager::install("clusterProfiler")
BiocManager::install("BiocManager")
BiocManager::install("AnnotationDbi")

#GO enrishment analysis, added SYMBOL column

genes_to_test <- rownames(df.top[df.top$log2FoldChange > 0.5,])

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)

#GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

df1.top <- as.data.frame(df.top)

as.data.frame(GO_results)
rownames(df1.top)
df1.top$SYMBOL <- rownames(df1.top)
str(df1.top)

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))

png("out1.png", res =250, width = 1200, height = 2500 )
print(fit)
dev.off()

fit


df1.top
write.csv(df1.top, file = "Deseq_results_SYMBOL.csv")

#Volcano plot


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

install.packages("textshaping")
library(EnhancedVolcano)


EnhancedVolcano(df1.top, x = "log2FoldChange", y = "padj", lab =df1.top$SYMBOL,
                pCutoff = 1e-4, FCcutoff = 1)


selected =c("COL1A1", "COL3A1", "SOX9", "NR4A1", "CXCL8", "MYC", "B2M", "CD74")

jpeg("EV2.jpg")


EnhancedVolcano(df1.top, x = "log2FoldChange", y = "padj", lab =df1.top$SYMBOL, title ="Differential Expression Analysis", subtitle = "CRC vs Healthy",
                pCutoff = 10e-3, FCcutoff = 1, pointSize = 3.0, labSize = 6.0, selectLab = selected )

dev.off()



