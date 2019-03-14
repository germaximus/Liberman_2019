library(magrittr)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(rstudioapi)
library(dplyr)
library(stringr)
library(pheatmap)
library(data.table)
library(dendextend)
setwd(dirname(getActiveDocumentContext()$path))


#------------------------------------------- Import data ----------------------------------------------------
SampleTable <- fread(file = "./geneCounts/SampleTable.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# import mRNA-seq data
mrnaCounts <- lapply(SampleTable[type == 'mRNAseq', file], function(x) {
                data <- read.table(file = paste0("./geneCounts/", x), skip=4, row.names = 1, stringsAsFactors = FALSE) %>% select(1)
}) %>% as.data.frame() %>% setNames(paste0('id', SampleTable$sampleID[1:12]))

# import Ribo-seq data
riboCounts <- lapply(SampleTable[type == 'riboSeq', file], function(x) {
  data <- read.table(file = paste0("./geneCounts/", x), skip=4, row.names = 1, stringsAsFactors = FALSE) %>% select(1)
  
}) %>% as.data.frame() %>% setNames(paste0('id', SampleTable$sampleID[1:12]))


#------------------------------------------ Matrix Scatterplots -----------------------------------------------------
png('MatrixPlot_mRNAseq.png', width = 2400, height = 2400, units = 'px', res = 300)
pairs(mrnaCounts, log = "yx", pch = '.', gap = 0.25, xaxt = "n", yaxt = "n", labels = sapply(SampleTable[type == 'mRNAseq', file], function(x) { str_extract(x, '[^_]+') }))
dev.off()

png('MatrixPlot_Riboseq.png', width = 2400, height = 2400, units = 'px', res = 300)
pairs(riboCounts, log = "yx", pch = '.', gap = 0.25, xaxt = "n", yaxt = "n", labels = sapply(SampleTable[type == 'riboSeq', file], function(x) { str_extract(x, '[^_]+')    }))
dev.off()

#-------------------------------------- PCA plots --------------------------------------------
PlotData <- mrnaCounts[rowMeans(mrnaCounts) >= 10, ] %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         group = paste0(SampleTable$strain[1:12], '-', SampleTable$temperature[1:12]),
                         name  = SampleTable[type == 'mRNAseq', file]
)

png(file = "PCA_mRNAseq.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
  theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
  theme(axis.title = element_text(size = rel(1.5)))+
  theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()


PlotData <- riboCounts[rowMeans(riboCounts) >= 10, ] %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         group = paste0(SampleTable$strain[1:12], '-', SampleTable$temperature[1:12]),
                         name  = SampleTable[type == 'riboSeq', file]
)

png(file = "PCA_riboSeq.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
  theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
  theme(axis.title = element_text(size = rel(1.5)))+
  theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()

#----------------------------------- Correlation Heatmaps ---------------------------------------
colors    <- colorRampPalette(brewer.pal(9, "GnBu"))(100) %>% rev()

PlotData  <- mrnaCounts[rowMeans(mrnaCounts) >= 10, ] %>% as.matrix() %>% rlog()
distance  <- t(PlotData) %>% dist()
attributes(distance)$Labels <- sapply(SampleTable[type == 'mRNAseq', file], function(x) { str_extract(x, '[^_]+')  }) %>% unname()
mat <- as.matrix(distance)
p_mrna <- pheatmap(mat, clustering_distance_rows = distance, clustering_distance_cols = distance, color = colors, border_color = 'white', fontsize = 14, silent = TRUE)

png("CorrelationHeatmap_mRNAseq.png", width = 2400, height = 2400, res = 300, unit = 'px')
p_mrna
dev.off()

PlotData  <- riboCounts[rowMeans(riboCounts) >= 10, ] %>% as.matrix() %>% rlog()
distance  <- t(PlotData) %>% dist()
attributes(distance)$Labels <- sapply(SampleTable[type == 'riboSeq', file], function(x) { str_extract(x, '[^_]+')  }) %>% unname()
mat <- as.matrix(distance)
p_ribo <- pheatmap(mat, clustering_distance_rows = distance, clustering_distance_cols = distance, color = colors, border_color = 'white', fontsize = 14, silent = TRUE)
p_ribo_ordered <- pheatmap(mat, cluster_rows = rotate(p_ribo$tree_row, p_mrna$tree_row$labels[p_mrna$tree_row$order]), cluster_col = rotate(p_ribo$tree_col, p_mrna$tree_col$labels[p_mrna$tree_col$order]) , color = colors, border_color = 'white', fontsize = 14, silent = TRUE)

png("CorrelationHeatmap_riboSeq.png", width = 2400, height = 2400, res = 300, unit = 'px')
p_ribo_ordered
dev.off()


#--------------------------- Heatmap genes vs samples --------------------------------------------
colors    <- colorRampPalette(c("blue","black","yellow"))(256)

PlotData  <- mrnaCounts[rowMeans(mrnaCounts) >= 10, ] %>% as.matrix() %>% rlog()
colnames(PlotData) <- sapply(SampleTable[type == 'mRNAseq', file], function(x) { str_extract(x, '[^_]+')  }) %>% unname()
distCol <- t(PlotData) %>% dist()
hclustCol <- hclust(distCol, method = 'complete')
annotationCol <- data.frame(Temperature = c(rep('20(\u00B0C)', 6), rep('37(\u00B0C)', 6)))
row.names(annotationCol) <- colnames(PlotData)[hclustCol$order]


p_mrna <- pheatmap(PlotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = hclustCol,
               color = colors,  border_color = NA, show_rownames = FALSE, fontsize = 14, silent = TRUE,
               annotation_col = annotationCol, annotation_names_col = FALSE,
)

png("GenesHeatmap_mRNAseq.png", width = 2400, height = 2400, res = 300, unit = 'px')
p_mrna
dev.off()


PlotData  <- riboCounts[rowMeans(riboCounts) >= 10, ] %>% as.matrix() %>% rlog()
colnames(PlotData) <- sapply(SampleTable[type == 'riboSeq', file], function(x) { str_extract(x, '[^_]+')  }) %>% unname()
distCol <- t(PlotData) %>% dist()
hclustCol <- hclust(distCol, method = 'complete')
annotationCol <- data.frame(Temperature = c(rep('20(\u00B0C)', 6), rep('37(\u00B0C)', 6)))
row.names(annotationCol) <- colnames(PlotData)[hclustCol$order]


p_ribo <- pheatmap(PlotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = hclustCol,
                   color = colors,  border_color = NA, show_rownames = FALSE, fontsize = 14, silent = TRUE,
                   annotation_col = annotationCol, annotation_names_col = FALSE,
)

p_ribo_ordered <- pheatmap(PlotData, clustering_distance_rows = 'correlation', clustering_method = "average", scale = 'row', cluster_cols = rotate(p_ribo$tree_col, p_mrna$tree_col$labels[p_mrna$tree_col$order]),
                   color = colors,  border_color = NA, show_rownames = FALSE, fontsize = 14, silent = TRUE,
                   annotation_col = annotationCol, annotation_names_col = FALSE,
                
)

png("GenesHeatmap_riboSeq.png", width = 2400, height = 2400, res = 300, unit = 'px')
p_ribo_ordered
dev.off()



###############################################################################################################################
###################################   Differential Gene Expression Analysis   #################################################
###############################################################################################################################

colData <- data.frame(sample = colnames(mrnaCounts),
                      strain = factor(c(rep("WT", 6), rep("CD", 6)), levels = c("WT","CD")),
                      temperature = factor(c(rep('20', 3), rep('37', 3), rep('20', 3), rep('37', 3)), levels = c("20","37")),
                      stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(countData = mrnaCounts, colData = colData, design =~ strain + temperature + temperature:strain)
dds <- dds[rowMeans(counts(dds)) >= 10, ]
dds <- DESeq(dds)


resultsNames(dds)

result1 <-  results(dds, name = 'temperature_37_vs_20', cooksCutoff = FALSE, independentFiltering = FALSE)
result2 <-  results(dds, name = 'strain_CD_vs_WT', cooksCutoff = FALSE, independentFiltering = FALSE)
result3 <-  results(dds, name = 'strainCD.temperature37', cooksCutoff = FALSE, independentFiltering = FALSE)
summary(result1)
summary(result2)
summary(result3)

result1 <- result1[order(result1$padj),]
result2 <- result2[order(result2$padj),]
result3 <- result3[order(result3$padj),]



write.csv(result1, file="temperature_37_vs_20.csv")
write.csv(result2, file="strain_CD_vs_WT.csv")








