# Stats and plot generation for computational experiment 8
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE8_ensemblePrecisionsNRecalls.tsv"
ensemblesAPR = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensemblesAPR) = c("Accuracy","Precision","Recall")

# Reformat Accuracy/Precision/Recall information
ensemblesAPR = ensemblesAPR[rev(rownames(ensemblesAPR)),]
accuracyMat = cbind(ensemblesAPR[13:16,1],ensemblesAPR[9:12,1],ensemblesAPR[5:8,1],ensemblesAPR[1:4,1])
colnames(accuracyMat) = c("100% GCs","80% GCs","60% GCs","30% GCs")
rownames(accuracyMat) = c("30% Genes","60% Genes","80% Genes","100% Genes")
precisionMat = cbind(ensemblesAPR[13:16,2],ensemblesAPR[9:12,2],ensemblesAPR[5:8,2],ensemblesAPR[1:4,2])
colnames(precisionMat) = c("100% GCs","80% GCs","60% GCs","30% GCs")
rownames(precisionMat) = c("30% Genes","60% Genes","80% Genes","100% Genes")
recallMat = cbind(ensemblesAPR[13:16,3],ensemblesAPR[9:12,3],ensemblesAPR[5:8,3],ensemblesAPR[1:4,3])
colnames(recallMat) = c("100% GCs","80% GCs","60% GCs","30% GCs")
rownames(recallMat) = c("30% Genes","60% Genes","80% Genes","100% Genes")

file_name = "..\\data\\CE8_ensembleJaccardSims.tsv"
ensemblesJS = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensemblesJS$rows = 1:16
colnames(ensemblesJS) = c("JaccardSimilarity","Rows")

# Reformat Accuracy/Precision/Recall information
ensemblesJS = ensemblesJS[16:1,]
jaccardSimMat = cbind(ensemblesJS[13:16,1],ensemblesJS[9:12,1],ensemblesJS[5:8,1],ensemblesJS[1:4,1])
colnames(jaccardSimMat) = c("100% GCs","80% GCs","60% GCs","30% GCs")
rownames(jaccardSimMat) = c("30% Genes","60% Genes","80% Genes","100% Genes")

#-----------------------------------------------------------------------
# Plot accuracy/precision/recall of various ensemble configurations
#-----------------------------------------------------------------------
library(pheatmap)
library(RColorBrewer)

col.pal = brewer.pal(9,"Blues")

hm.parameters <- list(accuracyMat, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 0,
                      treeheight_col = 0,
                      fontsize = 12, fontsize_row = 12,
                      fontsize_col = 12,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "CE8_accuracy_by_ensemble_diversity.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)

hm.parameters <- list(precisionMat, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 0,
                      treeheight_col = 0,
                      fontsize = 12, fontsize_row = 12,
                      fontsize_col = 12,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "CE8_precision_by_ensemble_diversity.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)

hm.parameters <- list(recallMat, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 0,
                      treeheight_col = 0,
                      fontsize = 12, fontsize_row = 12,
                      fontsize_col = 12,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "CE8_recall_by_ensemble_diversity.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)

hm.parameters <- list(jaccardSimMat, 
                      color = col.pal,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 0,
                      treeheight_col = 0,
                      fontsize = 12, fontsize_row = 12,
                      fontsize_col = 12,
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      main = "",
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = "CE8_jaccardSim_ensemble_diversity.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)