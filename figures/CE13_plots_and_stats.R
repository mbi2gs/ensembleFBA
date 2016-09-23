# Stats and plot generation for computational experiment 7
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

drugDist_df = data.frame(numDrugs=c(113,15,34,11,44,44),numSpecies=c(1,2,3,4,5,6))

p1 = ggplot(drugDist_df,aes(x=numSpecies,y=numDrugs)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("Number of Species") +
  ylab("Num. Small Molecules") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE13_distribution_small_molecules.tiff",width = 7, height = 5, units = "cm", dpi = 600)

#---------------------------------------------------------------------
# Calculate p-values and display subsystem enrichment in a heat map
#---------------------------------------------------------------------
file_name = "..\\data\\EssRxns_S.mitis_enrichment.tsv"
mitis_params = read.table(file_name, row.names=NULL, sep = "\t", header=F)
pvals = matrix(0,nrow = nrow(mitis_params)-1, ncol = 6)
N = mitis_params[1,2]
M = mitis_params[1,1] - N
mitis_params = mitis_params[-1,]

# Parameter definitions are from http://www.mathworks.com/help/stats/hygepdf.html
for(i in 1:nrow(mitis_params))
{
  k = mitis_params[i,1]
  x = mitis_params[i,2]
  pvals[i,1] = sum(dhyper(x:k,N,M,k))
}

file_name = "..\\data\\EssRxns_S.gallolyticus_enrichment.tsv"
gallolyticus_params = read.table(file_name, row.names=NULL, sep = "\t", header=F)
N = gallolyticus_params[1,2]
M = gallolyticus_params[1,1] - N
gallolyticus_params = gallolyticus_params[-1,]

# Parameter definitions are from http://www.mathworks.com/help/stats/hygepdf.html
for(i in 1:nrow(gallolyticus_params))
{
  k = gallolyticus_params[i,1]
  x = gallolyticus_params[i,2]
  pvals[i,2] = sum(dhyper(x:k,N,M,k))
}

file_name = "..\\data\\EssRxns_S.oralis_enrichment.tsv"
oralis_params = read.table(file_name, row.names=NULL, sep = "\t", header=F)
N = oralis_params[1,2]
M = oralis_params[1,1] - N
oralis_params = oralis_params[-1,]

# Parameter definitions are from http://www.mathworks.com/help/stats/hygepdf.html
for(i in 1:nrow(oralis_params))
{
  k = oralis_params[i,1]
  x = oralis_params[i,2]
  pvals[i,3] = sum(dhyper(x:k,N,M,k))
}

file_name = "..\\data\\EssRxns_S.equinus_enrichment.tsv"
equinus_params = read.table(file_name, row.names=NULL, sep = "\t", header=F)
N = equinus_params[1,2]
M = equinus_params[1,1] - N
equinus_params = equinus_params[-1,]

# Parameter definitions are from http://www.mathworks.com/help/stats/hygepdf.html
for(i in 1:nrow(equinus_params))
{
  k = equinus_params[i,1]
  x = equinus_params[i,2]
  pvals[i,4] = sum(dhyper(x:k,N,M,k))
}

file_name = "..\\data\\EssRxns_S.pneumoniae_enrichment.tsv"
pneumoniae_params = read.table(file_name, row.names=NULL, sep = "\t", header=F)
N = pneumoniae_params[1,2]
M = pneumoniae_params[1,1] - N
pneumoniae_params = pneumoniae_params[-1,]

# Parameter definitions are from http://www.mathworks.com/help/stats/hygepdf.html
for(i in 1:nrow(pneumoniae_params))
{
  k = pneumoniae_params[i,1]
  x = pneumoniae_params[i,2]
  pvals[i,5] = sum(dhyper(x:k,N,M,k))
}

file_name = "..\\data\\EssRxns_S.vestibularis_enrichment.tsv"
vestibularis_params = read.table(file_name, row.names=NULL, sep = "\t", header=F)
N = vestibularis_params[1,2]
M = vestibularis_params[1,1] - N
vestibularis_params = vestibularis_params[-1,]

# Parameter definitions are from http://www.mathworks.com/help/stats/hygepdf.html
for(i in 1:nrow(vestibularis_params))
{
  k = vestibularis_params[i,1]
  x = vestibularis_params[i,2]
  pvals[i,6] = sum(dhyper(x:k,N,M,k))
}

# Read in subsystem names
file_name = "..\\data\\CE13_uniqueSubsystems.txt"
subsystems = read.table(file_name, row.names=NULL, sep = "\t", header=F)
row.names(pvals) = subsystems$V1
colnames(pvals) = c("S.mitis","S.gallolyticus","S.oralis","S.equinus","S.pneumoniae","S.vestibularis")

#---------------------------------------------------------------------
# Draw Heat Map
#---------------------------------------------------------------------
library(pheatmap)
negLogPvals = -log(pvals)
negLogPvals = negLogPvals[order(-rowSums(negLogPvals)),]
negLogPvals = negLogPvals[rowSums(negLogPvals) > 1.5,]

write.table(negLogPvals, file = "..\\data\\CE13_subsystemEnrichment.tsv", append = FALSE, quote = FALSE, sep = "\t",row.names = TRUE,col.names = TRUE)

colors = colorRampPalette(c("white", 'blue'))(50)
fontsize = 10

hm.parameters <- list(negLogPvals, # breaks=bk2,
                      color = colors,
                      cellwidth = 15, cellheight = 15, scale = "none",
                      treeheight_row = 20,
                      treeheight_col = 0,
                      fontsize = fontsize, fontsize_row = fontsize,
                      fontsize_col = fontsize,
                      kmeans_k = NA,
                      Rowv=FALSE,na.rm=F, na.color="black",  Colv=FALSE, 
                      show_rownames = T, show_colnames = T,
                      main = NA,
                      cluster_rows = TRUE, cluster_cols = FALSE,
                      filename = "CE13_subsys_enrichment.tiff")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)


# Correlation analysis
correlations = matrix(0,nrow = nrow(negLogPvals), ncol = 6)
rownames(correlations) = rownames(negLogPvals)
colnames(correlations) = colnames(negLogPvals)
  
for(i in 1:nrow(negLogPvals))
{
  mitis = c(1,0,0,0,0,0)
  gallolyticus = c(0,1,0,0,0,0)
  oralis = c(0,0,1,0,0,0)
  equinus = c(0,0,0,1,0,0)
  pneumoniae = c(0,0,0,0,1,0)
  vestibularis = c(0,0,0,0,0,1)
  
  correlations[i,1] = cor(mitis,negLogPvals[i,])
  correlations[i,2] = cor(gallolyticus,negLogPvals[i,])
  correlations[i,3] = cor(oralis,negLogPvals[i,])
  correlations[i,4] = cor(equinus,negLogPvals[i,])
  correlations[i,5] = cor(pneumoniae,negLogPvals[i,])
  correlations[i,6] = cor(vestibularis,negLogPvals[i,])
}

write.table(correlations, file = "..\\data\\CE13_subsystemEnrichment_correlations.tsv", append = FALSE, quote = FALSE, sep = "\t",row.names = TRUE,col.names = TRUE)



