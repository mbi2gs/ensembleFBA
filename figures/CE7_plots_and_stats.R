# Stats and plot generation for computational experiment 7
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE7_networkPrecisionsNRecalls.tsv"
networkAPR = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(networkAPR) = c("Accuracy","Precision","Recall")

file_name = "..\\data\\CE7_ensemblePrecisionsNRecalls.tsv"
ensembleAPR = read.table(file_name, row.names=NULL, sep = "\t", header=F)
rownames(ensembleAPR) = c("t1","tHalf","tAll")
ensembleAPR = data.frame(ensembleAPR)
ensembleAPR = cbind(ensembleAPR,c("t1","tHalf","tAll"))
colnames(ensembleAPR) = c("Accuracy","Precision","Recall","Threshold")
ensembleAPR_df = melt(ensembleAPR)
ensembleAPR_df$Threshold = factor(ensembleAPR_df$Threshold, levels = c("t1","tHalf","tAll"))

# Summarize network data
tmpARange = c(mean(networkAPR$Accuracy),max(networkAPR$Accuracy),min(networkAPR$Accuracy))
tmpPRange = c(mean(networkAPR$Precision),max(networkAPR$Precision),min(networkAPR$Precision))
tmpRRange = c(mean(networkAPR$Recall),max(networkAPR$Recall),min(networkAPR$Recall))
networkSummaries = rbind(tmpARange,tmpPRange,tmpRRange)
networkSummaries = data.frame(networkSummaries)
networkSummaries = cbind(networkSummaries,c("Accuracy","Precision","Recall"))
rownames(networkSummaries) = c("Accuracy","Precision","Recall")
colnames(networkSummaries) = c("nMean","nMax","nMin","Measure")

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Precision)
#-----------------------------------------------------------------------
p1 = ggplot(networkSummaries,aes(x=factor(Measure),y=nMean)) +
  geom_pointrange(size=.7,aes(ymin=nMax, ymax=nMin)) +
  theme_bw() +
  geom_point(data=ensembleAPR_df,size=1.3,aes(x=factor(variable),y=value,color=Threshold)) +
  geom_line(data=ensembleAPR_df,size=1.2,aes(x=factor(variable),y=value,color=Threshold,group=Threshold)) +
  xlab("Measure") +
  ylab("") +
  ylim(0,1) +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE7_ensemble_gene_prediction_performance.tiff",width = 15, height = 6, units = "cm", dpi = 600)

