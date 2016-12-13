# Stats and plot generation for computational experiment 14
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE14_ensembleAccuracy.tsv"
ensembleAccuracy = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensembleAccuracy) = c("Accuracy_t1","Accuracy_tHalf","Accuracy_tAll","Fractions")
ensembleAccuracy_df = melt(ensembleAccuracy,id=c("Fractions"))

file_name = "..\\data\\CE14_ensemblePrecision.tsv"
ensemblePrecision = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensemblePrecision) = c("Precision_t1","Precision_tHalf","Precision_tAll","Fractions")
ensemblePrecision_df = melt(ensemblePrecision,id=c("Fractions"))

file_name = "..\\data\\CE14_ensembleRecall.tsv"
ensembleRecall = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensembleRecall) = c("Recall_t1","Recall_tHalf","Recall_tAll","Fractions")
ensembleRecall_df = melt(ensembleRecall,id=c("Fractions"))

#-----------------------------------------------------------------------
# Plot Accuracy
#-----------------------------------------------------------------------
p1 = ggplot(ensembleAccuracy_df,aes(x=Fractions,y=value,group=variable)) + 
  theme_bw() +
  geom_line(size=1,aes(color=variable,group=variable)) +
  geom_point(size=2.5,aes(color=variable,shape=variable)) +
  xlab("Fraction of Data") + 
  ylab("Accuracy") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE14_ensemble_accuracy.tiff",width = 15, height = 6, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot Precision
#-----------------------------------------------------------------------
p1 = ggplot(ensemblePrecision_df,aes(x=Fractions,y=value,group=variable)) + 
  theme_bw() +
  geom_line(size=1,aes(color=variable,group=variable)) +
  geom_point(size=2.5,aes(color=variable,shape=variable)) +
  xlab("Fraction of Data") + 
  ylab("Precision") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE14_ensemble_precision.tiff",width = 15, height = 6, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot Recall
#-----------------------------------------------------------------------
p1 = ggplot(ensembleRecall_df,aes(x=Fractions,y=value,group=variable)) + 
  theme_bw() +
  geom_line(size=1,aes(color=variable,group=variable)) +
  geom_point(size=2.5,aes(color=variable,shape=variable)) +
  xlab("Fraction of Data") + 
  ylab("Recall") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE14_ensemble_recall.tiff",width = 15, height = 6, units = "cm", dpi = 600)
