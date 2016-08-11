# Stats and plot generation for computational experiment 7
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE9_ensembleAccuracyBySize.tsv"
ensembleAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
sizes = ensembleAccuracies[,1]
ensembleAccuracies = ensembleAccuracies[,-1]
summaryEAs = data.frame(cbind(sizes,apply(ensembleAccuracies,1,mean),apply(ensembleAccuracies,1,sd)))
colnames(summaryEAs) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensemblePrecisionBySize.tsv"
ensemblePrecisions = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensemblePrecisions = ensemblePrecisions[,-1]
summaryEPs = data.frame(cbind(sizes,apply(ensemblePrecisions,1,mean),apply(ensemblePrecisions,1,sd)))
colnames(summaryEPs) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensembleRecallBySize.tsv"
ensembleRecalls = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensembleRecalls = ensembleRecalls[,-1]
summaryERs = data.frame(cbind(sizes,apply(ensembleRecalls,1,mean),apply(ensembleRecalls,1,sd)))
colnames(summaryERs) = c("Sizes","Mean","StDev")

#-----------------------------------------------------------------------
# Plot accuracy of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
p1 = ggplot(summaryEAs,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Accuracy") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_accuracy_by_size.tiff",width = 8, height = 10, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot precision of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
p1 = ggplot(summaryEPs,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Precision") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_precision_by_size.tiff",width = 8, height = 10, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot recall of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
p1 = ggplot(summaryERs,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Recall") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_recall_by_size.tiff",width = 8, height = 10, units = "cm", dpi = 600)

