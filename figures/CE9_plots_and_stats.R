# Stats and plot generation for computational experiment 7
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
# "Any"
file_name = "..\\data\\CE9_ensembleAccuracyBySize_any.tsv"
ensembleAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
sizes = ensembleAccuracies[,1]
ensembleAccuracies = ensembleAccuracies[,-1]
summaryEAs_any = data.frame(cbind(sizes,apply(ensembleAccuracies,1,mean),apply(ensembleAccuracies,1,sd)))
colnames(summaryEAs_any) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensemblePrecisionBySize_any.tsv"
ensemblePrecisions = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensemblePrecisions = ensemblePrecisions[,-1]
summaryEPs_any = data.frame(cbind(sizes,apply(ensemblePrecisions,1,mean),apply(ensemblePrecisions,1,sd)))
colnames(summaryEPs_any) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensembleRecallBySize_any.tsv"
ensembleRecalls = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensembleRecalls = ensembleRecalls[,-1]
summaryERs_any = data.frame(cbind(sizes,apply(ensembleRecalls,1,mean),apply(ensembleRecalls,1,sd)))
colnames(summaryERs_any) = c("Sizes","Mean","StDev")

# "Majority"
file_name = "..\\data\\CE9_ensembleAccuracyBySize_majority.tsv"
ensembleAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
sizes = ensembleAccuracies[,1]
ensembleAccuracies = ensembleAccuracies[,-1]
summaryEAs_majority = data.frame(cbind(sizes,apply(ensembleAccuracies,1,mean),apply(ensembleAccuracies,1,sd)))
colnames(summaryEAs_majority) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensemblePrecisionBySize_majority.tsv"
ensemblePrecisions = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensemblePrecisions = ensemblePrecisions[,-1]
summaryEPs_majority = data.frame(cbind(sizes,apply(ensemblePrecisions,1,mean),apply(ensemblePrecisions,1,sd)))
colnames(summaryEPs_majority) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensembleRecallBySize_majority.tsv"
ensembleRecalls = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensembleRecalls = ensembleRecalls[,-1]
summaryERs_majority = data.frame(cbind(sizes,apply(ensembleRecalls,1,mean),apply(ensembleRecalls,1,sd)))
colnames(summaryERs_majority) = c("Sizes","Mean","StDev")

# "Consensus"
file_name = "..\\data\\CE9_ensembleAccuracyBySize_consensus.tsv"
ensembleAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
sizes = ensembleAccuracies[,1]
ensembleAccuracies = ensembleAccuracies[,-1]
summaryEAs_consensus = data.frame(cbind(sizes,apply(ensembleAccuracies,1,mean),apply(ensembleAccuracies,1,sd)))
colnames(summaryEAs_consensus) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensemblePrecisionBySize_consensus.tsv"
ensemblePrecisions = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensemblePrecisions = ensemblePrecisions[,-1]
summaryEPs_consensus = data.frame(cbind(sizes,apply(ensemblePrecisions,1,mean),apply(ensemblePrecisions,1,sd)))
colnames(summaryEPs_consensus) = c("Sizes","Mean","StDev")

file_name = "..\\data\\CE9_ensembleRecallBySize_consensus.tsv"
ensembleRecalls = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensembleRecalls = ensembleRecalls[,-1]
summaryERs_consensus = data.frame(cbind(sizes,apply(ensembleRecalls,1,mean),apply(ensembleRecalls,1,sd)))
colnames(summaryERs_consensus) = c("Sizes","Mean","StDev")

#-----------------------------------------------------------------------
# Plot accuracy of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
p1 = ggplot(summaryEAs_any,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Accuracy") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_accuracy_by_size_any.tiff",width = 8, height = 7, units = "cm", dpi = 600)

p1 = ggplot(summaryEAs_majority,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Accuracy") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_accuracy_by_size_majority.tiff",width = 8, height = 7, units = "cm", dpi = 600)


p1 = ggplot(summaryEAs_consensus,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Accuracy") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_accuracy_by_size_consensus.tiff",width = 8, height = 7, units = "cm", dpi = 600)


#-----------------------------------------------------------------------
# Plot precision of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
p1 = ggplot(summaryEPs_any,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Precision") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_precision_by_size_any.tiff",width = 8, height = 7, units = "cm", dpi = 600)

p1 = ggplot(summaryEPs_majority,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Precision") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_precision_by_size_majority.tiff",width = 8, height = 7, units = "cm", dpi = 600)

p1 = ggplot(summaryEPs_consensus,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Precision") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_precision_by_size_consensus.tiff",width = 8, height = 7, units = "cm", dpi = 600)


#-----------------------------------------------------------------------
# Plot recall of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
p1 = ggplot(summaryERs_any,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Recall") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_recall_by_size_any.tiff",width = 8, height = 7, units = "cm", dpi = 600)

p1 = ggplot(summaryERs_majority,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Recall") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_recall_by_size_majority.tiff",width = 8, height = 7, units = "cm", dpi = 600)

p1 = ggplot(summaryERs_consensus,aes(x=Sizes,y=Mean)) +
  geom_errorbar(size=0.2,aes(ymax=Mean+StDev,ymin=Mean-StDev)) +
  geom_line(size=2) +
  geom_point(size=2) +
  theme_bw() +
  xlab("Ensemble Size") +
  ylab("Recall") +
  theme(text = element_text(size=12))
print(p1)

ggsave("CE9_ensemble_recall_by_size_consensus.tiff",width = 8, height = 7, units = "cm", dpi = 600)

