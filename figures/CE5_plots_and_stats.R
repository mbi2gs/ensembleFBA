# Stats and plot generation for computational experiment 5
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE5_networkPrecision.tsv"
networkPrecision = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(networkPrecision) = c("Precision","N_gcs")

file_name = "..\\data\\CE5_ensemblePrecision.tsv"
ensemblePrecision = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensemblePrecision) = c("Precision_t1","Precision_tHalf","Precision_tAll","N_gcs")
ensemblePrecision_df = melt(ensemblePrecision,id=c("N_gcs"))

file_name = "..\\data\\CE5_networkRecall.tsv"
networkRecall = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(networkRecall) = c("Recall","N_gcs")

file_name = "..\\data\\CE5_ensembleRecall.tsv"
ensembleRecall = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensembleRecall) = c("Recall_t1","Recall_tHalf","Recall_tAll","N_gcs")
ensembleRecall_df = melt(ensembleRecall,id=c("N_gcs"))

file_name = "..\\data\\CE5_networkAccuracy.tsv"
networkAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(networkAccuracies) = c("Accuracy","N_gcs")

file_name = "..\\data\\CE5_ensembleAccuracy.tsv"
ensembleAccuracy = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(ensembleAccuracy) = c("Accuracy_t1","Accuracy_tHalf","Accuracy_tAll","N_gcs")
ensembleAccuracy_df = melt(ensembleAccuracy,id=c("N_gcs"))

# Summarize network data
networkPrecisionSummaries = c()
for( i in c(2,5,10,15,20,25,30))
{
  data_i =  filter(networkPrecision, N_gcs == i)
  tmpRange = c(mean(data_i$Precision),max(data_i$Precision),min(data_i$Precision),i)
  networkPrecisionSummaries = rbind(networkPrecisionSummaries,tmpRange)
}
networkPrecisionSummaries = data.frame(networkPrecisionSummaries)
colnames(networkPrecisionSummaries) = c("nMean","nMax","nMin","N_gcs")

networkRecallSummaries = c()
for( i in c(2,5,10,15,20,25,30))
{
  data_i =  filter(networkRecall, N_gcs == i)
  tmpRange = c(mean(data_i$Recall),max(data_i$Recall),min(data_i$Recall),i)
  networkRecallSummaries = rbind(networkRecallSummaries,tmpRange)
}
networkRecallSummaries = data.frame(networkRecallSummaries)
colnames(networkRecallSummaries) = c("nMean","nMax","nMin","N_gcs")

networkAccuracySummaries = c()
for( i in c(2,5,10,15,20,25,30))
{
  data_i =  filter(networkAccuracies, N_gcs == i)
  tmpRange = c(mean(data_i$Accuracy),max(data_i$Accuracy),min(data_i$Accuracy),i)
  networkAccuracySummaries = rbind(networkAccuracySummaries,tmpRange)
}
networkAccuracySummaries = data.frame(networkAccuracySummaries)
colnames(networkAccuracySummaries) = c("nMean","nMax","nMin","N_gcs")

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Precision)
#-----------------------------------------------------------------------
p1 = ggplot(networkPrecisionSummaries,aes(x=N_gcs,y=nMean,group=N_gcs)) + 
  geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
  theme_bw() +
  geom_line(data=ensemblePrecision_df,size=1,aes(x=N_gcs,y=value,color=variable,group=variable)) +
  geom_point(data=ensemblePrecision_df,size=2.5,aes(x=N_gcs,y=value,color=variable,shape=variable)) +
  xlab("Num.Growth.Conditions") + 
  ylab("Precision") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE5_ensemble_precision.tiff",width = 15, height = 6, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Recall)
#-----------------------------------------------------------------------
p1 = ggplot(networkRecallSummaries,aes(x=N_gcs,y=nMean,group=N_gcs)) + 
  geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
  theme_bw() +
  geom_line(data=ensembleRecall_df,size=1,aes(x=N_gcs,y=value,color=variable,group=variable)) +
  geom_point(data=ensembleRecall_df,size=2.5,aes(x=N_gcs,y=value,color=variable,shape=variable)) +
  xlab("Num.Growth.Conditions") + 
  ylab("Recall") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE5_ensemble_recall.tiff",width = 15, height = 6, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Precision)
#-----------------------------------------------------------------------
p1 = ggplot(networkAccuracySummaries,aes(x=N_gcs,y=nMean,group=N_gcs)) + 
  geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
  theme_bw() +
  geom_line(data=ensembleAccuracy_df,size=1,aes(x=N_gcs,y=value,color=variable,group=variable)) +
  geom_point(data=ensembleAccuracy_df,size=2.5,aes(x=N_gcs,y=value,color=variable,shape=variable)) +
  xlab("Num.Growth.Conditions") + 
  ylab("Accuracy") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE5_ensemble_accuracy.tiff",width = 15, height = 6, units = "cm", dpi = 600)
