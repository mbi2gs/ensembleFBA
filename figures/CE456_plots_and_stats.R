# Stats and plot generation for computational experiments 4, 5 and 6
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Combine data for 20 growth conditions across CEs 4, 5 and 6
#-----------------------------------------------------------------
netPrecision_20gcs = data.frame()
ensPrecision_20gcs = data.frame()
netRecall_20gcs = data.frame()
ensRecall_20gcs = data.frame()
netAccuracies_20gcs = data.frame()
ensAccuracies_20gcs = data.frame()

#-----------------------------------------------------------------
# Read in CE 4, 5 and 6 data
#-----------------------------------------------------------------
for(i in c("CE4","CE5","CE6"))
{
  file_name = paste0("..\\data\\",i,"_networkPrecision.tsv")
  networkPrecision = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  networkPrecision = cbind(networkPrecision,CE = rep(i,nrow(networkPrecision)))
  colnames(networkPrecision) = c("Precision","N_gcs","CE")
  netPrecision_20gcs = rbind(netPrecision_20gcs,filter(networkPrecision,N_gcs == 20))
  
  file_name = paste0("..\\data\\",i,"_ensemblePrecision.tsv")
  ensemblePrecision = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  ensemblePrecision = cbind(ensemblePrecision,CE = rep(i,nrow(ensemblePrecision))) 
  colnames(ensemblePrecision) = c("Precision_t1","Precision_tHalf","Precision_tAll","N_gcs","CE")
  ensemblePrecision_df = melt(ensemblePrecision,id=c("N_gcs","CE"))
  ensPrecision_20gcs = rbind(ensPrecision_20gcs,filter(ensemblePrecision_df,N_gcs == 20))
  
  file_name = paste0("..\\data\\",i,"_networkRecall.tsv")
  networkRecall = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  networkRecall = cbind(networkRecall,CE = rep(i,nrow(networkRecall))) 
  colnames(networkRecall) = c("Recall","N_gcs","CE")
  netRecall_20gcs = rbind(netRecall_20gcs,filter(networkRecall,N_gcs == 20))
  
  file_name = paste0("..\\data\\",i,"_ensembleRecall.tsv")
  ensembleRecall = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  ensembleRecall = cbind(ensembleRecall,CE = rep(i,nrow(ensembleRecall))) 
  colnames(ensembleRecall) = c("Recall_t1","Recall_tHalf","Recall_tAll","N_gcs","CE")
  ensembleRecall_df = melt(ensembleRecall,id=c("N_gcs","CE"))
  ensRecall_20gcs = rbind(ensRecall_20gcs,filter(ensembleRecall_df,N_gcs == 20))
  
  file_name = paste0("..\\data\\",i,"_networkAccuracy.tsv")
  networkAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  networkAccuracies = cbind(networkAccuracies,CE = rep(i,nrow(networkAccuracies))) 
  colnames(networkAccuracies) = c("Accuracy","N_gcs","CE")
  netAccuracies_20gcs = rbind(netAccuracies_20gcs,filter(networkAccuracies,N_gcs == 20))
  
  file_name = paste0("..\\data\\",i,"_ensembleAccuracy.tsv")
  ensembleAccuracy = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  ensembleAccuracy = cbind(ensembleAccuracy,CE = rep(i,nrow(ensembleAccuracy))) 
  colnames(ensembleAccuracy) = c("Accuracy_t1","Accuracy_tHalf","Accuracy_tAll","N_gcs","CE")
  ensembleAccuracy_df = melt(ensembleAccuracy,id=c("N_gcs","CE"))
  ensAccuracies_20gcs = rbind(ensAccuracies_20gcs,filter(ensembleAccuracy_df,N_gcs == 20))
}

#-----------------------------------------------------------------
# Summarize individual network data
#-----------------------------------------------------------------
networkPrecisionSummaries = c()
for( i in c("CE4","CE5","CE6"))
{
  data_i =  filter(netPrecision_20gcs,  CE == i)
  tmpRange = c(mean(data_i$Precision),max(data_i$Precision),min(data_i$Precision))
  networkPrecisionSummaries = rbind(networkPrecisionSummaries,tmpRange)
}
networkPrecisionSummaries = data.frame(networkPrecisionSummaries)
networkPrecisionSummaries = cbind(networkPrecisionSummaries,CE = c("CE4","CE5","CE6"))
colnames(networkPrecisionSummaries) = c("nMean","nMax","nMin","CE")

networkRecallSummaries = c()
for( i in c("CE4","CE5","CE6"))
{
  data_i =  filter(netRecall_20gcs,  CE == i)
  tmpRange = c(mean(data_i$Recall),max(data_i$Recall),min(data_i$Recall))
  networkRecallSummaries = rbind(networkRecallSummaries,tmpRange)
}
networkRecallSummaries = data.frame(networkRecallSummaries)
networkRecallSummaries = cbind(networkRecallSummaries,CE = c("CE4","CE5","CE6"))
colnames(networkRecallSummaries) = c("nMean","nMax","nMin","CE")

networkAccuracySummaries = c()
for( i in c("CE4","CE5","CE6"))
{
  data_i =  filter(netAccuracies_20gcs,  CE == i)
  tmpRange = c(mean(data_i$Accuracy),max(data_i$Accuracy),min(data_i$Accuracy))
  networkAccuracySummaries = rbind(networkAccuracySummaries,tmpRange)
}
networkAccuracySummaries = data.frame(networkAccuracySummaries)
networkAccuracySummaries = cbind(networkAccuracySummaries,CE = c("CE4","CE5","CE6"))
colnames(networkAccuracySummaries) = c("nMean","nMax","nMin","CE")

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Precision)
#-----------------------------------------------------------------------
p1 = ggplot(networkPrecisionSummaries,aes(x=CE,y=nMean,group=CE)) + 
     geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
     theme_bw() +
     geom_point(data=ensPrecision_20gcs,size=2.5,aes(x=CE,y=value,color=variable,shape=variable)) +
     xlab("Ensemble Type") + 
     ylab("Precision") +
     theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE456_ensemble_precision.tiff",width = 5, height = 9, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Recall)
#-----------------------------------------------------------------------
p1 = ggplot(networkRecallSummaries,aes(x=CE,y=nMean,group=CE)) + 
  geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
  theme_bw() +
  geom_point(data=ensRecall_20gcs,size=2.5,aes(x=CE,y=value,color=variable,shape=variable)) +
  xlab("Ensemble Type") + 
  ylab("Recall") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE456_ensemble_recall.tiff",width = 5, height = 9, units = "cm", dpi = 600)

#-----------------------------------------------------------------------
# Plot comparision between individual networks and ensemble (Precision)
#-----------------------------------------------------------------------
p1 = ggplot(networkAccuracySummaries,aes(x=CE,y=nMean,group=CE)) + 
  geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
  theme_bw() +
  geom_point(data=ensAccuracies_20gcs,size=2.5,aes(x=CE,y=value,color=variable,shape=variable)) +
  xlab("Ensemble Type") + 
  ylab("Accuracy") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE456_ensemble_accuracy.tiff",width = 5, height = 9, units = "cm", dpi = 600)
