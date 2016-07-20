# Stats and plot generation for computational experiment 4
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE5_networkAccuracies.tsv"
networkAccuracies = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(networkAccuracies) = c("Accuracy","N_gcs")

file_name = "..\\data\\CE5_ensembleAccuracy.tsv"
ensembleAccuracy = read.table(file_name, row.names=NULL, sep = "\t", header=F)
ensembleAccuracy = cbind(ensembleAccuracy,c(1,1,1,1,1,1,1))
colnames(ensembleAccuracy) = c("nMean","N_gcs","group")

# Summarize network data
networkSummaries = c()
for( i in c(2,5,10,15,20,25,30))
{
  data_i =  filter(networkAccuracies, N_gcs == i)
  tmpRange = c(mean(data_i$Accuracy),max(data_i$Accuracy),min(data_i$Accuracy),i)
  networkSummaries = rbind(networkSummaries,tmpRange)
}
networkSummaries = data.frame(networkSummaries)
colnames(networkSummaries) = c("nMean","nMax","nMin","N_gcs")

#-----------------------------------------------------------------
# Plot difference in network sizes
#-----------------------------------------------------------------
p1 = ggplot(networkSummaries,aes(x=N_gcs,y=nMean,group=N_gcs)) + 
     geom_pointrange(size=1.2,aes(ymin=nMax, ymax=nMin)) +
     theme_bw() +
     geom_line(data=ensembleAccuracy,size=1.5,color="blue",aes(group=group)) +
     xlab("Num.Growth.Conditions") + 
     ylab("Accuracy") +
     theme(legend.position="none", text = element_text(size=12))
print(p1)

ggsave("CE5_ensemble_accuracy.tiff",width = 15, height = 9, units = "cm", dpi = 600)
