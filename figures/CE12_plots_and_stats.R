# Stats and plot generation for computational experiment 7
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------
file_name = "..\\data\\CE12_stochastic_rxn_dist.tsv"
rxnDistributionS = read.table(file_name, row.names=NULL, sep = "\t", header=F)
colnames(rxnDistributionS) = c("Correct","Incorrect","InitialDist")
rxnDistributionS$sum = rxnDistributionS$Correct + rxnDistributionS$Incorrect
rxnDistributionS$proportion = rxnDistributionS$Correct / rxnDistributionS$sum
rxnDistributionS$proportion[is.nan(rxnDistributionS$proportion)] = 0
rxnDistributionS$bins = 0:(nrow(rxnDistributionS)-1)

#-----------------------------------------------------------------------
# Plot accuracy of ensembles as a function of ensemble size
#-----------------------------------------------------------------------
totalRxns = sum(rxnDistributionS$sum)
rxnDist_df = melt(rxnDistributionS[2:101,c(1,2,6)],id=c("bins"))
rxnDist_df$value = rxnDist_df$value / totalRxns
initDist_df = melt(rxnDistributionS[2:101,c(3,6)],id=c("bins"))
initDist_df$value = initDist_df$value / totalRxns
p1 = ggplot(rxnDist_df,aes(x=bins,y=value,group=factor(variable))) +
  geom_line(size=1,alpha=0.7,aes(color=factor(variable))) +
  geom_line(data=initDist_df,size=0.35,color="black") +
  theme_bw() +
  xlab("Number of GENREs") +
  ylab("Density") +
  scale_color_brewer(palette="Set1") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE12_reaction_distribution.tiff",width = 8, height = 5, units = "cm", dpi = 600)

p1 = ggplot(rxnDistributionS[2:101,],aes(x=bins,y=proportion)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("Number of GENREs") +
  ylab("") +
  theme(text = element_text(size=12),legend.position="none")
print(p1)

ggsave("CE12_proportion_correct.tiff",width = 8, height = 3, units = "cm", dpi = 600)
