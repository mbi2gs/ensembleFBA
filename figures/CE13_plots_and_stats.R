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

ggsave("CE13_distribution_small_molecules.tiff",width = 5, height = 7, units = "cm", dpi = 600)
