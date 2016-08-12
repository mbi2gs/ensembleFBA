# Stats and plot generation for computational experiment 3
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Read in data and calulate empirical 95% confidence intervals
#-----------------------------------------------------------------
file_name = paste0("..\\data\\CE3_globalVsequential.tsv")
data = read.table(file_name, row.names=NULL, sep = "\t", header=F)

colnames(data) = c("Seq.Sim.2.iPAU","Glob.Sim.2.iPAU","Seq.Size","Glob.Size","Seq.Solve.Time","Glob.Solve.Time")

data2 = select(data,Seq.Sim.2.iPAU,Glob.Sim.2.iPAU)
data3 = melt(data2)

wilcox.test(data2[,1],data[,2])

#-----------------------------------------------------------------
# Plot difference in network sizes
#-----------------------------------------------------------------
p1 = ggplot(data3,aes(x=factor(variable),y=value)) + 
     geom_boxplot() +
     theme_bw() +
     ylab("Jaccard Similarity to iPAU1129") + 
     theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=12))
print(p1)

ggsave("CE3_glob_minus_sequential.tiff",width = 6, height = 9, units = "cm", dpi = 600)
