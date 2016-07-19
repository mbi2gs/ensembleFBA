# Stats and plot generation for computational experiment 2
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)

#-----------------------------------------------------------------
# Read in data and calulate empirical 95% confidence intervals
#-----------------------------------------------------------------
seq_gms = c()
seq_jacc = c()
for (i in c(2,5,10,15,20,25,30))
{
  file_name = paste0("..\\data\\CE2_globalVsequential_",i,".tsv")
  seqData = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  
  # Calculate the difference in model size (parsimony of the solution)
  seqData$glob_minus_sequential = seqData[,4] - seqData[,3]
  
  # Estimate a 95% confidence interval assuming a normal distribution
  mean_gms = mean(seqData$glob_minus_sequential)
  quant_5_95 = quantile(seqData$glob_minus_sequential,c(0.025,0.975))
  tmp_gms = c(i,mean_gms,quant_5_95[1],quant_5_95[2])
  seq_gms = rbind(seq_gms,tmp_gms)
  
  mean_jacc = mean(seqData[,1])
  quant_5_95 = quantile(seqData[,1],c(0.025,0.975))
  tmp_jacc = c(i,mean_jacc,quant_5_95[1],quant_5_95[2])
  seq_jacc = rbind(seq_jacc,tmp_jacc)
}
colnames(seq_gms) = c("Num.Growth.Conditions","Global.Minus.Seq","lower","upper")
colnames(seq_jacc) = c("Num.Growth.Conditions","Ave.Jaccard.Sim","lower","upper")

gms_df = data.frame(seq_gms)
jacc_df = data.frame(seq_jacc)

#-----------------------------------------------------------------
# Plot unique reactions per network
#-----------------------------------------------------------------
p1 = ggplot(gms_df,aes(x=Num.Growth.Conditions,y=Global.Minus.Seq)) + 
     geom_point(size=3) +
     geom_errorbar(aes(ymax=upper,ymin=lower)) +
     theme_bw() +
     xlab("Num.Growth.Conditions") + 
     ylab("Num.Rxns.(Glob.-Seq.)") +
     theme(legend.position="none", text = element_text(size=12))
print(p1)

ggsave("CE2_glob_minus_sequential.tiff",width = 15, height = 9, units = "cm", dpi = 600)

#-----------------------------------------------------------------
# Plot Jaccard similarities between networks
#-----------------------------------------------------------------
p1 = ggplot(jacc_df,aes(x=Num.Growth.Conditions,y=Ave.Jaccard.Sim)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymax=upper,ymin=lower)) +
  theme_bw() +
  xlab("Num.Growth.Conditions") + 
  ylab("Jaccard Similarity") +
  theme(legend.position="none", text = element_text(size=12))
print(p1)

ggsave("CE2_jaccard_similarity.tiff",width = 15, height = 9, units = "cm", dpi = 600)