# Stats and plot generation for computational experiment 1
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)

#-----------------------------------------------------------------
# Read in data and calulate empirical 95% confidence intervals
#-----------------------------------------------------------------
seq_urxns = c()
seq_jacc = c()
for (i in c(2,5,10,15,20,25,30))
{
  file_name = paste0("..\\data\\CE1_gapFillSequence_",i,"_gcs.tsv")
  seqData = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  
  # Estimate a 95% confidence interval assuming a normal distribution
  mean_urxns = mean(seqData[,2])
  quant_5_95 = quantile(seqData[,2],c(0.025,0.975))
  tmp_urxns = c(i,mean_urxns,quant_5_95[1],quant_5_95[2])
  seq_urxns = rbind(seq_urxns,tmp_urxns)
  
  mean_jacc = mean(seqData[,1])
  quant_5_95 = quantile(seqData[,1],c(0.025,0.975))
  tmp_jacc = c(i,mean_jacc,quant_5_95[1],quant_5_95[2])
  seq_jacc = rbind(seq_jacc,tmp_jacc)
}
colnames(seq_urxns) = c("Num.Growth.Conditions","Ave.Num.Unique.Rxns","lower","upper")
colnames(seq_jacc) = c("Num.Growth.Conditions","Ave.Jaccard.Sim.","lower","upper")

urxns_df = data.frame(seq_urxns)
jacc_df = data.frame(seq_jacc)

#-----------------------------------------------------------------
# Plot unique reactions per network
#-----------------------------------------------------------------
p1 = ggplot(urxns_df,aes(x=Num.Growth.Conditions,y=Ave.Num.Unique.Rxns)) + 
     geom_point(size=3) +
     geom_errorbar(aes(ymax=upper,ymin=lower)) +
     theme_bw() +
     xlab("Num.Growth.Conditions") + 
     ylab("Unique Rxns per Network") +
     theme(legend.position="none", text = element_text(size=12))
print(p1)

ggsave("CE1_average_unique_rxns.tiff",width = 15, height = 9, units = "cm", dpi = 600)

#-----------------------------------------------------------------
# Plot Jaccard similarities between networks
#-----------------------------------------------------------------
p1 = ggplot(jacc_df,aes(x=Num.Growth.Conditions,y=Ave.Jaccard.Sim.)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymax=upper,ymin=lower)) +
  theme_bw() +
  xlab("Num.Growth.Conditions") + 
  ylab("Jaccard Similarity") +
  theme(legend.position="none", text = element_text(size=12))
print(p1)

ggsave("CE1_jaccard_similarity.tiff",width = 15, height = 9, units = "cm", dpi = 600)