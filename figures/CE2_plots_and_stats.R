# Stats and plot generation for computational experiment 2
# 
# Written by Matt Biggs, 2016

library(reshape2)
library(ggplot2)

#-----------------------------------------------------------------
# Read in data and calulate empirical 95% confidence intervals
#-----------------------------------------------------------------
gms = c()
jacc = c()
times = c()
for (i in c(2,5,10,15,20,25,30))
{
  file_name = paste0("..\\data\\CE2_globalVsequential_",i,".tsv")
  seqData = read.table(file_name, row.names=NULL, sep = "\t", header=F)
  
  # Calculate the difference in model size (parsimony of the solution) and solution time
  seqData$glob_minus_sequential = seqData[,4] - seqData[,3]
  seqData$time_glob_div_seq = seqData[,6] / seqData[,5]
  
  # Estimate a 95% confidence interval assuming a normal distribution
  mean_gms = mean(seqData$glob_minus_sequential)
  quant_5_95 = quantile(seqData$glob_minus_sequential,c(0.025,0.975))
  tmp_gms = c(i,mean_gms,quant_5_95[1],quant_5_95[2])
  gms = rbind(gms,tmp_gms)
  
  mean_jacc = mean(seqData[,1])
  quant_5_95 = quantile(seqData[,1],c(0.025,0.975))
  tmp_jacc = c(i,mean_jacc,quant_5_95[1],quant_5_95[2])
  jacc = rbind(jacc,tmp_jacc)
  
  mean_timeFold = mean(seqData$time_glob_div_seq)
  quant_5_95 = quantile(seqData$time_glob_div_seq,c(0.025,0.975))
  tmp_timeFold = c(i,mean_timeFold,quant_5_95[1],quant_5_95[2])
  times = rbind(times,tmp_timeFold)
}
colnames(gms) = c("Num.Growth.Conditions","Global.Minus.Seq","lower","upper")
colnames(jacc) = c("Num.Growth.Conditions","Ave.Jaccard.Sim","lower","upper")
colnames(times) = c("Num.Growth.Conditions","Sol.Time.Fold.Change","lower","upper")
  
gms_df = data.frame(gms)
jacc_df = data.frame(jacc)
time_df = data.frame(times)

#-----------------------------------------------------------------
# Plot difference in network sizes
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

#-----------------------------------------------------------------
# Plot solution time fold change (Global/Sequential)
#-----------------------------------------------------------------
p1 = ggplot(time_df,aes(x=Num.Growth.Conditions,y=Sol.Time.Fold.Change)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymax=upper,ymin=lower)) +
  theme_bw() +
  scale_y_log10() +
  xlab("Num.Growth.Conditions") + 
  ylab("Solve Time Fold Change (Global/Seq.)") +
  theme(legend.position="none", text = element_text(size=12))
print(p1)

ggsave("CE2_solution_time_fold_change.tiff",width = 15, height = 9, units = "cm", dpi = 600)