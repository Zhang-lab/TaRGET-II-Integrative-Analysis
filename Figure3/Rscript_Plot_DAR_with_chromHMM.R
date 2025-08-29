# DAR with chromHMM state from control to exposure
#BPA10mg_M as example
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

data<-read.table("New_Ctrl_BPA10mg_M_chrHMM_Li_adt_M_BA_BPA10mg.txt")
data$V19<-paste(data$V1,data$V2,data$V3,sep=",")
da1<-data[,c(19,17)];colnames(da1)<-c("Region","CRE");da1$Type<-"Ctrl_M";
da2<-data[,c(19,18)];colnames(da2)<-c("Region","CRE");da2$Type<-"BPA10mg_M";
result<-rbind(da1,da2);
result$CRE<-factor(result$CRE,
                   levels = c("AP","Promoter","Enhancer1","Enhancer2","OCR","None","LCR"))
ggplot(result, aes(x=Type, stratum=CRE, alluvium =Region, 
                   fill=CRE, label = CRE)) +
  scale_fill_brewer(palette = "Set3") +
  geom_flow(stat = "alluvium",width=0.1)+
  geom_stratum(width = 0.1,alpha = .75) +
  scale_x_discrete(limits = c("Ctrl_M", "BPA10mg_M"),expand = c(.05, .05))+
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 8),axis.title=element_text(size=8),
        plot.title = element_text(size=8),strip.text.x = element_text(size=8))
