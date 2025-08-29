
#changes of temporal expression patterns 
#BPA10mg_M as example

library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

data<-read.table("Cluster_BPA10mg_M_.txt")
x<-"Ctrl_Male"; y<-"BPA10mg_Male"
dat<-data[,c(11,13)];dat$Gene<-rownames(dat)
da1<-dat[,c(3,1)];colnames(da1)<-c("Gene","Cluster");da1$Type<-x
da2<-dat[,c(3,2)];colnames(da2)<-c("Gene","Cluster");da2$Type<-y
result<-rbind(da1,da2);result$Cluster<-as.factor(result$Cluster)

dat<-data[,c(11,13)];dat$Gene<-rownames(dat);dat<-dat[dat$Cluster_control != dat$Cluster_Exposure_Final,]
da1<-dat[,c(3,1)];colnames(da1)<-c("Gene","Cluster");da1$Type<-x
da2<-dat[,c(3,2)];colnames(da2)<-c("Gene","Cluster");da2$Type<-y
result<-rbind(da1,da2);result$Cluster<-as.factor(result$Cluster)

ggplot(result, aes(x=Type, stratum=Cluster, alluvium =Gene, 
                   fill=Cluster, label = Cluster)) +
  scale_fill_brewer(palette = "Set3") +
  geom_flow(stat = "alluvium",width=0.1)+
  geom_stratum(width = 0.1,alpha = .75) +
  scale_x_discrete(limits = c(x, y),expand = c(.05, .05))+
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title = element_text(size = 8),axis.title=element_text(size=8),
        plot.title = element_text(size=8),strip.text.x = element_text(size=8))
