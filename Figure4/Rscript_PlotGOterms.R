#GO term
library(ggplot2)
data<-read.csv("Selected_F_control_Cluster_GO.txt",header = T,sep="\t")
data[is.na(data)]<-1
d<-gather(data,cluster,Pvalue,cluster1:cluster9)
name<-as.character(data$GO);name<-name[length(name):1]
ggplot(d, aes(x = cluster,y = GO,col=-log10(Pvalue)))+
  geom_point(size=3) + 
  theme_classic()+ylab("GO")+
  scale_color_gradient2(low="blue",mid="white", high="red" )+
  scale_y_discrete(limits=name)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8,vjust = 0.5), 
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),legend.title =element_text(size=7) )