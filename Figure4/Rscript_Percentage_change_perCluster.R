# Pie plot show the percentage of changes in each cluster
library(ggplot2)
library(tidyr);

data<-read.table("Percentage_change_perCluster.txt",header = T)

#female
dat<-data[1:9,];dat$cluster<-rownames(dat);dat<-dat[,c(7,1:6)]; da<-gather(dat,Exp,Percent,BPA10mg:TBT);da$Type<-"Yes"
da2<-da;da2$Percent<-1-da2$Percent;da2$Type<-"No";d<-rbind(da,da2)
ggplot(d, aes(x=1, y = Percent, fill = Type ))+geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=c("white","orange"))+facet_grid(cluster ~ Exp)+coord_polar("y")+ theme_classic()+ 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.text.y  = element_blank(),axis.title.y = element_blank(),strip.text.x = element_text(size=8))

#male
dat<-data[10:18,];dat$cluster<-rownames(dat);dat<-dat[,c(7,1:6)]; da<-gather(dat,Exp,Percent,BPA10mg:TBT);da$Type<-"Yes"
da2<-da;da2$Percent<-1-da2$Percent;da2$Type<-"No";d<-rbind(da,da2)
