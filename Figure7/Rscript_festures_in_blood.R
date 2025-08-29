
# 1. numner of DEGs, DARs and DMRs in different exposures of blood
library(ggplot2)
data<-read.table("Number_DEG_DAR_DMR_in_Blood.txt",header = T)

#number of up-regulated features in female and male blood
da<-data[(grep("up-",as.character(data$Direction))),]
ggplot(da, aes(x=Type, y=Count,fill=Direction)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_x_discrete(limits=c("DEG","DAR","DMR"))+
  theme_classic()+facet_grid(Sex~Exp,scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8,vjust = 0.5), 
        strip.text.x = element_text(size=6),
        legend.title = element_text(size = 8))

#number of down-regulated features in female and male blood
da<-data[(grep("down-",as.character(data$Direction))),]
ggplot(da, aes(x=Type, y=Count,fill=Direction)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_x_discrete(limits=c("DEG","DAR","DMR"))+
  theme_classic()+facet_grid(Sex~Exp,scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8,vjust = 0.5), 
        strip.text.x = element_text(size=6),
        legend.title = element_text(size = 8))



#2. GO term enrichment

library(clusterProfiler)
library(enrichplot)
library("ggnewscale")
library(DOSE)
library(org.Mm.eg.db)
mm<-org.Mm.eg.db

#female
data<-read.table("Adult_female_Exposure_specific_DEG.txt",header = T)
k<-data[rowSums(data!=0)>0,];group<-c()
for(i in 1:ncol(k)){
  kk<-k[k[,i]!=0,]
  my.symbols<-as.character(rownames(kk))
  IDs<-AnnotationDbi::select(mm, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
  Entrez<-unique(IDs[!is.na(IDs$ENTREZID),c(1,2)])
  group$x<-Entrez$ENTREZID
  names(group)<-colnames(k)[1:i]
}
xx <- compareCluster(group, fun="enrichGO",OrgDb=mm, pvalueCutoff=0.05,ont="BP")
xx <- pairwise_termsim(xx)
emapplot(xx,legend_n = 3,cex_line = 0.4,cex_label_category = 0.4)

result<-xx@compareClusterResult
for(x in 1:nrow(result)){
  symbols<- AnnotationDbi::select(mm, keys = strsplit(result[x,"geneID"], "/")[[1]], 
                                  columns = c("SYMBOL"), keytype = "ENTREZID")
  result[x,"geneSymbols"]<- paste(symbols$SYMBOL, collapse = ",")
}
write.table(result,"Adult_female_Exposure_specific_DEG_GO.txt",sep="\t",quote=F,row.names = F)

#male
data<-read.table("Adult_male_Exposure_specific_DEG.txt",header = T)
k<-data[rowSums(data!=0)>0,];group<-c()
for(i in 1:ncol(k)){
  kk<-k[k[,i]!=0,]
  my.symbols<-as.character(rownames(kk))
  IDs<-AnnotationDbi::select(mm, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
  Entrez<-unique(IDs[!is.na(IDs$ENTREZID),c(1,2)])
  group$x<-Entrez$ENTREZID
  names(group)<-colnames(k)[1:i]
}
xx <- compareCluster(group, fun="enrichGO",OrgDb=mm, pvalueCutoff=0.05,ont="BP")
xx <- pairwise_termsim(xx)
emapplot(xx,legend_n = 3,cex_line = 0.4,cex_label_category = 0.4)

result<-xx@compareClusterResult
for(x in 1:nrow(result)){
  symbols<- AnnotationDbi::select(mm, keys = strsplit(result[x,"geneID"], "/")[[1]], 
                                  columns = c("SYMBOL"), keytype = "ENTREZID")
  result[x,"geneSymbols"]<- paste(symbols$SYMBOL, collapse = ",")
}
write.table(result,"Adult_male_Exposure_specific_DEG_GO.txt",sep="\t",quote=F,row.names = F)



#3. distribution of DEGs shared by multiple exposures
library(pheatmap)

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(4/paletteLength, 2, length.out=floor(paletteLength/2)))

#female
data<-read.table("3.1 adult_female_common.txt",header = T)
#distribution
d<-data[,1:(ncol(data)-1)]
pheatmap(d,show_rownames = F,cluster_rows = T, ,fontsize = 8,breaks = myBreaks,color=myColor)
#shared by at least four exposures
d<-data[data$num>3,c(1:(ncol(data)-1))]
pheatmap(d,show_rownames = T,cluster_rows = T, ,fontsize = 8,breaks = myBreaks,color=myColor)
#shared by 2 or 3 exposure
library(UpSetR)
d<-data[data$num<4,c(1:(ncol(data)-1))]
d1<-d;d1[d1>0]<-1;d1[d1<0]<-0;colnames(d1)<-paste("Up",colnames(d1),sep="_")
d2<-d;d2[d2>0]<-0;d2[d2<0]<-1;colnames(d2)<-paste("Down",colnames(d2),sep="_")
k<-cbind(d1,d2)
upset(k,nsets = 18,nintersects = 1000,point.size = 1.25,
      line.size = 0.3,sets = colnames(k),keep.order=F,mb.ratio = c(0.6, 0.4))

#male
data<-read.table("3.2 adult_male_common.txt",header = T)
#distribution
d<-data[,1:(ncol(data)-1)]
pheatmap(d,show_rownames = F,cluster_rows = T, ,fontsize = 8,breaks = myBreaks,color=myColor)
#shared by at least four exposures
d<-data[data$num>3,c(1:(ncol(data)-1))]
pheatmap(d,show_rownames = T,cluster_rows = T, ,fontsize = 8,breaks = myBreaks,color=myColor)
#shared by 2 or 3 exposure
library(UpSetR)
d<-data[data$num<4,c(1:(ncol(data)-1))]
d1<-d;d1[d1>0]<-1;d1[d1<0]<-0;colnames(d1)<-paste("Up",colnames(d1),sep="_")
d2<-d;d2[d2>0]<-0;d2[d2<0]<-1;colnames(d2)<-paste("Down",colnames(d2),sep="_")
k<-cbind(d1,d2)
upset(k,nsets = 18,nintersects = 1000,point.size = 1.25,
      line.size = 0.3,sets = colnames(k),keep.order=F,mb.ratio = c(0.6, 0.4))









