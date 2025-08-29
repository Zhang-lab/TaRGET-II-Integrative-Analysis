#1. fold change of blood-liver common DEGs
library(pheatmap)
library(UpSetR)

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(3/paletteLength, 2, length.out=floor(paletteLength/2)))
#female
data<-read.table("Li_Bl_adt_F_share_gene_logFC.txt",header = T)
d<-data[,1:18]; ann<-as.data.frame(data[,19]);rownames(ann)<-rownames(d);colnames(ann)<-"Type"
pheatmap(d,show_rownames = F,cluster_rows = F,cluster_cols = F, fontsize = 8,
         breaks = myBreaks,color=myColor,annotation_row = ann)
common<-read.table("common_liver_blood_DEG_Female_table.txt",header = T)
upset(common,nsets = 9,nintersects = 1000,point.size = 1.25,line.size = 0.3,
      sets = colnames(common),keep.order=F,mb.ratio = c(0.6, 0.4))
#male
data<-read.table("Li_Bl_adt_M_share_gene_logFC.txt",header = T)
d<-data[,1:18]; ann<-as.data.frame(data[,19]);rownames(ann)<-rownames(d);colnames(ann)<-"Type"
pheatmap(d,show_rownames = F,cluster_rows = F,cluster_cols = F, fontsize = 8,
         breaks = myBreaks,color=myColor,annotation_row = ann)
common<-read.table("common_liver_blood_DEG_Male_table.txt",header = T)
upset(common,nsets = 9,nintersects = 1000,point.size = 1.25,line.size = 0.3,
      sets = colnames(common),keep.order=F,mb.ratio = c(0.6, 0.4))

#2. Gene expression of liver-blood common DEGs
data<-read.table("Gene_CPM_table_alltissue.txt",header = T)
list<-c("Ppp1r3b","Gbp8","Jdp2")
sample<-colnames(data);s<-sample[c(grep("Li_F_adt",sample),grep("Bl_F_adt",sample))];
s1<-s[grep("Ctrl",s)];s2<-s[grep("TCDD",s)];s<-c(s1,s2)
dat<-data[,s];da<-dat[list,];da$Gene<-rownames(da);f<-gather(da,sample,CPM,1:(ncol(da)-1))
n<-as.character(f$sample);lab<-c();trt<-c();sex<-c();tissue<-c();for(j in 1:length(n)){
  lab[j]<-strsplit(as.character(n[j]),"_")[[1]][1];trt[j]<-strsplit(as.character(n[j]),"_")[[1]][3]
  tissue[j]<-strsplit(as.character(n[j]),"_")[[1]][5];
  sex[j]<-strsplit(as.character(n[j]),"_")[[1]][6]};f$type<-trt;f$sex<-sex;f$tissue<-tissue

f$type<-factor(f$type,levels = c("Ctrl","TCDD"))
ggplot(f, aes(type, CPM,fill=type))+scale_fill_brewer(palette = "Set3") +
  geom_boxplot(outlier.size = 2)+theme_classic()+facet_wrap(Gene~tissue,scales = "free",ncol = 2,nrow = 3)+
  theme(axis.text.x = element_text(size=10,angle = 90, hjust = 1,vjust = 0.5),axis.text.y = element_text(size=10))

g1<-c();g2<-c();g3<-c();g4<-c()
for(l in list){for(i in c("Bl","Li")){
  c1<-as.numeric(f[f$Gene==l & f$tissue==i & f$type=="Ctrl",]$CPM);
  c2<-as.numeric(f[f$Gene==l & f$tissue==i & f$type=="TCDD",]$CPM)
  a<-t.test(c1,c2)$p.value
  g1<-c(g1,l);g2<-c(g2,i);g4<-c(g4,a)}}
g<-as.data.frame(cbind(g1,g2,g4));colnames(g)<-c("Gene","Tissue","TtestPvalue")
write.table(g,"T-test-Pvalue.txt",sep="\t",quote=F,row.names = F)





