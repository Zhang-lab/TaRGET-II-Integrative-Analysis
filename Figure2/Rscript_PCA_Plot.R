library(plotly)
library(ggplot2)
library(ggpubr)
library(RUVSeq)
library(RColorBrewer)

#1. transcriptome data from all tissues
data<-read.table("CountTable_across_allTissues.txt",header=T)
d<-cpm(data+1)
n<-colnames(d)

lab<-c()
trt<-c()
sex<-c()
tis<-c()
stage<-c()
for(j in 1:length(n)){
  lab[j]<-strsplit(as.character(n[j]),"_")[[1]][1]
  trt[j]<-strsplit(as.character(n[j]),"_")[[1]][3]
  if(trt[j] != "Ctrl"){
    if((lab[j] == "BI") || (lab[j] == "MU")){
      trt[j]<-paste(lab[j],trt[j],sep=".")
    }
  }
  sex[j]<-strsplit(as.character(n[j]),"_")[[1]][6]
  tis[j]<-strsplit(as.character(n[j]),"_")[[1]][5]
  stage[j]<-strsplit(as.character(n[j]),"_")[[1]][7]
}
lab<-paste0(lab,"_RNA")
tissue<-paste(tis,stage,sep="_")
l<-sum(!duplicated(trt))

m<-t(log2(d))
pca<-prcomp(m,scale. = T)
pca_plot<-as.data.frame(pca$x[,1:3])
pca_plot$lab<-lab
pca_plot$exp<-trt
pca_plot$sex<-sex
pca_plot$tissue<-tissue
pca_plot$name<-rownames(m)
sdev_df <- data.frame(components = paste0("PC", 1:length(pca$sdev)), value = pca$sdev ** 2)
sdev_df$value <- sdev_df$value / sum(sdev_df$value) * 100
pca_plot$pc1_perc <- sdev_df$value[1] %>% round(., digits = 2)
pca_plot$pc2_perc <- sdev_df$value[2] %>% round(., digits = 2)
pca_plot$pc3_perc <- sdev_df$value[3] %>% round(., digits = 2)

cols <- brewer.pal(length(unique(pca_plot$exp)), "Set3");
p<-plot_ly(pca_plot,x=~PC1,y=~PC2,z=~PC3,color=~tissue,colors=cols[1:l]
) %>% add_markers( text = ~name)

#2.PCA for liver and blood data considering the differnet exposures and sex
# 5 month liver ATAC-seq as example
data<-read.table("LAB_liver_adt_correct_readCount.txt",header=T)
d<-cpm(data+1)
n<-colnames(d)
lab<-c()
trt<-c()
sex<-c()
tis<-c()
stage<-c()
for(j in 1:length(n)){
  lab[j]<-strsplit(as.character(n[j]),"_")[[1]][1]
  trt[j]<-strsplit(as.character(n[j]),"_")[[1]][3]
  if(trt[j] != "Ctrl"){
    if((lab[j] == "BI") || (lab[j] == "MU")){
      trt[j]<-paste(lab[j],trt[j],sep=".")
    }
  }
  sex[j]<-strsplit(as.character(n[j]),"_")[[1]][6]
  tis[j]<-strsplit(as.character(n[j]),"_")[[1]][5]
  stage[j]<-strsplit(as.character(n[j]),"_")[[1]][7]
}
lab<-paste0(lab,"_ATAC")
tissue<-paste(tis,stage,sep="_")
l<-sum(!duplicated(trt))

m<-t(log2(d))
pca<-prcomp(m,scale. = T)
pca_plot<-as.data.frame(pca$x[,1:3])
pca_plot$lab<-lab
pca_plot$exp<-trt
pca_plot$sex<-sex
pca_plot$tissue<-tissue
pca_plot$name<-rownames(m)
sdev_df <- data.frame(components = paste0("PC", 1:length(pca$sdev)), value = pca$sdev ** 2)
sdev_df$value <- sdev_df$value / sum(sdev_df$value) * 100
pca_plot$pc1_perc <- sdev_df$value[1] %>% round(., digits = 2)
pca_plot$pc2_perc <- sdev_df$value[2] %>% round(., digits = 2)
pca_plot$pc3_perc <- sdev_df$value[3] %>% round(., digits = 2)
cols <- brewer.pal(length(unique(pca_plot$exp)), "Set3");
p<-plot_ly(pca_plot,x=~PC1,y=~PC2,z=~PC3,color=~exp,colors=cols[1:l],
           symbol =~sex ,symbols = c("cross","circle-dot")
) %>% add_markers( text = ~name)
htmlwidgets::saveWidget(as_widget(p), "PCA3D_Liver_adt_correction.html")

#3. PCA plots for each exposure
list<-c("AL","BA","BI","DO","MU","WK","ZB")
colours<-c("#F8766D","#00BFC4","#619CFF")
tp<-"liver_wl_RUVrk4"
k<-4
for(f in list){
  file<-paste0(f,"_RUVr_k",k,".bed")
  d<-read.table(file,header=T)
  n<-colnames(d)
  lab<-c()
  trt<-c()
  sex<-c()
  for(j in 1:length(n)){
    trt[j]<-strsplit(as.character(n[j]),"_")[[1]][3]
    sex[j]<-strsplit(as.character(n[j]),"_")[[1]][6]
  }
  l<-sum(!duplicated(trt))
  d<-cpm(d+1)
  m<-t(log2(d))
  pca<-prcomp(m,scale. = T)
  pca_plot<-as.data.frame(pca$x[,1:3])
  pca_plot$exp<-trt
  pca_plot$sex<-sex
  pca_plot$name<-rownames(m)
  sdev_df <- data.frame(components = paste0("PC", 1:length(pca$sdev)), value = pca$sdev ** 2)
  sdev_df$value <- sdev_df$value / sum(sdev_df$value) * 100
  pca_plot$pc1_perc <- sdev_df$value[1] %>% round(., digits = 2)
  pca_plot$pc2_perc <- sdev_df$value[2] %>% round(., digits = 2)
  pca_plot$pc3_perc <- sdev_df$value[3] %>% round(., digits = 2)
  
  out<-paste0(f,"_",tp,"_by_exp_sex.html")
  p<-plot_ly(pca_plot,x=~PC1,y=~PC2,z=~PC3,color=~exp,colors=colours[1:l],
             symbol =~sex ,symbols = c("cross","circle-dot")
  ) %>% add_markers( text = ~name)
  htmlwidgets::saveWidget(as_widget(p), out)
  
  out<-paste0(f,"_",tp,"_by_exp_sex.pdf")
  pdf(file=out,5,8)
  title<-paste0(f,"_",tp)
  p1 <- ggplot(pca_plot, aes(x = PC1, y = PC2, col = exp, shape = sex)) +
    geom_point(size = 3) +
    xlab(paste0("PC1 (", unique(pca_plot$pc1_perc), "%)")) +
    ylab(paste0("PC2 (", unique(pca_plot$pc2_perc), "%)")) +
    scale_shape_manual(values = c(4, 19)) +
    scale_color_manual(values=colours[1:l])+
    ggtitle(title)+theme_classic()
  theme(title=element_text(size=14,face="bold"))+
    theme(legend.text=element_text(size=14))
  p2 <- ggplot(pca_plot, aes(x = PC1, y = PC3, col = exp, shape = sex)) +
    geom_point(size = 3) +
    xlab(paste0("PC1 (", unique(pca_plot$pc1_perc), "%)")) +
    ylab(paste0("PC3 (", unique(pca_plot$pc3_perc), "%)")) +
    scale_shape_manual(values = c(4, 19)) +
    scale_color_manual(values=colours[1:l])+
    ggtitle(title)+theme_classic()+
    theme(title=element_text(size=14,face="bold"))+
    theme(legend.text=element_text(size=14))
  print(ggarrange(p1, p2, nrow = 2))
  dev.off()
  
}



