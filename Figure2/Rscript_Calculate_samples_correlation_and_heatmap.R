###correlation for each exposure per sex####
data<-read.table("LAB_liver_adt_correct_readCount.txt")
data<-cpm(data)
group<-as.data.frame(cbind(c("AL","BA","BA","BI","DO","DO","MU","WK","ZB"),
                           c("TCDD","BPA10mg","BPA10ug","PM2.5","DEHP","Pb","PM2.5","TBT","As")))
#female|male
dat<-data[,grep("_F_",colnames(data))];
sex<-"Female"

#correlation
for(i in 1:nrow(group)){
  lab<-group[i,1];exp<-group[i,2];
  da<-dat[,grep(lab,colnames(dat))];
  d<-da[,c(grep(exp,colnames(da)),grep("Ctrl",colnames(da)))]
  d<-d[rowSums(d>10)>=length(grep("Ctrl",colnames(d))),];
  
  sample<-colnames(d);cor<-c()
  for(i in 1:ncol(d)){a<-as.numeric(d[,i]);c<-c();
  for(j in 1:ncol(d)){b<-as.numeric(d[,j]);c[j]<-cor(a,b,method="spearman")};cor<-cbind(cor,c)}
  cor<-as.data.frame(cor);rownames(cor)<-sample;colnames(cor)<-sample
  
  list<-colnames(d);trt<-c();
  for(x in 1:length(list)){
    y<-strsplit(list[x],"_")[[1]];
    if(y[3] != "Ctrl"){if((y[1] == "BI") || (y[1] == "MU")){y[3]<-paste(y[1],y[3],sep=".")}}
    trt[x]<-y[3];
  }
  ann<-as.data.frame(trt);rownames(ann)<-list;
  
  pdf(paste(lab,exp,sex,"cor.pdf",sep="_"),6,6)
  pheatmap(cor,annotation_row=ann,annotation_col =ann,fontsize = 7,cellwidth = 8,cellheight = 8)
  dev.off()
  
}
