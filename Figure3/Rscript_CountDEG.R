
# count the the DEG number and across 3 stages
data<-read.table("20240130_liver_blood_treatment_deg_after_PseudoCtrl_using_all_control.txt",header = T)
dat<-data[data$sig=="TRUE",];dat<-dat[,c(7,9,13:17)]
z<-dat[,c(3,5:7)];z<-z[!(duplicated(z)),];result<-c()
for(i in 1:nrow(z)){
  ts<-as.character(z[i,2]);age<-as.character(z[i,3]);sx<-as.character(z[i,4]);dir<-as.character(z[i,1]);
  da<-dat[dat$tissue==ts & dat$stage==age & dat$sex==sx & dat$direction ==dir,];d<-da[,c(2,3,4)];d<-d[!(duplicated(d)),];
  a<-as.data.frame(table(d$gene));b<-a[a$Freq>1,]; d$Share<-d$gene %in% b$Var1;d$Type<-ifelse(d$Share=="TRUE","Multiple",d$group)
  d<-d[,c(1,5)];d<-d[!(duplicated(d)),];x<-as.data.frame(table(d$Type))
  x$Tissue<-ts; x$Stage<-age; x$Sex<-sx; x$Direction<-dir
  result<-rbind(result,x)
}

#different stages
dat<-data[data$tissue=="liver" & data$sig=="TRUE",]

#3weeks, 20weeks, 40weeks
#Female, Male
a<-c();s<-c()
for(st in c("3weeks","20weeks","40weeks")){
  for(sx in c("Female","Male")){
    da<-dat[dat$stage==st & dat$sex==sx,];d<-da[,c(2,9,14)];
    f<-spread(d,group,log2FoldChange);f[is.na(f)]<-0
    gene<-inf[f$gene,];f<-cbind(gene[,c(1,4)],f)
    k<-f[,c(3,2,4:ncol(f))]
    out<-paste("liver",st,sx,"DEG.txt",sep="_")
    write.table(k,out,sep="\t",quote=F,row.names = F)
    
    b<-paste(st,sx,sep="_")
    a<-c(a,b)
    s<-c(s,dim(k)[1])
  }
}
