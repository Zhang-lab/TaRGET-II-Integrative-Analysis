library(RUVSeq)
library(edgeR)
library(readr)
library(ggplot2)
library(DESeq2)
library("RColorBrewer")

#Rscript DAR.R "AL" "TCDD" "Li_M_adt" "Ctrl Li_M_adt"
#Liver adult TCDD male as example

args=commandArgs(T)

peak<-read.table("LAB_liver_correct.txt",header=T)
file<-paste0(args[1],"_Li_adt_peak.txt")
listp<-read.table(file)
lp<-as.character(listp$V1)
pk<-rownames(peak)
ll<-lp[lp %in% pk]
peak<-peak[ll,]

cond11<-grep(args[2],colnames(peak))
cond11_name<-colnames(peak)[cond11]
cond1<-grep(args[3],cond11_name)
cond1_name<-cond11_name[cond1]
cond1_name<-cond1_name[grep(args[1],cond1_name)]

cond22<-grep(args[4],colnames(peak))
cond22_name<-colnames(peak)[cond22]
cond2<-grep(args[5],cond22_name)
cond2_name<-cond22_name[cond2]

cn<-c(cond1_name,cond2_name)
countdata<-peak[,cn]
countdata<-round(countdata*30)

n<-ncol(countdata)
l<-(length(cond1_name))
con1<-paste0(args[2],"_",args[3])
con2<-paste0(args[4],"_",args[5])

keep=rowSums(cpm(countdata)>2)>=l
countdata=countdata[keep,]

x <- factor(c(rep(con1,length(cond1_name)), rep(con2,length(cond2_name))),levels = c(con2, con1))
set <- newSeqExpressionSet(as.matrix(countdata),phenoData=data.frame(x,row.names=colnames(countdata)))

wtest=DGEList(counts=countdata,group=x)
wtest=estimateCommonDisp(wtest, verbose=TRUE)
wtest=estimateTagwiseDisp(wtest)
wcpm=cpm(wtest$pseudo.counts )
wet=exactTest(wtest)
wntop=topTags(wet,n=dim(wtest)[1])
wstop=wntop$table[order(rownames(wntop$table)),]
wsexp=wcpm[order(rownames(wcpm)),]
wnftab=cbind(wstop,wsexp)
wnnftab=wnftab[order(wnftab[,4]),]
wftab=wnnftab
out=paste0(args[1],"_DAR_",con1,"_",con2,"_DAR.txt")
write.table(wftab,file=out,sep="\t",quote=FALSE)
wftab=wftab[abs(wftab[,1])>log2(1.5) & wftab[,4]<0.01,]
out=paste0(args[1],"_DAR_",con1,"_",con2,"_DAR_filter.txt")
write.table(wftab,file=out,sep="\t",quote=FALSE)

