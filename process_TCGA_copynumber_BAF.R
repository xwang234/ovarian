#!/usr/bin/env Rscript
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
njob=100
filelist=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/filelist.txt",header=T,sep="\t",stringsAsFactors=F)

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
allsamples=rownames(data_copynumber_filtered)
idxtumor=sapply(allsamples,function(mysample){
  res=which(grepl(paste0(mysample,"-01"),filelist$sampleid)==T & grepl("byallele",filelist$filename)==T)
  #for multiple files
  if (mysample=="TCGA-23-1023") res=res[2]
  if (mysample=="TCGA-29-2414") res=res[1]
  if (mysample=="TCGA-13-1817") res=res[2]
  if (mysample=="TCGA-13-1819") res=res[1]
  return(res)
})

idxnormal=sapply(allsamples,function(mysample){
  res=which(grepl(paste0(mysample,"-1"),filelist$sampleid)==T & grepl("byallele",filelist$filename)==T)
  return(res)  
})

idxnormal1=sapply(1:length(allsamples),function(i){
  if (length(idxnormal[[i]])==0)
  {
    res=NA
  }else
  {
    res=idxnormal[[i]]
  }
  names(res)=names(idxnormal)[i]
  return(res)
})
sum(is.na(idxnormal1))
#[1] 12 samples without normal allele files
idxnormal=idxnormal1

snp_anno=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno=snp_anno[,1:4]
library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
mpi.bcast.Robj2slave(snp_anno)
mpi.bcast.Robj2slave(filelist)
mpi.bcast.Robj2slave(idxnormal)
mpi.bcast.Robj2slave(idxtumor)
mpi.bcast.Robj2slave(allsamples)


library(DNAcopy)
mpi.bcast.cmd(library(DNAcopy))
CBS_allele=function(i)
{
  mysample=allsamples[i]
  tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[i]],"/",filelist$filename[idxtumor[i]]),
                        skip = 1, sep="\t",header=T,stringsAsFactors = F)
  tumortable$Signal_A[which(tumortable$Signal_A<0)]=0
  tumortable$Signal_B[which(tumortable$Signal_B<0)]=0
  
  colnames(tumortable)[1]="Probe.Set.ID"
  
  #if it has normal allele files, find heterogeneity probes first
  if (!is.na(idxnormal[i]))
  {
    normaltable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxnormal[i]],"/",filelist$filename[idxnormal[i]]),
                           skip = 1, sep="\t",header=T,stringsAsFactors = F)
    normaltable$Signal_A[which(normaltable$Signal_A<0)]=0
    normaltable$Signal_B[which(normaltable$Signal_B<0)]=0
    af=normaltable$Signal_A/(normaltable$Signal_A+normaltable$Signal_B)
  }else # if no normal allele file, only tumor allele file
  {
    af=tumortable$Signal_A/(tumortable$Signal_A+tumortable$Signal_B)
  }
  idxprobes=which(af>=0.2 & af<=0.8) #set the threshold
  tumortable=merge(tumortable[idxprobes,],snp_anno)
  dh=abs(tumortable$Signal_A/(tumortable$Signal_A+tumortable$Signal_B)-0.5)
  
  tumortable=cbind.data.frame(tumortable,dh=dh)
  tumortable$Physical.Position=as.integer(tumortable$Physical.Position)
  
  chrs = paste0('chr',c(1:22,'X','Y'))
  CNA.obj = CNA(tumortable$dh, 
                tumortable$Chromosome, 
                tumortable$Physical.Position, 
                data.type = "logratio",
                sampleid=mysample)
  smoothed.CNA.obj = smooth.CNA(CNA.obj)
  
  segsall=c()
  #segment on each chromsome
  for (j in 1:length(chrs))
  {
    chr= chrs[j]
    if (grepl('chr',chr)) chr=substr(chr,4,nchar(chr))
    CNAchr.smoothed=smoothed.CNA.obj[smoothed.CNA.obj[,1]==chr,]
    if (nrow(CNAchr.smoothed)>0)
    {
      #Do the segment:
      segs <- segment(CNAchr.smoothed,alpha=0.1)
      segsall=rbind(segsall,segs$output)
    }
  }
  #write segment results to tmp folder
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.cbs.segment.txt")
  #write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  dh_probe_file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.txt")
  write.table(tumortable,file=dh_probe_file,row.names=F,col.names=T,sep="\t",quote=F)
  return(nrow(tumortable))
}
mpi.bcast.Robj2slave(CBS_allele)
res=mpi.parSapply(X=1:length(allsamples),FUN=CBS_allele,job.num=njob)

#After parallel computing, read the segment results in tmp folder, and combine them into a segment file
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber.allele.cbs.txt")
for (mysample in allsamples)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.cbs.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#once the segment file was written, we can start from here by simply reading the above results.
segsall=read.table(file=output,header=T,sep="\t")

copynumber_allele=readcnasegs1(hg19=T,segments=segsall)
colnames(copynumber_allele)=gsub(".","-",colnames(copynumber_allele),fixed=T)
data_copynumber_allele=cbind.data.frame(sample=colnames(copynumber_allele),t(copynumber_allele))
data_copynumber_allele[,"sample"]=as.character(data_copynumber_allele[,"sample"])
data_copynumber_allele=merge(tcga_platinumclass,data_copynumber_allele,by="sample")
rownames(data_copynumber_allele)=data_copynumber_allele[,"sample"]
data_copynumber_allele=data_copynumber_allele[,-which(colnames(data_copynumber_allele)=="sample")]

#compute the p-values:
pvalues_copynumber_allele=sapply(2:ncol(data_copynumber_allele),function(i){
  fm=glm(as.factor(data_copynumber_allele[,"platinumclass"])~as.numeric(data_copynumber_allele[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=rownames(copynumber_allele)[i-1]
  return(res)
})

manhattan_allele=draw_manhattan(pvalues=pvalues_copynumber_allele,keepres=T,maxy=5)
#check peaks
idx=which(manhattan_allele$chromosome==1)
manhattan_allele1=manhattan_allele[idx,]
which.min(manhattan_allele1$pvalue)

idx=which(manhattan_allele$chromosome==6)
manhattan_allele1=manhattan_allele[idx,]
which.min(manhattan_allele1$pvalue)

smallfdrgenes=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)
smallfdrgenes_pvalues=sapply(smallfdrgenes$gene,function(gene){
  idx=which(names(pvalues_copynumber_allele)==gene)
  if (length(idx)>0)
  {
    res=pvalues_copynumber_allele[[idx]]
  }else
  {
    res=NA
    #
  }
  #names(res)=gene
  return(res)
})
smallfdrgenes_pvalues=cbind.data.frame(gene=names(smallfdrgenes_pvalues),chr=smallfdrgenes$chr,start=smallfdrgenes$start,pvalues=smallfdrgenes_pvalues)

which.min(pvalues_copynumber_allele)
#HPGDS 
#5994
hist(as.numeric(copynumber_allele[5994,]))
plot(as.numeric(copynumber_allele[5994,]))
colnames(data_copynumber_allele)[5995]
sum(data_copynumber_allele[,5995]<0.2 & data_copynumber_allele[,1]=="Sensitive")
#107
sum(data_copynumber_allele[,5995]<0.2 & data_copynumber_allele[,1]=="Resistant")
#69
sum(data_copynumber_allele[,5995]>0.2 & data_copynumber_allele[,1]=="Sensitive")
#157
sum(data_copynumber_allele[,5995]>0.2 & data_copynumber_allele[,1]=="Resistant")
#38

#draw plot using total copynumber
pvalues_copynumber=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
genes=rownames(pvalues_copynumber)
pvalues_copynumber=pvalues_copynumber[,1]
names(pvalues_copynumber)=genes
draw_manhattan(pvalues = pvalues_copynumber)

which(names(pvalues_copynumber)=="HPGDS")
#[1] 6000

range(as.numeric(data_copynumber[1:nrow(data_copynumber),6001]),na.rm=T)
plot(c(0,0),xlim=c(-1.35,1.5),ylim=c(0,0.5),type="n",xlab="tcn",ylab="dh")
for (i in 1:nrow(data_copynumber))
{
  mysample=rownames(data_copynumber)[i]
  mycol="blue"
  if (data_copynumber[i,1]=="Resistant")
  {
    mycol="red"
  }
  idx=which(rownames(data_copynumber_allele)==mysample)
  points(x=data_copynumber[i,6001],y=data_copynumber_allele[idx,5995],col=mycol)
}
legend("topright",legend=c("Res","Sen"),col=c("red","blue"),pch=1)

x=data_copynumber[,6001]
z=data_copynumber[,1]
y=c()
for (i in 1:nrow(data_copynumber))
{
  mysample=rownames(data_copynumber)[i]
  idx=which(rownames(data_copynumber_allele)==mysample)
  y=c(y,data_copynumber_allele[idx,5995])
}

#compute 2-df pvalues--------------

two_df_pvalues=function(x,y,z)
{
  fit1 <- glm(z~x,family=binomial,x=T,y=T)
  fit2 <- glm(z~y,family=binomial,x=T,y=T)
  score1 <- fit1$x*(fit1$y-fit1$fitted)
  score2 <- fit2$x*(fit2$y-fit2$fitted)
  score <- cbind(score1,score2)
  bread <- matrix(0,4,4)
  bread[1:2,1:2] <- t(fit1$x) %*% (fit1$x*(1-fit1$fitted)*fit1$fitted)
  bread[3:4,3:4] <- t(fit2$x) %*% (fit2$x*(1-fit2$fitted)*fit2$fitted)
  varmat <- solve(bread) %*% (t(score)%*% score) %*% solve(bread)
  bb <- c(fit1$coef[2] ,fit2$coef[2])
  bbmat <- varmat[c(2,4),c(2,4)]
  tt <- drop(bb %*% solve(bbmat) %*% bb)
  pval <- 1-pchisq(tt,2)
}
two_df_pvalues_data=function(data1=data_copynumber,data2=data_copynumber_allele)
{
  #make sure the two datasets have the same order in row and columns
  commongenes=intersect(colnames(data1)[2:ncol(data1)],colnames(data2)[2:ncol(data2)])
  commonsamples=intersect(rownames(data1)[1:nrow(data1)],rownames(data2)[1:nrow(data2)])
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data1)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data1)==mysample)
  })
  #output
  z=data1[idx_row,1]
  data1=data1[idx_row,idx_col]
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data2)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data2)==mysample)
  })
  data2=data2[idx_row,idx_col]
  pvalues=sapply(1:ncol(data1),function(i){
    x=data1[,i]
    y=data2[,i]
    res=two_df_pvalues(x=x,y=y,z=z)
  })
  names(pvalues)=colnames(data1)
  ##for some cases when x or y are mostly 0, which causes the problem of p=0
  idx_pvalue0=which(pvalues==0)
  if (length(idx_pvalue0)>0) pvalues=pvalues[-idx_pvalue0]
  
  return(pvalues)
}

pvalues_2df=two_df_pvalues_data()
draw_manhattan(pvalues=pvalues_2df,maxy =7)
which.min(pvalues_2df)
pvalues_2df_tangent_gistic=two_df_pvalues_data(data1=data_copynumber_tangent_gistic)
draw_manhattan(pvalues=pvalues_2df_tangent_gistic,maxy=7)
which.min(pvalues_2df_tangent_gistic)
pvalues_2df_tangent=two_df_pvalues_data(data1=data_copynumber_tangent)
draw_manhattan(pvalues=pvalues_2df_tangent,maxy=7)
which.min(pvalues_2df_tangent)
load("../data/platinum_classificationdata_stringent_tangent.RData")
pvalues_2df_tangent_rmcnv=two_df_pvalues_data(data1=data_copynumber_tangent_rmcnv)
draw_manhattan(pvalues=pvalues_2df_tangent_rmcnv,maxy=7)
which.min(pvalues_2df_tangent_rmcnv)

#calculate correlation between tcn and dh
cor_tcn_dh_data=function(data1=data_copynumber_tangent,data2=data_copynumber_allele)
{
  #make sure the two datasets have the same order in row and columns
  commongenes=intersect(colnames(data1)[2:ncol(data1)],colnames(data2)[2:ncol(data2)])
  commonsamples=intersect(rownames(data1)[1:nrow(data1)],rownames(data2)[1:nrow(data2)])
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data1)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data1)==mysample)
  })
  data1=data1[idx_row,idx_col]
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data2)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data2)==mysample)
  })
  data2=data2[idx_row,idx_col]
  corvar=sapply(1:ncol(data1),function(i){
    x=data1[,i]
    y=data2[,i]
    res=cor(x,y,use="complete")
  })
  names(corvar)=colnames(data1)
  ##for some cases when x or y are mostly 0, which causes the problem of p=0
  return(corvar)
}
df_pvalues_2df=data.frame(gene=names(pvalues_2df),pvalues=pvalues_2df)
df_corvar=data.frame(gene=names(corvar),cor=corvar)
df_pvalues_2df=merge(df_pvalues_2df,df_corvar)

tmp=draw_manhattan(pvalues=corvar,keepres = T,logscale=F )

#check heterogeneous probes using normals and tumors
checkhtprobes=function(idxnormal)
{
  res=rep(NA,length(idxnormal))
  for (i in 1:length(idxnormal))
  {
    if (!is.na(idxnormal[i]))
    {
      mysample=allsamples[i]
      tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[i]],"/",filelist$filename[idxtumor[i]]),
                            skip = 1, sep="\t",header=T,stringsAsFactors = F)
      tumortable$Signal_A[which(tumortable$Signal_A<0)]=0
      tumortable$Signal_B[which(tumortable$Signal_B<0)]=0
      af_tumor=tumortable$Signal_A/(tumortable$Signal_A+tumortable$Signal_B)
      normaltable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxnormal[i]],"/",filelist$filename[idxnormal[i]]),
                             skip = 1, sep="\t",header=T,stringsAsFactors = F)
      normaltable$Signal_A[which(normaltable$Signal_A<0)]=0
      normaltable$Signal_B[which(normaltable$Signal_B<0)]=0
      af_normal=normaltable$Signal_A/(normaltable$Signal_A+normaltable$Signal_B)
      if (sum(normaltable$Composite.Element.REF==tumortable$Composite.Element.REF)<nrow(normaltable))
      {
        warning(paste0("probes not consistent for ",i))
      }
      cutoff=0.2
      idxprobes_normal=which(af_normal>=cutoff & af_normal<=1-cutoff)
      idxprobes_tumor=which(af_tumor>=cutoff & af_tumor<=1-cutoff)
      res[i]=length(intersect(idxprobes_normal,idxprobes_tumor))/length(idxprobes_normal)
    }
  }
  return(res)
}
save()