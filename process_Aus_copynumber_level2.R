#!/usr/bin/env Rscript

#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
#I use parallel computing when process the raw copy ratio data and apply CBS segment.

#include the function to draw manhattan plot:
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")

#process the raw data----------------------------
#get primary tumor samples
#section1: relate to CBS segment start
samplefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/sample.OV-AU.1474317952273.tsv"
allsamples=read.table(file=samplefile,header=T,sep="\t",stringsAsFactors=F)
primarysamples=allsamples[allsamples$specimen_type=="Primary tumour - solid tissue",]
primarysamples=primarysamples[primarysamples$sequencing_strategy=="non-NGS",]
primarysamples=primarysamples[primarysamples$repository=="GEO",]
uniq_primarysamples=unique(primarysamples$submitted_sample_id)

#ID convert: they used different IDs (GSM,AOCS,GEO) for the same set of samples. We need to know the mapping among these IDs.
GSM_AOCS_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSM_AOCS_ID.txt",sep="\t",header=F,stringsAsFactors = F)
AOCS_GEO_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSE65821_GEO_to_DCC_IDs.txt",header=T,sep="\t",stringsAsFactors = F)
#section1: relate to CBS segment end

# samples are from two platforms, these are the annotations:
omni_1.0=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni2.5-8v1_C_Gene_Annotation.txt",header=T,sep="\t",stringsAsFactors=F)
omni_1.1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni25M-8v1-1_B.annotated.txt",header=T,sep="\t",fill=T,stringsAsFactors=F)
#how many common probes on these two platfomrs:
commonprobes=intersect(omni_1.0$Name,omni_1.1$Name)
length(commonprobes)
#[1] 2338047

###---which probes should be used------------------------------------------------------------------------
#find out probeset with GC score >0.7 across all the samples. The same GC score filter was used by the Aus group
#We will have different probes for different samples if we apply the filter. We want to find the common ones, which satisfy the filter for all samples.
#parallel computing part:
njob=100
require(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
mpi.bcast.Robj2slave(commonprobes)
#apply filter to remove pobres.
removedprobes=function(myfile)
{
  #read each level 2 data
  mytable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                     stringsAsFactors=F)
  print(nrow(mytable))
  mytable=mytable[mytable$SNP.Name %in%commonprobes,]
  mytable=mytable[! is.na(mytable$GC.Score),]
  #keep results in tmp folder
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".removedpobes.txt")
  #filter:
  res=mytable$SNP.Name[which(mytable$GC.Score<=0.7)]
  res1=data.frame(matrix(NA,nrow=length(res)))
  colnames(res1)="Name"
  res1[,1]=res
  write.table(res1,file=output,row.names=F,quote=F)
  return(0)
}
mpi.bcast.Robj2slave(removedprobes)

#the folder contains copy ratio data for all samples, we only need to get those for primary tumor samples:
#section2:relate to CBS segment start
allfiles=list.files("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/")
primarysamplefiles=c()
for (mysample in uniq_primarysamples)
{
  idx=which(AOCS_GEO_ID$submitted_sample_id==mysample)
  if (length(idx)==1)
  {
    mygeo=AOCS_GEO_ID$GEO_sample_id[idx]
    idx1=which(GSM_AOCS_ID[,2]==mygeo)
    mygsm=GSM_AOCS_ID[idx1,1]
    mygsm=gsub(" ","",mygsm,fixed=T)
    idx2=which(grepl(paste0(mygsm,"\\w*.txt"),allfiles,perl=T)==T)
    myfile=allfiles[idx2]
    primarysamplefiles=c(primarysamplefiles,myfile)
  }
}
#section2:relate to CBS segment end
#write the probes to be removed for each sample in tmp folder
res=mpi.parSapply(X=primarysamplefiles,FUN=removedprobes,job.num=njob)

#read the output of above parallel computing, and find all the probes need to be removed 
probes2remove=c()
for (myfile in primarysamplefiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".removedpobes.txt"),header=T)
  probes2remove=c(probes2remove,tmp$Name)
  probes2remove=unique(probes2remove)
}
length(probes2remove)
#[1] 538087

#the pobes to keep:
probes2keep=commonprobes[!commonprobes %in% probes2remove]
length(probes2keep)
#[1] 2338047
#remove probes in cnvs:
cnvfile="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/examplefiles/SNP6.merged.151117.hg19.CNV.txt"
cnvs=read.table(file=cnvfile,header=T,sep="\t")
cnvs$Chromosome=gsub("23","X",cnvs$Chromosome)
cnvs$Chromosome=gsub("24","Y",cnvs$Chromosome)
library(GenomicRanges)
mpi.bcast.cmd(library(GenomicRanges))
gr_cnvs=GRanges(seqnames=cnvs$Chromosome,ranges=IRanges(start=cnvs$Start,end=cnvs$End))
df_probes2keep=data.frame(Name=probes2keep)
df_probes2keep=merge(df_probes2keep,omni_1.1[,1:3])
gr_probes2keep=GRanges(seqnames=df_probes2keep$Chr,ranges=IRanges(start=df_probes2keep$MapInfo,width=1))
mpi.bcast.Robj2slave(gr_cnvs)
mpi.bcast.Robj2slave(gr_probes2keep)
probes_in_cnv1=function(i)
{
  olap=subsetByOverlaps(gr_probes2keep[i,],gr_cnvs)
  res=F
  if (length(olap)>0) res=T
  return(res)
}
mpi.bcast.Robj2slave(probes_in_cnv1)
idx_incnv=mpi.parSapply(X=1:length(gr_probes2keep),FUN=probes_in_cnv1,job.num=njob)
probes2keep_rmcnv=probes2keep[!idx_incnv]



#save(probes2keep,probes2keep_rmcnv,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/probes2keep.RData")


#CBS segmentation----------------------------------------------------------------
library(DNAcopy)
mpi.bcast.cmd(library(DNAcopy))
chrs = paste0('chr',c(1:22,'X','Y'))
#only work on probes meet the GC filter
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/probes2keep.RData")

mpi.bcast.Robj2slave(probes2keep)
mpi.bcast.Robj2slave(probes2keep_rmcnv)
#parallel computing part:
#function to apply CBS on each level 2 copy ratio file:
#To see how it works, you can run CBS(myfile=primarysamplefiles[1])
CBS=function(myfile,type="tcn",removecnv=F)
{
  #to find sample id:
  GSM_AOCS_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSM_AOCS_ID.txt",sep="\t",header=F,stringsAsFactors = F)
  AOCS_GEO_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSE65821_GEO_to_DCC_IDs.txt",header=T,sep="\t",stringsAsFactors = F)
  gsm=unlist(strsplit(myfile,"_",fixed=T))[1]
  GSM_AOCS_ID$V1=gsub(" ","",GSM_AOCS_ID$V1,fixed=T)
  idx1=which(GSM_AOCS_ID$V1==gsm)
  aocs=GSM_AOCS_ID$V2[idx1]
  idx2=which(AOCS_GEO_ID$GEO_sample_id==aocs)
  #sample id:
  mysample=AOCS_GEO_ID$submitted_sample_id[idx2]
  
  chrs = paste0('chr',c(1:22,'X','Y'))
  #read the log ratio data:
  mytable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                     stringsAsFactors=F)
  #Apply CBS on the data: only keep a few columns, and probes after filtering
  if (removecnv==T)
  {
    mytable=mytable[mytable$SNP.Name %in% probes2keep_rmcnv,c("SNP.Name","Chr","Position" ,"B.Allele.Freq","Log.R.Ratio")]
  }else
  {
    mytable=mytable[mytable$SNP.Name %in% probes2keep,c("SNP.Name","Chr","Position" ,"B.Allele.Freq","Log.R.Ratio")]
  }
  
  if (type=="tcn") #total copy number
  {
    CNA.obj = CNA(mytable$Log.R.Ratio, 
                  mytable$Chr, 
                  mytable$Position, 
                  data.type = "logratio",
                  sampleid=mysample)
  }else #DH
  {
    cutoff=0.2
    mytable=mytable[mytable$B.Allele.Freq>=cutoff & mytable$B.Allele.Freq<=1-cutoff,]
    CNA.obj = CNA(abs(0.5-mytable$B.Allele.Freq), 
                  mytable$Chr, 
                  mytable$Position, 
                  data.type = "logratio",
                  sampleid=mysample)
    
  }
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
      if (type=="tcn")
      {
        segs <- segment(CNAchr.smoothed,alpha=0.0001)
      }else
      {
        segs <- segment(CNAchr.smoothed,alpha=0.05)
      }
      segsall=rbind(segsall,segs$output)
    }
  }
  #write segment results to tmp folder
  if (removecnv==T)
  {
    if (type=="tcn")
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.rmcnv.segment.txt")
    }else
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.rmcnv.segment.txt")
    }
  }else
  {
    if (type=="tcn")
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.segment.txt")
    }else
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.segment.txt")
      dh_probe_file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.txt")
      write.table(mytable,file=dh_probe_file,row.names=F,col.names=T,sep="\t",quote=F)
    }
  }
  write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  
  return(nrow(mytable))
}

mpi.bcast.Robj2slave(CBS)
#Do segment on all samples, write segment results to tmp folder
res=mpi.parSapply(X=primarysamplefiles,removecnv=T,FUN=CBS,job.num=njob)
res=mpi.parSapply(X=primarysamplefiles,type="allele",FUN=CBS,job.num=njob)
res=mpi.parSapply(X=primarysamplefiles,removecnv=T,type="allele",FUN=CBS,job.num=njob)


#quit parallel computing
mpi.close.Rslaves()
mpi.quit()

#After parallel computing, read the segment results in tmp folder, and combine them into a segment file
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.cbs.txt")
for (myfile in primarysamplefiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#once the segment file was written, we can start from here by simply reading the above results.
segsall=read.table(file=output,header=T,sep="\t")

#form the marker file for GISTIC------------------
res=data.frame(matrix(NA,nrow=length(probes2keep),ncol=3))
colnames(res)=c('unitName','chromosome','position')
res[,1]=probes2keep
df_probes2keep=data.frame(Name=probes2keep)
tmp=merge(df_probes2keep,omni_1.0)
colnames(tmp)[2]="chr"
colnames(tmp)[3]="start"
res=sortgenetable(tmp)
res=res[,1:3]
colnames(res)=c('unitName','chromosome','position')
output="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.markers.txt"
#write.table(res,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#form copynumber data at gene level based on segment result from DNAcopy-------------------
readcnasegs1=function(hg19=T,segments,opt=1)
{
  library(CNTools)
  colnames(segments)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  dim(segments)
  length(unique(segments$ID))
  segments=segments[!is.na(segments$seg.mean),]
  cnseg_snp6=CNSeg(segments)
  if (opt==1)
  {
#     hg19_snp6=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/copynumber_geneposition_biomart.txt",header=T,sep="\t")
#     geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
#     hg19=merge(hg19_snp6,geneposition,by="gene")
#     hg19=hg19[,c(1,5,6,7)]
    #use the one from firehose
    hg19=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
  }
  if (opt==2) #use the genesets defined by CNTools
  {
    load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")
    hg19=geneInfo_hg19[,c(4,1,2,3)]
  }
  colnames(hg19)=c("gene","chr","start","end")
  hg19=sortgenetable(hg19)
  
  hg19$chr=gsub("23","X",hg19$chr)
  hg19$chr=gsub("24","Y",hg19$chr)
  idxdup=duplicated(hg19$gene)
  hg19=hg19[! idxdup,]
  #change the column name of chr
  colnames(hg19)[2]="chrom"
  nstart=5
  snp6_max<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "max")
  snp6_max=rs(snp6_max)
  snp6_min<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "min")
  snp6_min=rs(snp6_min)
  snp6=snp6_max
  
  #pick the extreme one
  for (i in nstart:ncol(snp6))
  {
    pickmin=which(abs(snp6_min[,i])-abs(snp6_max[,i])>0)
    if (length(pickmin)>0) snp6[pickmin,i]=snp6_min[pickmin,i]
  }
  
  rownames(snp6)=snp6$gene
  snp6=snp6[,nstart:ncol(snp6)]
  
  #remove all 0 genes
  idxkeep=rep(T,nrow(snp6))
  idxkeep=sapply(1:nrow(snp6),function(i){
    res=F
    num0=sum(snp6[i,]==0)
    if (num0<ncol(snp6))
    {
      res=T
    }
    return(res)
  })
  snp6=snp6[idxkeep,]
  return(snp6)
}

#the data was generated based on mean value of probes
aocs_mean_copynumber=readcnasegs1(hg19=T,segments=segsall)
colnames(aocs_mean_copynumber)=gsub(".","-",colnames(aocs_mean_copynumber),fixed=T)


#check p-values--------------------------------------
#read resistance def
outcomedeffile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/nature14410-Supplementary Table 5.xlsx"
library(gdata)
outcomedef=read.xls(outcomedeffile,1,skip=1)
#the last two lines are comments
outcomedef=outcomedef[1:(nrow(outcomedef)-2),]
outcomedef$Chemotherapy.response=as.character(outcomedef$Chemotherapy.response)
outcomedef$Case.ID=as.character(outcomedef$Case.ID)
outcomedef$Case.ID=gsub("_","-",outcomedef$Case.ID,fixed=T)
outcomedef=outcomedef[outcomedef$Sample.type=="7PrimaryTumour",]
sample_platinumstatus=data.frame(matrix(NA,nrow=nrow(outcomedef),ncol=3))
colnames(sample_platinumstatus)=c("patient","sample","status")
sample_platinumstatus$patient=outcomedef$Case.ID
sample_platinumstatus$status=outcomedef$Chemotherapy.response
for (i in 1:nrow(sample_platinumstatus))
{
  idx=which(grepl(sample_platinumstatus$patient[i],uniq_primarysegsample)==T)
  if (length(idx)>0)
  {
    sample_platinumstatus$sample[i]=uniq_primarysegsample[idx]
  }
}
#make the sample order of resistance (output) the same as the one used in copynumber data
aocs_mean_platinumclass=rep(NA,ncol(aocs_mean_copynumber))
for (i in 1:ncol(aocs_mean_copynumber))
{
  idx=which(sample_platinumstatus$sample==colnames(aocs_mean_copynumber)[i])
  if (length(idx)>0)
  {
    aocs_mean_platinumclass[i]=sample_platinumstatus$status[idx]
  }
}

aocs_mean_platinumclass=gsub("refractory","resistant",aocs_mean_platinumclass)
#combine data and output:
data_aocs_mean_copynumber=cbind(platinumclass=aocs_mean_platinumclass,as.data.frame(t(aocs_mean_copynumber)))

#save(data_aocs_copynumber,data_aocs_mean_copynumber,data_aocs_mean_copynumber_rmcnv,data_aocs_mean_copynumber_allele,data_aocs_mean_copynumber_allele_rmcnv,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_copynumber.RData")

#remove cnvs
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.rmcnv.cbs.txt")
for (myfile in primarysamplefiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.rmcnv.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#once the segment file was written, we can start from here by simply reading the above results.
segsall=read.table(file=output,header=T,sep="\t")
aocs_mean_copynumber_rmcnv=readcnasegs1(hg19=T,segments=segsall)
colnames(aocs_mean_copynumber_rmcnv)=gsub(".","-",colnames(aocs_mean_copynumber_rmcnv),fixed=T)
#combine data and output:
data_aocs_mean_copynumber_rmcnv=cbind(platinumclass=aocs_mean_platinumclass,as.data.frame(t(aocs_mean_copynumber_rmcnv)))

#work with the allele data
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.allele.cbs.txt")
for (myfile in primarysamplefiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#once the segment file was written, we can start from here by simply reading the above results.
segsall=read.table(file=output,header=T,sep="\t")
aocs_mean_copynumber_allele=readcnasegs1(hg19=T,segments=segsall)
colnames(aocs_mean_copynumber_allele)=gsub(".","-",colnames(aocs_mean_copynumber_allele),fixed=T)
sum(colnames(aocs_mean_copynumber_allele)==aocs_platinumclass$sample)
#combine data and output:
data_aocs_mean_copynumber_allele=cbind(platinumclass=aocs_mean_platinumclass,as.data.frame(t(aocs_mean_copynumber_allele)))

#remove cnv
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.allele.rmcnv.cbs.txt")
for (myfile in primarysamplefiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.rmcnv.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#once the segment file was written, we can start from here by simply reading the above results.
segsall=read.table(file=output,header=T,sep="\t")
aocs_mean_copynumber_allele_rmcnv=readcnasegs1(hg19=T,segments=segsall)
colnames(aocs_mean_copynumber_allele_rmcnv)=gsub(".","-",colnames(aocs_mean_copynumber_allele_rmcnv),fixed=T)
sum(colnames(aocs_mean_copynumber_allele_rmcnv)==aocs_platinumclass$sample)
#combine data and output:
data_aocs_mean_copynumber_allele_rmcnv=cbind(platinumclass=aocs_mean_platinumclass,as.data.frame(t(aocs_mean_copynumber_allele_rmcnv)))



#compute 2df-pvalue
pvalues_aocs_mean_2df=two_df_pvalues_data(data1=data_aocs_mean_copynumber,data2=data_aocs_mean_copynumber_allele)
draw_manhattan(pvalues=pvalues_aocs_mean_2df)

pvalues_aocs_2df=two_df_pvalues_data(data1=data_aocs_copynumber,data2=data_aocs_mean_copynumber_allele)
draw_manhattan(pvalues=pvalues_aocs_2df)

aocs_mean_copynumber_pvalues=sapply(2:ncol(data_aocs_mean_copynumber),function(i){
  fm=glm(as.factor(data_aocs_mean_copynumber[,"platinumclass"])~as.numeric(data_aocs_mean_copynumber[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_aocs_mean_copynumber)[i-1]
  return(res)
})
draw_manhattan(pvalues=aocs_mean_copynumber_pvalues[!is.na(aocs_mean_copynumber_pvalues)])

aocs_copynumber_pvalues=sapply(2:ncol(data_aocs_copynumber),function(i){
  fm=glm(as.factor(data_aocs_copynumber[,"platinumclass"])~as.numeric(data_aocs_copynumber[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_aocs_copynumber)[i-1]
  return(res)
})
draw_manhattan(pvalues=aocs_copynumber_pvalues[!is.na(aocs_copynumber_pvalues)])

pvalues_aocs_mean_copynumber_allele=sapply(2:ncol(data_aocs_mean_copynumber_allele),function(i){
  fm=glm(as.factor(data_aocs_mean_copynumber_allele[,"platinumclass"])~as.numeric(data_aocs_mean_copynumber_allele[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_aocs_mean_copynumber_allele)[i-1]
  return(res)
})
draw_manhattan(pvalues=pvalues_aocs_mean_copynumber_allele[!is.na(pvalues_aocs_mean_copynumber_allele)])

aocs_mean_copynumber_rmcnv_pvalues=sapply(2:ncol(data_aocs_mean_copynumber_rmcnv),function(i){
  fm=glm(as.factor(data_aocs_mean_copynumber_rmcnv[,"platinumclass"])~as.numeric(data_aocs_mean_copynumber_rmcnv[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_aocs_mean_copynumber_rmcnv)[i-1]
  return(res)
})
draw_manhattan(pvalues=aocs_mean_copynumber_rmcnv_pvalues[!is.na(aocs_mean_copynumber_rmcnv_pvalues)])

aocs_mean_copynumber_allele_rmcnv_pvalues=sapply(2:ncol(data_aocs_mean_copynumber_allele_rmcnv),function(i){
  fm=glm(as.factor(data_aocs_mean_copynumber_allele_rmcnv[,"platinumclass"])~as.numeric(data_aocs_mean_copynumber_allele_rmcnv[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_aocs_mean_copynumber_allele_rmcnv)[i-1]
  return(res)
})
draw_manhattan(pvalues=aocs_mean_copynumber_allele_rmcnv_pvalues[!is.na(aocs_mean_copynumber_allele_rmcnv_pvalues)])

pvalues_aocs_mean_rmcnv_2df=two_df_pvalues_data(data1=data_aocs_mean_copynumber_rmcnv,data2=data_aocs_mean_copynumber_allele_rmcnv)
draw_manhattan(pvalues=pvalues_aocs_mean_rmcnv_2df)

aocs_cor_tcn_dh=cor_tcn_dh_data(data1=data_aocs_mean_copynumber,data2=data_aocs_mean_copynumber_allele)
draw_manhattan(pvalues=aocs_cor_tcn_dh,logscale = F)

aocs_rmcnv_cor_tcn_dh=cor_tcn_dh_data(data1=data_aocs_mean_copynumber_rmcnv,data2=data_aocs_mean_copynumber_allele_rmcnv)
draw_manhattan(pvalues=aocs_rmcnv_cor_tcn_dh,logscale=F)

#the generated gistic data:
copynumber_aocs_gistic=read.table(file="../../Tools/GISTIC/AOCS_rx0_conf0.99_armpeel1_brlen0.7_broad1/all_data_by_genes.txt",header=T,
                                     sep="\t",quote="")
rownames(copynumber_aocs_gistic)=copynumber_aocs_gistic$Gene.Symbol
copynumber_aocs_gistic=copynumber_aocs_gistic[,4:ncol(copynumber_aocs_gistic)]
colnames(copynumber_aocs_gistic)=gsub(".","-",colnames(copynumber_aocs_gistic),fixed=T)

data_copynumber_aocs_gistic=cbind.data.frame(sample=colnames(copynumber_aocs_gistic),t(copynumber_aocs_gistic))
data_copynumber_aocs_gistic[,"sample"]=as.character(data_copynumber_aocs_gistic[,"sample"])
aocs_platinumclass=rep(NA,ncol(aocs_copynumber))
for (i in 1:ncol(aocs_copynumber))
{
  idx=which(sample_platinumstatus$sample==colnames(aocs_copynumber)[i])
  if (length(idx)>0)
  {
    aocs_platinumclass[i]=sample_platinumstatus$status[idx]
  }
}
aocs_platinumclass=gsub("refractory","resistant",aocs_platinumclass)
aocs_platinumclass=data.frame(sample=colnames(aocs_copynumber),platinumclass=aocs_platinumclass)
data_copynumber_aocs_gistic=merge(aocs_platinumclass,data_copynumber_aocs_gistic,by="sample")

rownames(data_copynumber_aocs_gistic)=data_copynumber_aocs_gistic[,"sample"]
data_copynumber_aocs_gistic=data_copynumber_aocs_gistic[,-which(colnames(data_copynumber_aocs_gistic)=="sample")]

copynumber_aocs_gistic_pvalues=sapply(2:ncol(data_copynumber_aocs_gistic),function(i){
  fm=glm(as.factor(data_copynumber_aocs_gistic[,"platinumclass"])~as.numeric(as.character(data_copynumber_aocs_gistic[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_aocs_gistic)[i-1]
  return(res)
})
test=copynumber_aocs_gistic_pvalues[!is.na(copynumber_aocs_gistic_pvalues)]
draw_manhattan(test)

