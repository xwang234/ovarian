#!/usr/bin/env Rscript
# revise the code by James Dai July 24, 2017
#resistance def:
rm(list=ls())
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")

snp_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno <- snp_anno[,1:4]

allsamples=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_samplelist.txt",stringsAsFactors = F)
allsamples=allsamples[,1]
failedsamples=c("TCGA-13-1481","TCGA-23-1122","TCGA-25-1323","TCGA-29-2414","TCGA-61-2111","TCGA-23-1119")

allsamples=allsamples[!allsamples %in% failedsamples]
tumorsamples=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_tumoronlysamplelist.txt",stringsAsFactors = F)
allsamples=c(allsamples,tumorsamples[,1])

LOHseg <- matrix(NA,length(allsamples),6)

for (i in 1:length(allsamples)){
  mysample=allsamples[i]
  
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.segment")
  dat <- read.table(output,header=T)
  dat$tcn=dat$nMajor+dat$nMinor
  dat$lohCall=F
  dat$cnloh=F

  for (j in 1:nrow(dat))
  {
    if (dat$nMajor[j]==0 | dat$nMinor[j]==0) dat$lohCall[j]=T
    if(dat$tcn[j]==2 & dat$lohCall[j]==T) dat$cnloh[j] = T 
  }
  LOHseg[i,1] <- sum(dat$endpos[dat$lohCall]-dat$startpos[dat$lohCall])/10^6
  LOHseg[i,2] <- length(dat$endpos[!is.na(dat$lohCall) & dat$lohCall])
  LOHseg[i,3] <- sum(as.numeric(dat$endpos[dat$tcn!=2]) -as.numeric(dat$startpos[dat$tcn!=2]))/10^6
  LOHseg[i,4] <- length(dat$endpos[dat$tcn!=2 & !dat$lohCall])
  
  LOHseg[i,5] <- sum(dat$endpos[dat$cnloh]-dat$startpos[dat$cnloh])/10^6
  LOHseg[i,6] <- length(dat$endpos[!is.na(dat$cnloh) & dat$cnloh])
  
}

loh <- data.frame(cbind(allsamples,LOHseg))
loh[,1] <- as.character(loh[,1])
loh <- loh[order(loh[,1]),]
names(loh) <- c("sampleID","loh_length","n_lohsegs","tcn_length","n_tcnsegs","cnloh_length","n_cnlohsegs")
resist <- data_copynumber_tangent_filtered[,1]
sampleID <- rownames(data_copynumber_tangent_filtered)
resist <- data.frame(cbind(sampleID,as.character(resist)))
names(resist) <- c("sampleID","resistance")
resist <- resist[order(resist$sampleID),]
resist <- merge(resist,loh,by="sampleID")
for (i in 3:8)  resist[,i] <- as.numeric(as.character(resist[,i]))

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_purity_ploidy.RData")
info <- info[order(info$sampleID),]
resist <- merge(resist,info,by="sampleID")
resist$avelohlength <- resist$loh_length/resist$n_lohsegs

summary(glm(I(resistance=="Sensitive")~loh_length,family=binomial,data=resist)) 
summary(glm(I(resistance=="Sensitive")~n_lohsegs,family=binomial,data=resist)) 
summary(glm(I(resistance=="Sensitive")~avelohlength,family=binomial,data=resist)) 
summary(glm(I(resistance=="Sensitive")~tcn_length,family=binomial,data=resist)) 
summary(glm(I(resistance=="Sensitive")~n_tcnsegs,family=binomial,data=resist)) 
summary(glm(I(resistance=="Sensitive")~cnloh_length,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~n_cnlohsegs,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~purity,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~I(purity>0.8),family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~ploidy,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~I(ploidy>=2),family=binomial,data=resist))



library(beeswarm)
postscript("~/SCHARPHOME/ScharpFile/Article/TCGA_ovarian/ASCAT_LOH_platinum.ps",width=8,height=6)
par(mfrow=c(1,1))
beeswarm(resist$loh_length/3000 ~resist$resistance,pch=19,cex=1.2,col=3:4,xlab="Platinum response", ylab="% genome has LOH")
bxplot(resist$loh_length/3000 ~resist$resistance,probs=0.5,col=2,add=T,lwd=3)
dev.off()

postscript("~/SCHARPHOME/ScharpFile/Article/TCGA_ovarian/ASCAT_LOHsegs_platinum.ps",width=8,height=6)
beeswarm(resist$n_lohsegs ~resist$resistance,pch=19,cex=1.2,col=3:4,xlab="Platinum response", ylab="number of LOH segments",ylim=c(0,250))
bxplot(resist$n_lohsegs ~resist$resistance,probs=0.5,col=2,add=T)
dev.off()

postscript("~/SCHARPHOME/ScharpFile/Article/TCGA_ovarian/ASCAT_ploidy_platinum.ps",width=8,height=6)
par(mfrow=c(1,1))
beeswarm(resist$ploidy ~resist$resistance,pch=19,cex=1.2,col=3:4,xlab="Platinum response", ylab="Ploidy")
abline(h=2,lty=2,lwd=3,col=2)
dev.off()



segsall=c()
for (mysample in allsamples)
{
  dat=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.segment"),header=T)
  names(dat)[1] <- "ID"
  dat$tcn=dat$nMajor+dat$nMinor
  dat$loh=0
  for (j in 1:nrow(dat))
  {
    if ((dat$nMajor[j]==0 | dat$nMinor[j]==0) & dat$nMajor[j]+dat$nMinor[j]>0) dat$loh[j]=1
  }
  
  segsall=rbind(segsall,dat[,c(1:4,8)])
}

segsall=c()
for (mysample in allsamples)
{
  dat=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.segment"),header=T)
  names(dat)[1] <- "ID"
  dat$tcn=dat$nMajor+dat$nMinor
  dat$loh=0
  for (j in 1:nrow(dat))
  {
    if ((dat$nMajor[j]==0 | dat$nMinor[j]==0) & dat$nMajor[j]+dat$nMinor[j]>0) dat$loh[j]=1
  }
  
  segsall=rbind(segsall,dat)
}
segsall_cn=segsall[,c(1:4,7)]
segsall_cn$tcn=segsall_cn$tcn-2
segsall_loh=segsall[,c(1:4,8)]
segsall_loh=segsall_loh[segsall_loh$loh==1,]
save(segsall_cn,segsall_loh,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_segment.RData")


sortgenetable=function(genetable)
{
  genetable$chr=as.character(genetable$chr)
  genetable$chr=gsub("chr","",genetable$chr)
  genetable$chr=gsub("23","X",genetable$chr)
  genetable$chr=gsub("24","Y",genetable$chr)
  if (class(genetable$start[1])=="factor")
  {
    genetable$start=as.numeric(as.character(genetable$start))
  }
  if (class(genetable$start[1])=="character")
  {
    genetable$start=as.numeric(genetable$start)
  }
  
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genetable)))
  for (chr in chrs)
  {
    tmptable=genetable[which(genetable$chr==chr),]
    tmptable=tmptable[order(tmptable$start),]
    res=rbind(res,tmptable)
  }
  return(res)
}


readcnasegs2=function(hg19=T,snp6copynumber,opt=1)
{
  library(CNTools)
  colnames(snp6copynumber)=c("ID","chrom","loc.start","loc.end","seg.mean")
  dim(snp6copynumber)
  length(unique(snp6copynumber$ID))
  snp6copynumber=snp6copynumber[!is.na(snp6copynumber$seg.mean),]
  cnseg_snp6=CNSeg(snp6copynumber)
  if (opt==1)
  {
    hg19_snp6=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/copynumber_geneposition_biomart.txt",header=T,sep="\t")
    geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
    hg19=merge(hg19_snp6,geneposition,by="gene")
    hg19=hg19[,c(1,5,6,7)]
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
    if (length(pickmin)>0) {
      cat(i,"..") 
      snp6[pickmin,i]=snp6_min[pickmin,i]
    }  
  }
  rownames(snp6)=snp6$gene
  snp6=snp6[,nstart:ncol(snp6)]
  return(snp6)
}
TCGA_gene_loh=readcnasegs2(hg19=T,snp6copynumber=segsall)
TCGA_gene_loh <- TCGA_gene_loh[,order(names(TCGA_gene_loh))]

outp <- matrix(NA,nrow(TCGA_gene_loh),1)

for (i in 1:nrow(TCGA_gene_loh)){
  if ((i %% 1000)==0) cat(i,"..")
  temp <- matrix(NA,ncol(TCGA_gene_loh),2)
  temp[,1] <- as.vector(names(TCGA_gene_loh))
  temp[,2] <- as.character(TCGA_gene_loh[i,])
  
  temp2 <- data.frame(temp)
  names(temp2) <- c("sampleID","loh")
  temp2$loh <- as.numeric(as.character(temp2$loh))
  
  if (mean(temp2$loh!=0)>0.05) {
    resist1 <- merge(resist,temp2,by="sampleID")
    fit2 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
    outp[i,1] <- summary(fit2)$coeff[2,4]
  }
}

lohpvalue=outp[,1]
names(lohpvalue)=rownames(TCGA_gene_loh)

save(lohpvalue,segsall,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_pvalue.Rdata")
source("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
draw_manhattan(pvalues=lohpvalue,ylab="log10(p)")

outp2=outp
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_pvalues_qvalues.RData")
lohpvalue1=outp[,5]
names(lohpvalue1)=rownames(outp)

lohpvalue=outp2[,1]
names(lohpvalue)=rownames(TCGA_gene_loh)

par(mfrow=c(2,1))
draw_manhattan(pvalues=lohpvalue,ylab="-log10(p)",main="ASCAT")
draw_manhattan(pvalues=lohpvalue1,ylab="-log10(p)",main="PSCBS")

#check purity and ploidy
info=data.frame(sampleID=allsamples,purity=0,ploidy=0)
for (i in 1:length(allsamples))
{
  mysample=allsamples[i]
  dat=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.res"),header=F)
  info$purity[i]=dat$V2[2]
  info$ploidy[i]=dat$V2[3]
}
save(info,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_purity_ploidy.RData")


resist2=merge(resist1,info,by="sampleID")
summary(glm(I(resistance=="Sensitive")~purity,family=binomial,data=resist2)) #p=0.264
summary(glm(I(resistance=="Sensitive")~ploidy,family=binomial,data=resist2)) #p=0.28115

#work on cytobands-----------------
library(GenomicRanges)
cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband$chrom=as.character(cytoband$chrom)
cytoband=cytoband[cytoband$chrom!="chrY",]
cytoband$chrom=gsub("chr","",cytoband$chrom)
gr_cytoband=GRanges(seqnames = cytoband$chrom,ranges=IRanges(start=cytoband$start,end=cytoband$end))
  
  
extractcytoband=function(seg=segsall,cutoff=0.25)
{
  allsamples=unique(as.character(seg$ID))
  res=data.frame(matrix(0,nrow=nrow(cytoband),ncol=length(allsamples)))
  rownames(res)=paste0(cytoband$chrom,"_",cytoband$name)
  colnames(res)=allsamples
  for (i in 1:length(allsamples))
  {
    if (i %%10 ==0) cat(i,"..")
    sampleseg=seg[seg$ID==allsamples[i] & seg$loh==1,]
    
    gr_sampleseg=GRanges(seqnames = sampleseg$chr,ranges=IRanges(start=sampleseg$startpos,end=sampleseg$endpos),loh=sampleseg$loh)
    tmp=countOverlaps(gr_cytoband,gr_sampleseg)
    idxs=which(tmp>0)
    for (j in idxs)
    {
      olap=intersect(gr_cytoband[j],gr_sampleseg)
      if (length(olap)>0)
      {
        if (sum(width(olap))>=cutoff*width(gr_cytoband[j]))
        {
          res[j,i]=1
        }
      }
    }
  }
}

extractcytoband_mpi=function(i,cutoff=0.25)
{
  allsamples=unique(as.character(seg$ID))
  res=data.frame(matrix(0,ncol=nrow(cytoband),nrow=1))
  colnames(res)=paste0(cytoband$chrom,"_",cytoband$name)
  rownames(res)=allsamples[i]
  sampleseg=seg[seg$ID==allsamples[i] & seg$loh==1,]
    
  gr_sampleseg=GRanges(seqnames = sampleseg$chr,ranges=IRanges(start=sampleseg$startpos,end=sampleseg$endpos),loh=sampleseg$loh)
  tmp=countOverlaps(gr_cytoband,gr_sampleseg)
  idxs=which(tmp>0)
  for (j in idxs)
  {
    olap=intersect(gr_cytoband[j],gr_sampleseg)
    if (length(olap)>0)
    {
        if (sum(width(olap))>=cutoff*width(gr_cytoband[j]))
        {
          res[1,j]=1
        }
    }
  }
  return(res)
}

#salloc -t 0-5 -n 49 mpirun -n 1 R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_pvalue.Rdata")
seg=segsall
mpi.bcast.Robj2slave(seg)
mpi.bcast.Robj2slave(cytoband)
mpi.bcast.Robj2slave(gr_cytoband)
mpi.bcast.cmd(library(GenomicRanges))

allsamples=unique(as.character(seg$ID))
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband_mpi,job.num=njobs)
res1=matrix(res,byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25.RData")
#TCGA_cytoband_loh=extractcytoband()
TCGA_cytoband_loh=res1
outp <- matrix(NA,nrow(TCGA_cytoband_loh),1)

for (i in 1:nrow(TCGA_cytoband_loh)){
  if ((i %% 100)==0) cat(i,"..")
  temp <- matrix(NA,ncol(TCGA_cytoband_loh),2)
  temp[,1] <- as.vector(colnames(TCGA_cytoband_loh))
  temp[,2] <- as.character(TCGA_cytoband_loh[i,])
  
  temp2 <- data.frame(temp)
  names(temp2) <- c("sampleID","loh")
  temp2$loh <- as.numeric(as.character(temp2$loh))
  
  if (mean(temp2$loh!=0)>0.05) {
    resist1 <- merge(resist,temp2,by="sampleID")
    fit2 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
    outp[i,1] <- summary(fit2)$coeff[2,4]
  }
}

cytobandres=cbind(cytoband,p=outp[,1])

save(cytobandres,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25_pvalue.RData")
