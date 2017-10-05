#!/usr/bin/env Rscript


#resistance def:
samplefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/sample.OV-AU.1474317952273.tsv"
allsamples=read.table(file=samplefile,header=T,sep="\t",stringsAsFactors=F)
primarysamples=allsamples[allsamples$specimen_type=="Primary tumour - solid tissue",]
primarysamples=primarysamples[primarysamples$sequencing_strategy=="non-NGS",]
primarysamples=primarysamples[primarysamples$repository=="GEO",]
uniq_primarysamples=unique(primarysamples$submitted_sample_id)

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
  idx=which(grepl(sample_platinumstatus$patient[i],uniq_primarysamples)==T)
  if (length(idx)>0)
  {
    sample_platinumstatus$sample[i]=uniq_primarysamples[idx]
  }
}
resist=sample_platinumstatus[,2:3]
names(resist) <- c("sampleID","resistance")

allsamples=resist$sampleID


LOHseg <- matrix(NA,length(allsamples),4)

for (i in 1:length(allsamples)){
  mysample=allsamples[i]
  
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.segment")
  dat <- read.table(output,header=T)
  dat$tcn=dat$nMajor+dat$nMinor
  dat$lohCall=F
  for (j in 1:nrow(dat))
  {
    if ((dat$nMajor[j]==0 | dat$nMinor[j]==0) & dat$nMajor[j]+dat$nMinor[j]>0) dat$lohCall[j]=T
  }
  LOHseg[i,1] <- sum(dat$endpos[dat$lohCall]-dat$startpos[dat$lohCall])/10^6
  LOHseg[i,2] <- length(dat$endpos[!is.na(dat$lohCall) & dat$lohCall])
  LOHseg[i,3] <- sum(as.numeric(dat$endpos[dat$tcn!=2 & !dat$lohCall]) -as.numeric(dat$startpos[dat$tcn!=2 & !dat$lohCall]))/10^6
  LOHseg[i,4] <- length(dat$endpos[dat$tcn!=2 & !dat$lohCall])
}
loh <- data.frame(cbind(allsamples,LOHseg))
loh[,1] <- as.character(loh[,1])
loh <- loh[order(loh[,1]),]
names(loh) <- c("sampleID","loh_length","n_lohsegs","tcn_length","n_tcnsegs")

resist <- resist[order(resist$sampleID),]
resist <- merge(resist,loh,by="sampleID")
for (i in 3:6)  resist[,i] <- as.numeric(as.character(resist[,i]))
table(resist$resistance)

summary(glm(I(resistance=="sensitive")~loh_length,family=binomial,data=resist)) #0.379
summary(glm(I(resistance=="sensitive")~n_lohsegs,family=binomial,data=resist)) #0.931
summary(glm(I(resistance=="sensitive")~tcn_length,family=binomial,data=resist)) #0.702
summary(glm(I(resistance=="sensitive")~n_tcnsegs,family=binomial,data=resist)) #0.871


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
save(segsall_cn,segsall_loh,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_segment.RData")


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

AUS_gene_loh=readcnasegs2(hg19=T,snp6copynumber=segsall)
AUS_gene_loh <- AUS_gene_loh[,order(names(AUS_gene_loh))]

outp <- matrix(NA,nrow(AUS_gene_loh),1)

for (i in 1:nrow(AUS_gene_loh)){
  if ((i %% 1000)==0) cat(i,"..")
  temp <- matrix(NA,ncol(AUS_gene_loh),2)
  temp[,1] <- as.vector(names(AUS_gene_loh))
  temp[,2] <- as.character(AUS_gene_loh[i,])
  
  temp2 <- data.frame(temp)
  names(temp2) <- c("sampleID","loh")
  temp2$loh <- as.numeric(as.character(temp2$loh))
  
  if (mean(temp2$loh!=0)>0.05) {
    resist1 <- merge(resist,temp2,by="sampleID")
    fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
    outp[i,1] <- summary(fit2)$coeff[2,4]
  }
}

lohpvalue=outp[,1]
names(lohpvalue)=rownames(AUS_gene_loh)

source("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
draw_manhattan(pvalues=lohpvalue,ylab="log10(p)")
save(lohpvalue,segsall,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_pvalue.Rdata")
# outp2=outp
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_pvalues_qvalues.RData")
# lohpvalue1=outp[,5]
# names(lohpvalue1)=rownames(outp)
# 
# lohpvalue=outp2[,1]
# names(lohpvalue)=rownames(TCGA_gene_loh)

lohpvalue1=lohpvalue
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_pvalue.Rdata")

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_platinum_cnloh.RData")
lohpvalue2=outp1[,5]
names(lohpvalue2)=rownames(outp1)

par(mfrow=c(2,1))
draw_manhattan(pvalues=lohpvalue,ylab="-log10(p)",main="TCGA,ASCAT")
draw_manhattan(pvalues=lohpvalue1,ylab="-log10(p)",main="AUS,ASCAT")
draw_manhattan(pvalues=lohpvalue2,ylab="-log10(p)",main="AUS,PSCBS")

#check purity and ploidy
allsamples=as.character(resist$sampleID)
info=data.frame(sampleID=resist$sampleID,purity=0,ploidy=0)
for (i in 1:length(allsamples))
{
  mysample=allsamples[i]
  dat=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.res"),header=F)
  info$purity[i]=dat$V2[2]
  info$ploidy[i]=dat$V2[3]
}
save(info,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_purity_ploidy.RData")
info_aus=info
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_purity_ploidy.RData")
t.test(info$purity,info_aus$purity)
t.test(info$ploidy,info_aus$ploidy)

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
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_pvalue.Rdata")
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
save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_cytoband_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_cytoband_0.25.RData")
#AUS_cytoband_loh=extractcytoband()
AUS_cytoband_loh=res1

#AUS_cytoband_loh=extractcytoband()
outp <- matrix(NA,nrow(AUS_cytoband_loh),1)

for (i in 1:nrow(AUS_cytoband_loh)){
  if ((i %% 100)==0) cat(i,"..")
  temp <- matrix(NA,ncol(AUS_cytoband_loh),2)
  temp[,1] <- as.vector(colnames(AUS_cytoband_loh))
  temp[,2] <- as.numeric(unlist(AUS_cytoband_loh[i,]))
  
  temp2 <- data.frame(temp)
  names(temp2) <- c("sampleID","loh")
  temp2$loh <- as.numeric(as.character(temp2$loh))
  
  if (mean(temp2$loh!=0)>0.05) {
    resist1 <- merge(resist,temp2,by="sampleID")
    fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
    outp[i,1] <- summary(fit2)$coeff[2,4]
  }
}

cytobandres_aus=cbind(cytoband,p=outp[,1])
save(cytobandres_aus,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_cytoband_0.25_pvalue.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_cytoband_0.25.RData")

cytobandres_aus=cytobandres

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25_pvalue.RData")
par(mfrow=c(2,1))
plot(-log10(cytobandres$p))
plot(-log10(cytobandres_aus$p))
sum(cytobandres$p<=0.05 & cytobandres_aus$p<=0.05,na.rm=T)

cytobandres$chrom=factor(cytobandres$chrom,levels=c(1:22,"X"))
cytobandres=cytobandres[order(cytobandres$chrom,cytobandres$start),]

plotgenome=function(chr=cytobandres$chrom,pos=cytobandres$start,value=-log10(cytobandres$p),idx=1:nrow(cytobandres),ylab="",main="",pvaluecutoff=NULL,fdrthreshold=0.05)
{
  chrs=c(1:22,"X")
  dictfile="/fh/fast/dai_j/CancerGenomics/Tools/database/reference/compact/ucsc.hg19.compact.dict"
  chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
  chrlen=chrlen$V3
  chrlen=gsub("LN:","",chrlen,fixed=T)
  chrlen=as.numeric(chrlen)
  names(chrlen)=chrs
  chrstart=rep(0,length(chrlen))
  for (i in 2:length(chrlen))
  {
    tmp=0
    for (j in 1:(i-1))
    {
      tmp=tmp+chrlen[j]
    }
    chrstart[i]=tmp
  }
  names(chrstart)=names(chrlen)
  chrend=rep(0,length(chrlen))
  for (i in 1:length(chrlen))
  {
    tmp=0
    for (j in 1:i)
    {
      tmp=tmp+chrlen[j]
    }
    chrend[i]=tmp
  }
  names(chrend)=names(chrlen)
  
  chr=chr[idx]
  pos=pos[idx]
  value=value[idx]
  res=data.frame(chr=chr,pos=pos,posall=rep(NA,length(chr)),value=value)
  res$chr=gsub(23,"X",res$chr)
  res$chr=gsub(24,"Y",res$chr)
  res$chr=factor(res$chr,levels=chrs)
  res=res[order(res$chr,res$pos),]
  res=res[res$chr %in% chrs,]
  for (mychr in chrs)
  {
    idx1=which(res$chr %in% mychr)
    idx2=which(names(chrstart)==mychr)
    res$posall[idx1]=res$pos[idx1]+chrstart[idx2]
  }
  
  ymax=1.1*max(res$value,na.rm=T)
  
  res1=res[complete.cases(res),]

  par(mar=c(4.1,4.5,2.1,2.1))
  plot(c(0,res$posall[nrow(res)]),c(-0.5,ymax),xaxt="n",yaxt="n",type="n",
       xlab="",ylab=ylab,cex=1.6,cex.lab=1.6,cex.axis=1.6,frame=F,main=main)

  x.poly <- c(res1$posall, res1$posall[length(res1$posall)], res1$posall[1])         # Adjoin two x-coordinates
  y.poly <- c(res1$value, 0, 0)                     #
  polygon(x.poly, y.poly, col=gray(0.95), border=NA)          # Show the polygon fill only
  lines(res1$posall,res1$value)
  axis(side=2, at=seq(0,ceiling(ymax),0.5),cex.axis=1.6,cex.lab=1.6,pos=0,lwd.ticks=1.1)
  #points(res1$posall,res1$value)
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  for (i in 1:length(chrs))
  {
    mychr=chrs[i]
    # rect(chrstart[i],1.05*max(res$dist2,na.rm=T),chrend[i],1.15*max(res$dist2,na.rm=T),
    #      col=color2[i],border=F)
    # text(0.5*(chrstart[i]+chrend[i]),1.15*max(res$dist2,na.rm=T),labels=chrs[i],col=color2[i],cex=1.3)
    rect(chrstart[i],-0.1,chrend[i],0,col=color2[i],border=F)
    if (i<=15 | (i>15 & i %% 2==1))
    text(0.5*(chrstart[i]+chrend[i]),-0.35,labels=chrs[i],col=color2[i],cex=1.6)
  }
  if (!is.null(pvaluecutoff))
  {
    segments(0,-log10(pvaluecutoff),par('usr')[2],-log10(pvaluecutoff),col="red")
    #text(sum(chrlen1)/2,-log10(pvaluecutoff)+0.5,paste0("FDR=",fdrthreshold))
    text(sum(chrlen)*0.9,-log10(pvaluecutoff)+0.25,paste0("FDR=",fdrthreshold),cex=1.1)
  }
  return(res)
}
par(mfrow=c(2,1))
tmp1=plotgenome(ylab="-log10(p)",main="TCGA")
tmp2=plotgenome(chr=cytobandres_aus$chrom,pos=cytobandres_aus$start,value=-log10(cytobandres_aus$p)
               ,idx=1:nrow(cytobandres_aus),ylab="-log10(p)",main="AUS")

#after running Platinum_LOH_association_permutation.R
fdr=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_qvalues_tcga_1000permutations.txt",header=T)
qcutoff=0.05
qcutoff=0.1
idx=which(fdr$qvalues_permut<=qcutoff)
pvaluecutoff=max(fdr$pvalues[idx])

par(mfrow=c(2,1))
tmp1=plotgenome(ylab="-log10(p)",main="TCGA",pvaluecutoff = pvaluecutoff,fdrthreshold = 0.1)
tmp2=plotgenome(chr=cytobandres_aus$chrom,pos=cytobandres_aus$start,value=-log10(cytobandres_aus$p)
                ,idx=1:nrow(cytobandres_aus),ylab="-log10(p)",main="AUS")

idx1=which(cytobandres_aus$p[idx]<=0.05)
cytobandres_aus[idx[idx1],]

