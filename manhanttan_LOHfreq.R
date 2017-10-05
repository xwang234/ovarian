#!/usr/bin/env Rscript
source("functions.R")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_gene_copynumber_loh.RData")
resist$resistance=as.character(resist$resistance)
resist$sampleID=as.character(resist$sampleID)
samples_sens=resist$sampleID[resist$resistance=="Sensitive"]
samples_resi=resist$sampleID[resist$resistance=="Resistant"]

sens=TCGA_gene_loh[,colnames(TCGA_gene_loh) %in% samples_sens]
resi=TCGA_gene_loh[,colnames(TCGA_gene_loh) %in% samples_resi]

draw_loh_manhanttan=function(data1=sens,main="",chrs=NULL)
{
  tmp=data1>0
  freq=rowSums(tmp)
  freq=freq/ncol(data1)
  names(freq)=rownames(data1)
  draw_manhattan(pvalues = freq,maxy=6,chrs=chrs,keepres=F,logscale=F,main=main,ylab="Frequency of LOH")
}

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_pvalues_qvalues.RData")
numit=1000
outputprefix=paste0("TCGA",5)
fdrs=read.table(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"_",numit,"permutations.txt"))

par(mfrow=c(3,1))
tcga=draw_manhattan(fdrs=fdrs,fdrthreshold=fdrcutoff,maxy=6,chrs=NULL,keepres=T,logscale=T,main=paste0("TCGA"))
draw_loh_manhanttan(data1=sens,main="Sensative",chrs=1:22)
draw_loh_manhanttan(data1=resi,main="Resistant",chrs=1:22)
par(mfrow=c(1,1))

#cytoband result----------------
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_loh_cytoband_0.25.RData")
cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband$chrom=as.character(cytoband$chrom)
cytoband=cytoband[cytoband$chrom!="chrY",]
cytoband$chrom=gsub("chr","",cytoband$chrom)
sum(rownames(res1)==paste0(cytoband$chrom,"_",cytoband$name))
resist$resistance=as.character(resist$resistance)
resist$sampleID=as.character(resist$sampleID)

samples_sens=resist$sampleID[resist$resistance=="Sensitive"]
samples_resi=resist$sampleID[resist$resistance=="Resistant"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]

#Only use TCGA data
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
resist=merge(resist,clinicaltable,by="sampleID")
colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
#training
train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
train_sens=res1[,colnames(res1) %in% train_samples_sens]
train_resi=res1[,colnames(res1) %in% train_samples_resi]
tmp=plotgenome_freq(sens=train_sens,resi=train_resi)
#testing
test_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==0]
test_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==0]
test_sens=res1[,colnames(res1) %in% test_samples_sens]
test_resi=res1[,colnames(res1) %in% test_samples_resi]
tmp=plotgenome_freq(sens=test_sens,resi=test_resi)
par(mfrow=c(2,1))
tmp1=plotgenome_freq(sens=train_sens,resi=train_resi)
tmp2=plotgenome_freq(sens=test_sens,resi=test_resi)

#AUS
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_loh_cytoband_0.25.RData")
sum(rownames(res1)==paste0(cytoband$chrom,"_",cytoband$name))
resist$resistance=as.character(resist$resistance)
resist$sampleID=as.character(resist$sampleID)
samples_sens=resist$sampleID[resist$resistance=="sensitive"]
samples_resi=resist$sampleID[resist$resistance!="sensitive"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]

t.test(unlist(sens[208,]),unlist(resi[208,]))
rownames(sens)[208]
#[1] "13_q22.3"
#TCGA:
# Welch Two Sample t-test
# 
# data:  unlist(sens[208, ]) and unlist(resi[208, ])
# t = 1.946, df = 188.3, p-value = 0.05314
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.001502417  0.221092362
# sample estimates:
#   mean of x mean of y 
# 0.6705426 0.5607477 
#AUS:
# Welch Two Sample t-test
# 
# data:  unlist(sens[208, ]) and unlist(resi[208, ])
# t = 0.81566, df = 64.416, p-value = 0.4177
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1440065  0.3427859
# sample estimates:
#   mean of x mean of y 
# 0.6129032 0.5135135 


plotgenome_freq=function(sens,resi)
{
  tmp=sens>0
  freq=rowSums(tmp)
  freq=freq/ncol(sens)
  names(freq)=rownames(sens)
  cytoband_freq=cbind(cytoband,freq_sens=freq)
  
  tmp=resi>0
  freq=rowSums(tmp)
  freq=freq/ncol(resi)
  names(freq)=rownames(resi)
  cytoband_freq=cbind(cytoband_freq,freq_resi=freq)
  
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
  
  chr=cytoband_freq$chrom
  pos=cytoband_freq$start
 
  res=data.frame(chr=chr,pos=pos,posall=rep(NA,length(chr)),value1=cytoband_freq$freq_sens,value2=cytoband_freq$freq_resi)
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
  
  ymax=1.3*max(res$value1,na.rm=T)
  
  res2=res[complete.cases(res),]
  
  par(mar=c(1.1,4.5,1.1,1.1))
  plot(c(0,res$posall[nrow(res)]),c(-0.5,ymax),xaxt="n",yaxt="n",type="n",
       xlab="",ylab="LOH Frequency",cex=1.6,cex.lab=1.6,cex.axis=1.6,frame=F,main="")
  chrs=as.character(unique(res2$chr))
  for (i in 1:length(chrs))
  {
    idx=which(res2$chr==chrs[i])
    lines(res2$posall[idx],res2$value1[idx],col="green",cex=1.2)
    lines(res2$posall[idx],res2$value2[idx],col="red",cex=1.2)
  }
  
  axis(side=2, at=seq(0,ceiling(ymax),0.5),cex.axis=1.6,cex.lab=1.6,pos=0,lwd.ticks=1.1)
  #points(res2$posall,res2$value)
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  for (i in 1:length(chrs))
  {
    mychr=chrs[i]
    # rect(chrstart[i],1.05*max(res$dist2,na.rm=T),chrend[i],1.15*max(res$dist2,na.rm=T),
    #      col=color2[i],border=F)
    # text(0.5*(chrstart[i]+chrend[i]),1.15*max(res$dist2,na.rm=T),labels=chrs[i],col=color2[i],cex=1.3)
    rect(chrstart[i],-0.05,chrend[i],0,col=color2[i],border=F)
    #if (i<=15 | (i>15 & i %% 2==1))
      text(0.5*(chrstart[i]+chrend[i]),-0.15,labels=chrs[i],col=color2[i],cex=0.8)
  }
  legend("top",legend=c("Sensitive","Resistant"),col=c("green","red"),lty=1,cex=1.2)
  return(res)
}


load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_0.25.RData")

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.25.RData")
cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband$chrom=as.character(cytoband$chrom)
cytoband=cytoband[cytoband$chrom!="chrY",]
cytoband$chrom=gsub("chr","",cytoband$chrom)
sum(rownames(res1)==paste0(cytoband$chrom,"_",cytoband$name))
resist$resistance=as.character(resist$resistance)
resist$sampleID=as.character(resist$sampleID)
samples_sens=resist$sampleID[resist$resistance=="Sensitive"]
samples_resi=resist$sampleID[resist$resistance=="Resistant"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]

#AUS
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_nocateg_0.25.RData")

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_0.25.RData")
sum(rownames(res1)==paste0(cytoband$chrom,"_",cytoband$name))
resist$resistance=as.character(resist$resistance)
resist$sampleID=as.character(resist$sampleID)
samples_sens=resist$sampleID[resist$resistance=="sensitive"]
samples_resi=resist$sampleID[resist$resistance=="resistant"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]

plotgenome_cn_freq=function(opt="gain")
{
  if (opt=="gain")
  {
    tmp=sens>0
    freq=rowSums(tmp)
    freq=freq/ncol(sens)
    names(freq)=rownames(sens)
    cytoband_freq=cbind(cytoband,freq_sens=freq)
    
    tmp=resi>0
    freq=rowSums(tmp)
    freq=freq/ncol(resi)
    names(freq)=rownames(resi)
    cytoband_freq=cbind(cytoband_freq,freq_resi=freq)
  }else
  {
    tmp=sens<0
    freq=rowSums(tmp)
    freq=freq/ncol(sens)
    names(freq)=rownames(sens)
    cytoband_freq=cbind(cytoband,freq_sens=freq)
    
    tmp=resi<0
    freq=rowSums(tmp)
    freq=freq/ncol(resi)
    names(freq)=rownames(resi)
    cytoband_freq=cbind(cytoband_freq,freq_resi=freq)
  }
  
  
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
  
  chr=cytoband_freq$chrom
  pos=cytoband_freq$start
  
  res=data.frame(chr=chr,pos=pos,posall=rep(NA,length(chr)),value1=cytoband_freq$freq_sens,value2=cytoband_freq$freq_resi)
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
  
  ymax=1.6*max(res$value1,na.rm=T)
  
  res2=res[complete.cases(res),]
  
  par(mar=c(1.1,4.5,1.1,1.1))
  plot(c(0,res$posall[nrow(res)]),c(-0.5,ymax),xaxt="n",yaxt="n",type="n",
       xlab="",ylab="TCN Frequency",cex=1.6,cex.lab=1.6,cex.axis=1.6,frame=F,main="")
  chrs=as.character(unique(res2$chr))
  for (i in 1:length(chrs))
  {
    idx=which(res2$chr==chrs[i])
    lines(res2$posall[idx],res2$value1[idx],col="green",cex=1.2)
    lines(res2$posall[idx],res2$value2[idx],col="red",cex=1.2)
  }
  
  axis(side=2, at=seq(0,ceiling(ymax),0.5),cex.axis=1.6,cex.lab=1.6,pos=0,lwd.ticks=1.1)
  #points(res2$posall,res2$value)
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  for (i in 1:length(chrs))
  {
    mychr=chrs[i]
    # rect(chrstart[i],1.05*max(res$dist2,na.rm=T),chrend[i],1.15*max(res$dist2,na.rm=T),
    #      col=color2[i],border=F)
    # text(0.5*(chrstart[i]+chrend[i]),1.15*max(res$dist2,na.rm=T),labels=chrs[i],col=color2[i],cex=1.3)
    rect(chrstart[i],-0.05,chrend[i],0,col=color2[i],border=F)
    #if (i<=15 | (i>15 & i %% 2==1))
    text(0.5*(chrstart[i]+chrend[i]),-0.15,labels=chrs[i],col=color2[i],cex=0.8)
  }
  legend("top",legend=c("Sensitive","Resistant"),col=c("green","red"),lty=1,cex=1.1)
  return(res)
}