#!/usr/bin/env Rscript
ggd.qqplot = function(pvector, main=NULL, xlab=NULL, ylab=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  idx=is.finite(o)
  o=o[idx]
  e = -log10(1:length(o)/length(o))
  if (is.null(xlab))
  {
    xlab="Expected"
  }else
  {
    xlab=paste0(xlab," expected")
  }
  if (is.null(ylab))
  {
    ylab="Observed"
  }else
  {
    ylab=paste0(ylab," observed")
  }
  
  plot(e,o,pch=19,cex=1, main=main, ...,
       #xlab=expression(Expected~~-log[10](italic(p))),
       #ylab=expression(Observed~~-log[10](italic(p))),
       xlab=xlab,
       ylab=ylab,
       xlim=c(0,max(e)), ylim=c(0,max(o)+1),cex.axis=1.5,cex.lab=1.5)
  #lines(e,e,col="gray55")
  abline(a=0,b=1,col="gray33",lwd=3,lty=3)
}

library(gap)
#tmp=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/mrna_geneposition_knowntable.txt",header=T,sep="\t",stringsAsFactors=F)
tmp=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_geneposition_knowntable.txt",header=T,sep="\t",stringsAsFactors=F)
genes_chr=tmp[,"chr"]
genes_pos=tmp[,"start"]
idxkeep=which(genes_chr %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
pvalues_copynumber_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
pvalues_mrna_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum.txt",header=T)
#genedata=data.frame(chromosome=genes_chr[idxkeep],position=genes_pos[idxkeep],pvalues=pvalues_mrna_platinum[idxkeep,1],stringsAsFactors=F)
genedata=data.frame(chromosome=genes_chr[idxkeep],position=genes_pos[idxkeep],pvalues=pvalues_copynumber_platinum[idxkeep,1],stringsAsFactors=F)
sortgenedata=function(genedata)
{
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genedata)))
  for (chr in chrs)
  {
    tmptable=genedata[which(genedata$chromosome==chr),]
    tmptable=tmptable[order(tmptable$position),]
    res=rbind(res,tmptable)
  }
  return(res)
}
genedata=sortgenedata(genedata)
plotpos=function(genedata)
{
  par(mfrow=c(2,1))
  ggd.qqplot(genedata$pvalues, "Q-Q plot" )
  
  #color2<-rep(c("brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
  #              "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4","azure4"),2)
  
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  par(las=2, xpd=TRUE,cex.axis=1.25,cex=1,cex.lab=2.5)
  ops<-mht.control(colors=color2,yline=2,xline=3,usepos=FALSE,srt=0)
  mhtplot(genedata,ops,pch=10,bg=color2,ylab=expression(Observed~~-log[10](italic(p))),
          cex.axis=1,cex.lab=1)
  axis(2,at=c(0:10))
  axis(1,pos=0,labels=FALSE,tick=FALSE)
  abline(0,0)
  #abline(h=7,col="black")
}
plotpos(genedata)
par(mfrow=c(2,1))
hist(pvalues_mrna_platinum[,1],main="mrna")
hist(pvalues_copynumber_platinum[,1],main="copynumber")
qvalues_mrna=p.adjust(pvalues_mrna_platinum[,1],method = "BH")
qvalues_copynumber=p.adjust(pvalues_copynumber_platinum[,1],method = "BH")