#!/usr/bin/env Rscript

library(GenomicRanges)
cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband$chrom=as.character(cytoband$chrom)
cytoband=cytoband[cytoband$chrom!="chrY",]
cytoband$chrom=gsub("chr","",cytoband$chrom)
gr_cytoband=GRanges(seqnames = cytoband$chrom,ranges=IRanges(start=cytoband$start,end=cytoband$end))

#loh: 0,1
extractcytoband_mpi=function(i,colID=1,colchr=2,colstart=3,colend=4,colloh=5,cutoff=0.25)
{
  allsamples=unique(as.character(seg[,colID]))
  res=data.frame(matrix(0,ncol=nrow(cytoband),nrow=1))
  colnames(res)=paste0(cytoband$chrom,"_",cytoband$name)
  rownames(res)=allsamples[i]
  sampleseg=seg[seg[,colID]==allsamples[i] & seg[,colloh]==1,]
  
  gr_sampleseg=GRanges(seqnames = sampleseg[,colchr],ranges=IRanges(start=sampleseg[,colstart],end=sampleseg[,colend]),loh=sampleseg[,colloh])
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

#copynumber -1,0,1
extractcytoband1_mpi=function(i,colID=1,colchr=2,colstart=3,colend=4,colloh=5,cutoff=0.25)
{
  allsamples=unique(as.character(seg[,colID]))
  res=data.frame(matrix(0,ncol=nrow(cytoband),nrow=1))
  colnames(res)=paste0(cytoband$chrom,"_",cytoband$name)
  rownames(res)=allsamples[i]
  sampleseg=seg[seg[,colID]==allsamples[i] & seg[,colloh]!=0,]
  
  gr_sampleseg=GRanges(seqnames = sampleseg[,colchr],ranges=IRanges(start=sampleseg[,colstart],end=sampleseg[,colend]),loh=sampleseg[,colloh])
  idx=which(mcols(gr_sampleseg)$loh==1)
  gr_sampleseg_gain=gr_sampleseg[idx]
  idx=which(mcols(gr_sampleseg)$loh==-1)
  gr_sampleseg_loss=gr_sampleseg[idx]
  for (j in 1:length(gr_cytoband))
  {
    olap1=intersect(gr_cytoband[j],gr_sampleseg_gain)
    olap2=intersect(gr_cytoband[j],gr_sampleseg_loss)
    
    
    if ( (length(olap1)>0 & length(olap2)==0) | (length(olap1)>0 & length(olap2)>0 & sum(width(olap1))> sum(width(olap2))) )
    {
      if (sum(width(olap1))>=cutoff*width(gr_cytoband[j]))
      {
        res[1,j]=1
      }
    }
    if ( (length(olap2)>0 & length(olap1)==0) | (length(olap1)>0 & length(olap2)>0 & sum(width(olap2))>= sum(width(olap1))) )
    {
      if (sum(width(olap2))>=cutoff*width(gr_cytoband[j]))
      {
        res[1,j]=-1
      }
    }
  }
  return(res)
}

#copynumber -1,0,1, cutoff and no categorization
extractcytoband2_mpi=function(i,colID=1,colchr=2,colstart=3,colend=4,colloh=5,cutoff=0.25)
{
  allsamples=unique(as.character(seg[,colID]))
  res=data.frame(matrix(0,ncol=nrow(cytoband),nrow=1))
  colnames(res)=paste0(cytoband$chrom,"_",cytoband$name)
  rownames(res)=allsamples[i]
  sampleseg=seg[seg[,colID]==allsamples[i] & seg[,colloh]!=0,]
  
  gr_sampleseg=GRanges(seqnames = sampleseg[,colchr],ranges=IRanges(start=sampleseg[,colstart],end=sampleseg[,colend]),loh=sampleseg[,colloh])
  idx=which(mcols(gr_sampleseg)$loh>0)
  gr_sampleseg_gain=gr_sampleseg[idx]
  idx=which(mcols(gr_sampleseg)$loh<0)
  gr_sampleseg_loss=gr_sampleseg[idx]
  for (j in 1:length(gr_cytoband))
  {
    olap1=intersect(gr_cytoband[j],gr_sampleseg_gain)
    olap2=intersect(gr_cytoband[j],gr_sampleseg_loss)
    
    
    if ( (length(olap1)>0 & length(olap2)==0) | (length(olap1)>0 & length(olap2)>0 & sum(width(olap1))> sum(width(olap2))) )
    {
      if (sum(width(olap1))>=cutoff*width(gr_cytoband[j]))
      {
        olap11=subsetByOverlaps(gr_sampleseg_gain,olap1)
        if (length(olap11)==1)
        {
          res[1,j]=mcols(olap11)$loh
        }else
        {
          if (length(olap11)>1)
          {
            totlen=sum(width(olap1))
            tmp=0
            for (k in 1:length(olap1))
            {
              tmp=tmp+width(olap1)[k]*mcols(olap11)$loh[k]
            }
            res[1,j]=tmp/totlen
          }
        }
      }
    }
    if ( (length(olap2)>0 & length(olap1)==0) | (length(olap1)>0 & length(olap2)>0 & sum(width(olap2))>= sum(width(olap1))) )
    {
      if (sum(width(olap2))>=cutoff*width(gr_cytoband[j]))
      {
        olap22=subsetByOverlaps(gr_sampleseg_loss,olap2)
        if (length(olap22)==1)
        {
          res[1,j]=mcols(olap22)$loh
        }else
        {
          if (length(olap22)>1)
          {
            totlen=sum(width(olap2))
            tmp=0
            for (k in 1:length(olap2))
            {
              tmp=tmp+width(olap2)[k]*mcols(olap22)$loh[k]
            }
            res[1,j]=tmp/totlen
          }
        }
      }
    }
  }
  return(res)
}


compute_cytopvalue=function(cytodata)
{
  outp <- rep(NA,nrow(cytodata))
  names(outp)=rownames(cytodata)
  for (i in 1:nrow(cytodata)){
    if ((i %% 100)==0) cat(i,"..")
    temp <- matrix(NA,ncol(cytodata),2)
    temp[,1] <- as.vector(colnames(cytodata))
    temp[,2] <- as.numeric(unlist(cytodata[i,]))
    
    temp2 <- data.frame(temp)
    names(temp2) <- c("sampleID","loh")
    temp2$loh <- as.numeric(as.character(temp2$loh))
    
    if (mean(temp2$loh!=0)>0.05) {
      resist1 <- merge(resist,temp2,by="sampleID")
      resist1$resistance=as.character(resist1$resistance)
      resist1$resistance=tolower(resist1$resistance)
      fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
      outp[i] <- summary(fit2)$coeff[2,4]
    }
  }
  
  return(outp)
}

computepermutation=function(jobn)
{
  set.seed(jobn+1000)
  #shuffle the outcome
  resist1=resist
  resist1$resistance=resist$resistance[sample(length(resist$resistance))]
  res=data.frame(matrix(NA,nrow=nrow(cytodata),ncol=1))
  colnames(res)="p"
  rownames(res)=rownames(cytodata)
  
  for (i in 1:nrow(cytodata)){
    temp <- matrix(NA,ncol(cytodata),2)
    temp[,1] <- as.vector(colnames(cytodata))
    temp[,2] <- as.character(cytodata[i,])
    
    temp2 <- data.frame(temp)
    names(temp2) <- c("sampleID","loh")
    temp2$loh <- as.numeric(as.character(temp2$loh))
    
    if (mean(temp2$loh!=0)>0.05) {
      resist2 <- merge(resist1,temp2,by="sampleID")
      resist2$resistance=as.character(resist2$resistance)
      resist2$resistance=tolower(resist2$resistance)
      fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist2,x=T,y=T)
      res[i,1] <- summary(fit2)$coeff[2,4]
    }
  }
  return(res)
}

compute_aqvalue_permutation=function(pvalue,pvalues,pvalues_permutation,numit=1000)
{
  #numit=nrow(pvalues_permutation)/nrow(pvalues)
  qvalue=sum(pvalues_permutation<=pvalue,na.rm=T)/numit/sum(pvalues<=pvalue,na.rm=T)
}

mpi_compute_qvalue_permutation=function(pvalues,pvalues_permutation,numit,outputprefix=NULL)
{
  #pvalues is a named vector
  mpi.bcast.Robj2slave(pvalues_permutation)
  #keep the orignial order of pvalues
  namespvalues=names(pvalues)
  #the decreasingly ordered pvalues
  pvalues1=pvalues[order(pvalues,decreasing =T)]
  mpi.bcast.Robj2slave(pvalues1)
  mpi.bcast.Robj2slave(numit)
  #work on nonNA pvalues
  idxnoNA=which(!is.na(pvalues1))
  qvalues1=rep(NA,length(pvalues1))
  names(qvalues1)=names(pvalues1)
  res1 <- NULL
  nrun <- ceiling(length(idxnoNA)/1000)
  print(paste0("total num of run: ",nrun))
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):length(idxnoNA)
    z=pvalues1[idxnoNA[cseq]]
    res=mpi.parSapply(X=z,FUN=compute_aqvalue_permutation,pvalues=pvalues1,pvalues_permutation=pvalues_permutation,numit=numit,job.num=njobs)
    res1=c(res1,res)
  }
  qvalues1[idxnoNA]=res1
  #make monotone
  qvalues=rep(NA,length(pvalues1))
  for (i in idxnoNA)
  {
    
    #make q values monotone
    if (i==idxnoNA[1])
    {
      #keep the first q value associated with the highest p-value
      qvalues[i]=qvalues1[i]
    }else
    {
      #make sure q-value is non-increasing
      qvalues[i]=min(qvalues[i-1],qvalues1[i],na.rm=T)
    }
    if (i %% 10000==0) cat(i,"..")
  }
  
  #use the original order of p-values
  idx=match(namespvalues,names(qvalues1))
  qvalues=qvalues[idx]
  names(qvalues)=namespvalues
  
  res2=data.frame(pvalues=pvalues,qvalues_permut=qvalues)
  rownames(res2)=namespvalues
  #save result
  if (!is.null(outputprefix))
  {
    filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"_",numit,"permutations.txt")
    write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
  }
  print("done")
  return(qvalues)
}

plotgenome=function(chr=cytobandres$chrom,pos=cytobandres$start,value=-log10(cytobandres$p),idx=1:nrow(cytobandres),ylab="",main="",pvaluecutoff=NULL,qcutoff=0.05)
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
    #text(sum(chrlen1)/2,-log10(pvaluecutoff)+0.5,paste0("FDR=",qcutoff))
    text(sum(chrlen)*0.9,-log10(pvaluecutoff)+0.25,paste0("FDR=",qcutoff),cex=1.1)
  }
  return(res)
}

plot2genomes=function(cytobandres1=cytobandres_tcga_pscbs_loh,cytobandres2=cytobandres_aus_pscbs_loh,
                      main1="",main2="",qvalues = NULL,qcutoff =NULL,pvalues=NULL)
{
  par(mfrow=c(2,1))
  if (!is.null(qvalues) & !is.null(qcutoff))
  {
    idx=which(qvalues<=qcutoff)
    pvaluecutoff=max(pvalues[idx])
    tmp1=plotgenome(chr=cytobandres1$chrom,pos=cytobandres1$start,value=-log10(cytobandres1$p)
                    ,idx=1:nrow(cytobandres1),ylab="-log10(p)",main=main1,pvaluecutoff = pvaluecutoff,qcutoff = qcutoff)
    idx=which(tmp1$value>=-log10(pvaluecutoff))
    tmp2=plotgenome(chr=cytobandres2$chrom,pos=cytobandres2$start,value=-log10(cytobandres2$p)
                    ,idx=1:nrow(cytobandres2),ylab="-log10(p)",main=main2)
    idx1=which(tmp2$value[idx]>=-log10(0.05))
    if (length(idx1)>0)
    {
      print(paste0(length(idx1)," cytobands were validated"))
      print(tmp2[idx[idx1],])
      for (i in 1:length(idx1))
      {
        arrows(tmp2$posall[idx[idx1[i]]], tmp2$value[idx[idx1[i]]]+0.3,tmp2$posall[idx[idx1[i]]], tmp2$value[idx[idx1[i]]]+0.1,length=0.05, col="red")
      }
    }
  }
}
#TCGA---COPYNUMBER-------
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_purity_ploidy.RData")
#load segmentation
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_segment.RData")
seg=segsall_cn
allsamples=unique(as.character(seg$ID))
#categorize
if (opt==1)
{
  idx1=which(seg$tcn>0)
  idx2=which(seg$tcn<0)
  seg$tcn[idx1]=1
  seg$tcn[idx2]=-1
}

if (opt==2)
{
  for (i in 1:length(allsamples))
  {
    idx=which(seg$ID==allsamples[i])
    ploidy=info$ploidy[info$sampleID==allsamples[i]]
    idx1=seg$tcn[idx]>=ploidy-2+1
    idx2=seg$tcn[idx]<=ploidy-2-1
    seg$tcn[idx[idx1]]=1
    seg$tcn[idx[idx2]]=-1
    seg$tcn[idx[(!idx1) & (!idx2)]]=0
  }
}

if (opt==3)
{
  #no cagegorization
}

#salloc -t 1-5 -n 49 mpirun -n 1 R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)

mpi.bcast.Robj2slave(seg)
mpi.bcast.Robj2slave(cytoband)
mpi.bcast.Robj2slave(gr_cytoband)
mpi.bcast.Robj2slave(allsamples)
mpi.bcast.cmd(library(GenomicRanges))

res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband1_mpi,cutoff=0.9,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.9.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.5.RData")
#opt=2
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_ploidy0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_ploidy0.25.RData")

cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))

#save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_permutationp.RData")

#opt=2
#save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_ploidy_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_ploidy_permutationp.RData")

pvalues=compute_cytopvalue(cytodata)
pvalues_tcga=pvalues
cytobandres_tcga_ascat_tcn=cbind(cytoband,p=pvalues)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_ascat_tcn",numit=1000)
#save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_qvalues0.35.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_qvalues0.35.RData")

#opt=2
#save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_ploidy_qvalues.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_ploidy_qvalues.RData")

qcutoff=0.05
qcutoff=0.2
idx=which(res2<=qcutoff)
pvaluecutoff=max(pvalues[idx])
tmp1=plotgenome(chr=cytobandres_tcga_ascat_tcn$chrom,pos=cytobandres_tcga_ascat_tcn$start,value=-log10(cytobandres_tcga_ascat_tcn$p)
                ,idx=1:nrow(cytobandres_tcga_ascat_tcn),ylab="-log10(p)",main="TCGA_ASCAT_TCN",pvaluecutoff = pvaluecutoff,qcutoff = qcutoff)


#opt=3
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband2_mpi,cutoff=0.25,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_0.25.RData")
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
#save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_0.25.permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_0.25.permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
pvalues_tcga=pvalues
cytobandres_tcga_ascat_tcn=cbind(cytoband,p=pvalues)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_ascat_tcn_cytoband_nocateg_0.25",numit=1000)
#save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_qvalues0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_nocateg_qvalues0.25.RData")



#AUS---COPYNUMBER-----
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_purity_ploidy.RData")
#load segmentation
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_segment.RData")
seg=segsall_cn
allsamples=unique(as.character(seg$ID))
if (opt==1)
{
  idx1=which(seg$tcn>0)
  idx2=which(seg$tcn<0)
  seg$tcn[idx1]=1
  seg$tcn[idx2]=-1
}


if (opt==2)
{
  for (i in 1:length(allsamples))
  {
    idx=which(seg$ID==allsamples[i])
    ploidy=info$ploidy[info$sampleID==allsamples[i]]
    idx1=seg$tcn[idx]>=ploidy-2+1
    idx2=seg$tcn[idx]<=ploidy-2-1
    seg$tcn[idx[idx1]]=1
    seg$tcn[idx[idx2]]=-1
    seg$tcn[idx[(!idx1) & (!idx2)]]=0
  }
}


mpi.bcast.Robj2slave(seg)
mpi.bcast.Robj2slave(allsamples)
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband1_mpi,cutoff=0.85,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples

#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_0.85.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_0.85.RData")
#opt=2
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_ploidy0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_ploidy0.25.RData")

#opt=3
if (opt==3)
{
  #no cagegorization
}
mpi.bcast.Robj2slave(seg)
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband2_mpi,cutoff=0.25,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_nocateg_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_nocateg_0.25.RData")



cytodata=res1

pvalues=compute_cytopvalue(cytodata)
cytobandres_aus_ascat_tcn=cbind(cytoband,p=pvalues)

qcutoff=0.3
plot2genomes(cytobandres1=cytobandres_tcga_ascat_tcn,cytobandres2=cytobandres_aus_ascat_tcn,
             main1="TCGA_ASCAT_TCN",main2="AUS_ASCAT_TCN",qvalues=res2,qcutoff =qcutoff,pvalues = pvalues_tcga)

idx=which(res2<=qcutoff)
idx1=which(cytobandres_aus_ascat_tcn$p[idx]<=0.05)
cytobandres_aus_ascat_tcn[idx[idx1],]






#TCGA---LOH-----
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25_pvalue.RData") #Platinum_LOH_ASCAT_assocation.R



#AUS---LOH----
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_cytoband_0.25_pvalue.RData")  #Platinum_AUS_LOH_ASCAT_assocation.R


res2=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_qvalues_tcga_1000permutations.txt")
pvalues_tcga=res2[,1]
names(pvalues_tcga)=rownames(res2)
res2=res2[,2]
qcutoff=0.05
qcutoff=0.1

plot2genomes(cytobandres1=cytobandres,cytobandres2=cytobandres_aus,
                      main1="TCGA_ASCAT_LOH",main2="AUS_ASCAT_LOH",qvalues=res2,qcutoff =qcutoff,pvalues = pvalues_tcga)



idx=which(res2<=qcutoff)
idx1=which(cytobandres_aus_pscbs_loh$p[idx]<=0.05)
cytobandres_aus_pscbs_loh[idx[idx1],]


#compare segment results on TCGA tumor-nomal pairs with tumor-only
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_segment.RData")
segsall_loh_pair=segsall_loh
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tumor_segment.RData")
samples_pair=unique(as.character(segsall_loh_pair$ID))
samples_tumor=unique(as.character(segsall_loh$ID))
allsamples=intersect(samples_pair,samples_tumor)
res=data.frame(matrix(NA,nrow=length(allsamples),ncol=5))
colnames(res)=c("pair_len","pair_num","tumor_len","tumor_num","overlap_len")
rownames(res)=allsamples
for (i in 1:length(allsamples))
{
  mysample=allsamples[i]
  tmp1=segsall_loh_pair[segsall_loh_pair$ID==mysample,]
  res$pair_len[i]=sum(tmp1$endpos-tmp1$startpos)/10^6
  res$pair_num[i]=nrow(tmp1)
  gr_tmp1=GRanges(seqnames = tmp1$chr,ranges = IRanges(start=tmp1$startpos,end=tmp1$endpos))
  tmp2=segsall_loh[segsall_loh$ID==mysample,]
  res$tumor_len[i]=sum(tmp2$endpos-tmp2$startpos)/10^6
  res$tumor_num[i]=nrow(tmp2)
  gr_tmp2=GRanges(seqnames = tmp2$chr,ranges = IRanges(start=tmp2$startpos,end=tmp2$endpos))
  res$overlap_len[i]=sum(width(intersect(gr_tmp1,gr_tmp2)))/10^6
}
res$ratio_olap_pair=round(res$overlap_len/res$pair_len*100,digits=2)
res$ratio_olap_tumor=round(res$overlap_len/res$tumor_len*100,digits = 2)
sum(res$ratio_olap_pair<50 & res$ratio_olap_tumor<50)
#[1] 20
sum(res$ratio_olap_pair>50 & res$ratio_olap_tumor>50)
#[1] 289
plot(res$pair_len,res$tumor_len,xlab="Paried-samples (Mb)",ylab="Tumor-only (Mb)",main="Total length of segments")
abline(1,1,col="red")
cor(res$pair_len,res$tumor_len)
t.test(res$pair_len,res$tumor_len,paired = T)

plot(res$pair_num,res$tumor_num,xlab="Paried-samples",ylab="Tumor-only",main="Number of segments")
abline(1,1,col="red")
cor(res$pair_num,res$tumor_num)
t.test(res$pair_num,res$tumor_num,paired = T)

#one similar case
#one discrepant case
which.max(res$tumor_num)
samplename="TCGA-23-1110"
rownames(res)[132]
samplename="TCGA-23-1110"
seg_pair=segsall_loh_pair[segsall_loh_pair$ID==samplename,]
seg_tumor=segsall_loh[segsall_loh$ID==samplename,]
gr_seg_pair=GRanges(seqnames = seg_pair$chr,ranges = IRanges(start=seg_pair$startpos,end=seg_pair$endpos))
gr_seg_tumor=GRanges(seqnames = seg_tumor$chr,ranges = IRanges(start=seg_tumor$startpos,end=seg_tumor$endpos))
olap=intersect(gr_seg_pair,gr_seg_tumor)
length(gr_seg_tumor)
length(gr_seg_pair)
length(olap)
sum(width(gr_seg_tumor))
sum(width(gr_seg_pair))
sum(width(olap))





