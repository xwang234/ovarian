#!/usr/bin/env Rscript

library(GenomicRanges)
cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband$chrom=as.character(cytoband$chrom)
cytoband=cytoband[cytoband$chrom!="chrY",]
cytoband$chrom=gsub("chr","",cytoband$chrom)
gr_cytoband=GRanges(seqnames = cytoband$chrom,ranges=IRanges(start=cytoband$start,end=cytoband$end))

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
                      main1="",main2="",qvalues = NULL,qcutoff =NULL)
{
  par(mfrow=c(2,1))
  if (!is.null(qvalues) & is.null(qcutoff))
  {
    idx=which(qvalues<=qcutoff)
    pvaluecutoff=max(pvalues[idx])
    tmp1=plotgenome(chr=cytobandres1$chrom,pos=cytobandres1$start,value=-log10(cytobandres1$p)
                    ,idx=1:nrow(cytobandres1),ylab="-log10(p)",main=main1,pvaluecutoff = pvaluecutoff,qcutoff = qcutoff)
    idx=which(tmp1$value>=-log10(pvaluecutoff))
    tmp2=plotgenome(chr=cytobandres_aus_pscbs_loh$chrom,pos=cytobandres_aus_pscbs_loh$start,value=-log10(cytobandres_aus_pscbs_loh$p)
                    ,idx=1:nrow(cytobandres_aus_pscbs_loh),ylab="-log10(p)",main="AUS_PSCBS_LOH")
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
#TCGA---LOH-----
#load segmentation
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_segment.RData") #Platinum_LOH_association_w.R
# resist_pscbs=resist
# resist_pscbs$sampleID=as.character(resist_pscbs$sampleID)
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25.RData")
# resist_ascat=resist
# resist_ascat$sampleID=as.character(resist_ascat$sampleID)
# comsamples=resist_ascat$sampleID[resist_ascat$sampleID %in% resist_pscbs$sampleID]
# resist2=merge(resist_ascat,resist_pscbs,by="sampleID")
# cor(resist2$loh_length.x,resist2$loh_length.y)
# #[1] 0.8058196
# cor(resist2$n_lohsegs.x,resist2$n_lohsegs.y)
# #[1] 0.5894469
# cor(resist2$tcn_length.x,resist2$tcn_length.y)
# #[1] -0.0155721

seg=segsall_loh[,-5]
#categorize
seg$dhMean=1*(seg$dhMean>=0.2)
#salloc -t 0-5 -n 49 mpirun -n 1 R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)

mpi.bcast.Robj2slave(seg)
mpi.bcast.Robj2slave(cytoband)
mpi.bcast.Robj2slave(gr_cytoband)
mpi.bcast.cmd(library(GenomicRanges))

allsamples=unique(as.character(seg$ID))
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband_mpi,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_0.25.RData")
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
cytobandres_tcga_pscbs_loh=cbind(cytoband,p=pvalues)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_pscbs_loh",numit=1000)
#save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_qvalues.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_qvalues.RData")
qcutoff=0.05
qcutoff=0.1
idx=which(res2<=qcutoff)
pvaluecutoff=max(pvalues[idx])

tmp1=plotgenome(chr=cytobandres_tcga_pscbs_loh$chrom,pos=cytobandres_tcga_pscbs_loh$start,value=-log10(cytobandres_tcga_pscbs_loh$p)
                ,idx=1:nrow(cytobandres_tcga_pscbs_loh),ylab="-log10(p)",main="TCGA_PSCBS_LOH",pvaluecutoff = pvaluecutoff)

#TCGA---Copynumber
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_segment.RData") #Platinum_LOH_association_w.R
seg=segsall_cn[,-5]
idx1=seg$tcnMean[idx]>=0.5
idx2=seg$tcnMean[idx]<= -0.4
seg$tcnMean[idx[idx1]]=1
seg$tcnMean[idx[idx2]]=-1
seg$tcnMean[idx[(!idx1) & (!idx2)]]=0
mpi.bcast.Robj2slave(seg)

allsamples=unique(as.character(seg$ID))
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband1_mpi,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_0.25.RData")
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_pscbs_loh_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
cytobandres_tcga_pscbs_loh=cbind(cytoband,p=pvalues)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_pscbs_loh",numit=1000)




#AUS---LOH----
allsamples=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/allsamples_aus_loh_pscbs.txt")
segsall=c()
for (mysample in allsamples[,1])
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.PSCBS.segment.txt"),header=T)
  names(tmp)[1] <- "ID"
  tmp <- tmp[,c(1,2,5,6,7,8,11,12,10,14,34)]
  segsall=rbind(segsall,tmp)
}
segsall=segsall[complete.cases(segsall),]
segsall_cn <- segsall[!is.na(segsall$tcnMean),c(1:6)]
segsall_cn$tcnMean <- segsall_cn$tcnMean -2 
segsall_loh <- segsall[segsall$lohCall & !is.na(segsall$dhMean),c(1:2,7:10)]
save(segsall_cn,segsall_loh,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_pscbs_segment.RData")
seg=segsall_loh[,-5]
#categorize
seg$dhMean=1*(seg$dhMean>=0.2)
mpi.bcast.Robj2slave(seg)
allsamples=unique(as.character(seg$ID))
res=mpi.parSapply(X=1:length(allsamples),FUN=extractcytoband_mpi,job.num=njobs)
res1=matrix(unlist(res),byrow = F,ncol = length(allsamples))
res1=as.data.frame(res1)
rownames(res1)=paste0(cytoband$chrom,"_",cytoband$name)
colnames(res1)=allsamples
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_pscbs_loh_resist.RData")
#save(res1,resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_pscbs_loh_cytoband_0.25.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_pscbs_loh_cytoband_0.25.RData")
cytodata=res1
pvalues_aus=compute_cytopvalue(cytodata)
cytobandres_aus_pscbs_loh=cbind(cytoband,p=pvalues_aus)

par(mfrow=c(2,1))
tmp1=plotgenome(chr=cytobandres_tcga_pscbs_loh$chrom,pos=cytobandres_tcga_pscbs_loh$start,value=-log10(cytobandres_tcga_pscbs_loh$p)
                ,idx=1:nrow(cytobandres_tcga_pscbs_loh),ylab="-log10(p)",main="TCGA_PSCBS_LOH",pvaluecutoff = pvaluecutoff,qcutoff = qcutoff)
tmp2=plotgenome(chr=cytobandres_aus_pscbs_loh$chrom,pos=cytobandres_aus_pscbs_loh$start,value=-log10(cytobandres_aus_pscbs_loh$p)
                ,idx=1:nrow(cytobandres_aus_pscbs_loh),ylab="-log10(p)",main="AUS_PSCBS_LOH")

plot2genomes(cytobandres1=cytobandres_tcga_pscbs_loh,cytobandres2=cytobandres_aus_pscbs_loh,
                      main1="TCGA_PSCBS_LOH",main2="AUS_PSCBS_LOH",qvalues=res2,qcutoff =qcutoff)

idx=which(res2<=qcutoff)
idx1=which(cytobandres_aus_pscbs_loh$p[idx]<=0.05)
cytobandres_aus_pscbs_loh[idx[idx1],]




