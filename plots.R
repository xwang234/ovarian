#!/usr/bin/env Rscript
rm(list=ls())
#functions for p-values----------------------
#cmp 1 vs non-1
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
    
    #if (mean(temp2$loh!=0)>0.05) {
    if (mean(temp2$loh==1)>0.05) {  
      resist1 <- merge(resist,temp2,by="sampleID")
      resist1$resistance=as.character(resist1$resistance)
      resist1$resistance=tolower(resist1$resistance)
      #fit2 <- glm(I(resistance=="sensitive")~loh+residual_disease+grade+stage+age,family=binomial,data=resist1,x=T,y=T)
      #fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
      fit2 <- glm(I(resistance=="sensitive")~I(loh==1),family=binomial,data=resist1,x=T,y=T)
      outp[i] <- summary(fit2)$coeff[2,4]
    }
  }
  
  return(outp)
}

#include -1,0,1
compute_cytopvalue_all=function(cytodata)
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
    #temp[,2] <- as.character(cytodata[i,])
    temp[,2] <- unlist(cytodata[i,])
    
    temp2 <- data.frame(temp)
    names(temp2) <- c("sampleID","loh")
    temp2$loh <- as.numeric(as.character(temp2$loh))
    
    if (mean(temp2$loh==1)>0.05) {
      resist2 <- merge(resist1,temp2,by="sampleID")
      resist2$resistance=as.character(resist2$resistance)
      resist2$resistance=tolower(resist2$resistance)
      #fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist2,x=T,y=T)
      fit2 <- glm(I(resistance=="sensitive")~I(loh==1),family=binomial,data=resist2,x=T,y=T)
      #fit2 <- glm(I(resistance=="sensitive")~loh+residual_disease+grade+stage+age,family=binomial,data=resist2,x=T,y=T)
      res[i,1] <- summary(fit2)$coeff[2,4]
    }
  }
  return(res)
}

#include -1,0,1
computepermutation_all=function(jobn)
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

plotgenome=function(cytobandres=cytobandres,ylab="",main="",qvalues=NULL,qcutoff=0.05)
{
  chr=cytobandres$chrom
  pos=cytobandres$start
  value=-log10(cytobandres$p)
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
  
  par(mar=c(1.1,4.1,1.1,1.1))
  plot(c(0,res$posall[nrow(res)]),c(-0.3,ymax),xaxt="n",yaxt="n",type="n", mgp=c(1,1,0),
       xlab="",ylab=ylab,cex=2.2,cex.lab=2.2,cex.axis=2.2,frame=F,main=main)
  
  x.poly <- c(res1$posall, res1$posall[length(res1$posall)], res1$posall[1])         # Adjoin two x-coordinates
  y.poly <- c(res1$value, 0, 0)                     #
  polygon(x.poly, y.poly, col=gray(0.95), border=NA)          # Show the polygon fill only
  lines(res1$posall,res1$value)
  axis(side=2, at=seq(0,ceiling(ymax),0.5),cex.axis=2.2,cex.lab=2.2,pos=0,lwd.ticks=1.1)
  #points(res1$posall,res1$value)
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
    text(0.5*(chrstart[i]+chrend[i]),-0.12,labels=chrs[i],col=color2[i],cex=1.6,font=2)
    lines(c(chrend[i],chrend[i]),c(0,ymax),lty=2,col="gray",lwd=1)  
  }
  pvaluescutoff=NULL
  if (!is.null(qvalues))
  {
    idx=which(qvalues<=qcutoff)
    if (length(idx)>0) pvaluecutoff=max(10^(-value[idx]))
  }
  if (!is.null(pvaluecutoff))
  {
    segments(0,-log10(pvaluecutoff),par('usr')[2],-log10(pvaluecutoff),col="red")
    #text(sum(chrlen1)/2,-log10(pvaluecutoff)+0.5,paste0("FDR=",qcutoff))
    text(sum(chrlen)*0.75,-log10(pvaluecutoff)+0.1,paste0("FDR=",qcutoff),cex=2.2)
  }else
  {
    if (!is.null(qvalues))
    {
      text(sum(chrlen)*0.75,ymax-0.5,paste0("Minimum FDR=",round(min(qvalues,na.rm=T),2)),cex=2.2)
    }
  }
  return(res)
}
#cytoband result----------------

cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband$chrom=as.character(cytoband$chrom)
cytoband=cytoband[cytoband$chrom!="chrY",]
cytoband$chrom=gsub("chr","",cytoband$chrom)



#Only use TCGA data

#LOH freq plot
plotgenome_freq=function(sens,resi,ylab="LOH frequency")
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
  plot(c(0,res$posall[nrow(res)]),c(-0.1,1),xaxt="n",yaxt="n",type="n",mgp=c(2,1,0),
       xlab="",ylab=ylab,cex=2.2,cex.lab=2.2,cex.axis=2.2,frame=F,main="")
  chrs=as.character(unique(res2$chr))
  for (i in 1:length(chrs))
  {
    idx=which(res2$chr==chrs[i])
    lines(res2$posall[idx],res2$value1[idx],col="green",lwd=2)
    lines(res2$posall[idx],res2$value2[idx],col="red",lwd=2)
    lines(c(chrend[i],chrend[i]),c(0,1),lty=2,col="gray",lwd=1)
  }
  
  axis(side=2, at=seq(0,1,0.2),cex.axis=2.2,cex.lab=2.2,pos=0,lwd.ticks=1.1)
  #points(res2$posall,res2$value)
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  for (i in 1:length(chrs))
  {
    mychr=chrs[i]
    # rect(chrstart[i],1.05*max(res$dist2,na.rm=T),chrend[i],1.15*max(res$dist2,na.rm=T),
    #      col=color2[i],border=F)
    # text(0.5*(chrstart[i]+chrend[i]),1.15*max(res$dist2,na.rm=T),labels=chrs[i],col=color2[i],cex=1.3)
    rect(chrstart[i],-0.025,chrend[i],0,col=color2[i],border=F)
    #if (i<=15 | (i>15 & i %% 2==1))
    text(0.5*(chrstart[i]+chrend[i]),-0.05,labels=chrs[i],col=color2[i],cex=1.6,font=2)
  }
  #legend("topleft",legend=c("Sensitive","Resistant"),col=c("green","red"),lty=1,cex=1.2,bty="n")
  legend(c(0,1),legend=c("Sensitive","Resistant"),col=c("green","red"),lty=1,lwd=2,cex=2.5,bty="n")
  return(res)
}

rm(res1,resist)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_loh_cytoband_0.25.RData") #res1 is the data
sum(rownames(res1)==paste0(cytoband$chrom,"_",cytoband$name))
resist$resistance=as.character(resist$resistance)
resist$sampleID=as.character(resist$sampleID)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
clinicaltable$residual_disease <- ifelse(clinicaltable$residual_disease=="No Macroscopic disease",1,0)
clinicaltable$residual_disease[is.na(clinicaltable$residual_disease)] <- 0
clinicaltable$stage <- ifelse(clinicaltable$stage =="Stage IV",1,0)
clinicaltable$stage[is.na(clinicaltable$stage)] <- 0
clinicaltable$grade=as.character(clinicaltable$grade)
clinicaltable$grade <- ifelse(clinicaltable$grade %in% c("G3","G4"),1,0)
clinicaltable$grade[is.na(clinicaltable$grade)] <- 0
clinicaltable$grade=factor(clinicaltable$grade)

resist=merge(resist,clinicaltable,by="sampleID")
colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
#training
train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
train_sens=res1[,colnames(res1) %in% train_samples_sens]
train_resi=res1[,colnames(res1) %in% train_samples_resi]
#testing
test_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==0]
test_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==0]
test_sens=res1[,colnames(res1) %in% test_samples_sens]
test_resi=res1[,colnames(res1) %in% test_samples_resi]

#png("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_training_LOH_freq.png",res=300,width=1800,height=600,pointsize = 1)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_LOH_training_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=train_sens,resi=train_resi)
dev.off()
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_LOH_testing_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp2=plotgenome_freq(sens=test_sens,resi=test_resi)
dev.off()

#genome-wide p-value------
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)

cytodata=cbind(train_sens,train_resi)
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_loh_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_loh_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_train_ascat_loh_cytoband_nocateg_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_loh_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_loh_cytoband_qvalues0.25.RData")
cytobandres_tcga_train_ascat_loh=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_LOH_training_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_train_ascat_loh,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

pvalues_test=compute_cytopvalue(cbind(test_sens,test_resi))
cytobandres_tcga_test_ascat_loh=cbind(cytoband,p=pvalues_test)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_LOH_testing_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_test_ascat_loh,ylab="-log10(p)",main="")
dev.off()

#combined TCGA data
samples_sens=resist$sampleID[resist$resistance=="Sensitive"]
samples_resi=resist$sampleID[resist$resistance=="Resistant"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_LOH_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=sens,resi=resi)
dev.off()
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_loh_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_loh_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_ascat_loh_cytoband_nocateg_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_loh_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_loh_cytoband_qvalues0.25.RData")
cytobandres_tcga_ascat_loh=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_LOH_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_ascat_loh,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

#Australian data
rm(resist)
rm(res1)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_loh_cytoband_0.25.RData")
samples_sens=resist$sampleID[resist$resistance=="sensitive"]
samples_resi=resist$sampleID[resist$resistance!="sensitive"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_LOH_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=sens,resi=resi)
dev.off()
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_loh_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_loh_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="aus_ascat_loh_cytoband_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_loh_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_loh_cytoband_qvalues0.25.RData")
cytobandres_aus_ascat_loh=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_LOH_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_aus_ascat_loh,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

#for TCN gains------------------
rm(res1)
rm(resist)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.25.RData") 
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
clinicaltable$residual_disease <- ifelse(clinicaltable$residual_disease=="No Macroscopic disease",1,0)
clinicaltable$residual_disease[is.na(clinicaltable$residual_disease)] <- 0
clinicaltable$stage <- ifelse(clinicaltable$stage =="Stage IV",1,0)
clinicaltable$stage[is.na(clinicaltable$stage)] <- 0
clinicaltable$grade=as.character(clinicaltable$grade)
clinicaltable$grade <- ifelse(clinicaltable$grade %in% c("G3","G4"),1,0)
clinicaltable$grade[is.na(clinicaltable$grade)] <- 0
clinicaltable$grade=factor(clinicaltable$grade)

resist=merge(resist,clinicaltable,by="sampleID")
colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
#training
train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
train_sens=res1[,colnames(res1) %in% train_samples_sens]
train_resi=res1[,colnames(res1) %in% train_samples_resi]

#testing
test_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==0]
test_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==0]
test_sens=res1[,colnames(res1) %in% test_samples_sens]
test_resi=res1[,colnames(res1) %in% test_samples_resi]

#png("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_training_LOH_freq.png",res=300,width=1800,height=600,pointsize = 1)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Gain_training_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=train_sens,resi=train_resi,ylab="TCN gain frequency")
dev.off()
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Gain_testing_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp2=plotgenome_freq(sens=test_sens,resi=test_resi,ylab="TCN gain frequency")
dev.off()

cytodata=cbind(train_sens,train_resi)
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Gain_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Gain_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_train_ascat_TCN_Gain_cytoband_nocateg_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Gain_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Gain_cytoband_qvalues0.25.RData")
cytobandres_tcga_train_ascat_tcn_gain=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Gain_training_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_train_ascat_tcn_gain,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

pvalues_test=compute_cytopvalue(cbind(test_sens,test_resi))
cytobandres_tcga_test_ascat_tcn_gain=cbind(cytoband,p=pvalues_test)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Gain_testing_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_test_ascat_tcn_gain,ylab="-log10(p)",main="")
dev.off()

#combined TCGA data
samples_sens=resist$sampleID[resist$resistance=="Sensitive"]
samples_resi=resist$sampleID[resist$resistance=="Resistant"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Gain_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=sens,resi=resi,ylab="TCN gain frequency")
dev.off()
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_gain_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_gain_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_ascat_tcn_gain_cytoband_nocateg_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_gain_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_gain_cytoband_qvalues0.25.RData")
cytobandres_tcga_ascat_tcn_gain=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Gain_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_ascat_tcn_gain,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.1)
dev.off()

#Australian data
rm(resist)
rm(res1)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_0.25.RData")
samples_sens=resist$sampleID[resist$resistance=="sensitive"]
samples_resi=resist$sampleID[resist$resistance!="sensitive"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_TCN_Gain_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=sens,resi=resi,ylab="TCN gain frequency")
dev.off()
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Gain_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Gain_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="aus_ascat_TCN_Gain_cytoband_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Gain_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Gain_cytoband_qvalues0.25.RData")
cytobandres_aus_ascat_tcn_gain=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_TCN_Gain_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_aus_ascat_tcn_gain,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()


#for TCN loss--------------------------
rm(res1)
rm(resist)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.25.RData")
res1=-res1 #work on loss
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
clinicaltable$residual_disease <- ifelse(clinicaltable$residual_disease=="No Macroscopic disease",1,0)
clinicaltable$residual_disease[is.na(clinicaltable$residual_disease)] <- 0
clinicaltable$stage <- ifelse(clinicaltable$stage =="Stage IV",1,0)
clinicaltable$stage[is.na(clinicaltable$stage)] <- 0
clinicaltable$grade=as.character(clinicaltable$grade)
clinicaltable$grade <- ifelse(clinicaltable$grade %in% c("G3","G4"),1,0)
clinicaltable$grade[is.na(clinicaltable$grade)] <- 0
clinicaltable$grade=factor(clinicaltable$grade)

resist=merge(resist,clinicaltable,by="sampleID")
colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
#training
train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
train_sens=res1[,colnames(res1) %in% train_samples_sens]
train_resi=res1[,colnames(res1) %in% train_samples_resi]

#testing
test_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==0]
test_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==0]
test_sens=res1[,colnames(res1) %in% test_samples_sens]
test_resi=res1[,colnames(res1) %in% test_samples_resi]

postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Loss_training_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=train_sens,resi=train_resi,ylab="TCN loss frequency")
dev.off()
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Loss_testing_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp2=plotgenome_freq(sens=test_sens,resi=test_resi,ylab="TCN loss frequency")
dev.off()

cytodata=cbind(train_sens,train_resi)
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Loss_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Loss_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_train_ascat_TCN_Loss_cytoband_nocateg_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Loss_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_Loss_cytoband_qvalues0.25.RData")
cytobandres_tcga_train_ascat_tcn_loss=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Loss_training_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_train_ascat_tcn_loss,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.1)
dev.off()

pvalues_test=compute_cytopvalue(cbind(test_sens,test_resi))
cytobandres_tcga_test_ascat_tcn_loss=cbind(cytoband,p=pvalues_test)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Loss_testing_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_test_ascat_tcn_loss,ylab="-log10(p)",main="")
dev.off()

#combined TCGA data
samples_sens=resist$sampleID[resist$resistance=="Sensitive"]
samples_resi=resist$sampleID[resist$resistance=="Resistant"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Loss_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=sens,resi=resi,ylab="TCN loss frequency")
dev.off()
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_Loss_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_Loss_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga_ascat_tcn_Loss_cytoband_nocateg_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_Loss_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_tcn_Loss_cytoband_qvalues0.25.RData")
cytobandres_tcga_ascat_tcn_loss=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_Loss_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_ascat_tcn_loss,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.1)
dev.off()

#Australian data
rm(resist)
rm(res1)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_0.25.RData")
res1=-res1 #work on loss
samples_sens=resist$sampleID[resist$resistance=="sensitive"]
samples_resi=resist$sampleID[resist$resistance!="sensitive"]
sens=res1[,colnames(res1) %in% samples_sens]
resi=res1[,colnames(res1) %in% samples_resi]
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_TCN_Loss_freq.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp1=plotgenome_freq(sens=sens,resi=resi,ylab="TCN loss frequency")
dev.off()
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Loss_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Loss_cytoband_permutationp.RData")
pvalues=compute_cytopvalue(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="aus_ascat_TCN_Loss_cytoband_0.25",numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Loss_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_Loss_cytoband_qvalues0.25.RData")
cytobandres_aus_ascat_tcn_loss=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_TCN_Loss_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_aus_ascat_tcn_loss,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

#redraw the TCN pvalue plot using -1,0,1
#TCGA training---
rm(res1)
rm(resist)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.25.RData") 

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
clinicaltable$residual_disease <- ifelse(clinicaltable$residual_disease=="No Macroscopic disease",1,0)
clinicaltable$residual_disease[is.na(clinicaltable$residual_disease)] <- 0
clinicaltable$stage <- ifelse(clinicaltable$stage =="Stage IV",1,0)
clinicaltable$stage[is.na(clinicaltable$stage)] <- 0
clinicaltable$grade=as.character(clinicaltable$grade)
clinicaltable$grade <- ifelse(clinicaltable$grade %in% c("G3","G4"),1,0)
clinicaltable$grade[is.na(clinicaltable$grade)] <- 0
clinicaltable$grade=factor(clinicaltable$grade)

resist=merge(resist,clinicaltable,by="sampleID")
colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
#training
train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
train_sens=res1[,colnames(res1) %in% train_samples_sens]
train_resi=res1[,colnames(res1) %in% train_samples_resi]

#testing
test_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==0]
test_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==0]
test_sens=res1[,colnames(res1) %in% test_samples_sens]
test_resi=res1[,colnames(res1) %in% test_samples_resi]

cytodata=cbind(train_sens,train_resi)
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation_all,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_cytoband_permutationp.RData")
pvalues=compute_cytopvalue_all(cytodata)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix=NULL,numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_train_ascat_TCN_cytoband_qvalues0.25.RData")
cytobandres_tcga_train_ascat_tcn=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_training_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_train_ascat_tcn_gain,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

pvalues_test=compute_cytopvalue_all(cbind(test_sens,test_resi))
cytobandres_tcga_test_ascat_tcn=cbind(cytoband,p=pvalues_test)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_testing_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_test_ascat_tcn,ylab="-log10(p)",main="")
dev.off()

#combined TCGA
cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation_all,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_TCN_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_TCN_cytoband_permutationp.RData")
pvalues=compute_cytopvalue_all(cytodata)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix=NULL,numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_TCN_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat_TCN_cytoband_qvalues0.25.RData")
cytobandres_tcga_train_ascat_tcn=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_TCN_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_tcga_train_ascat_tcn,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

#Australian data
rm(resist)
rm(res1)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_ascat_tcn_cytoband_0.25.RData")

cytodata=res1
mpi.bcast.Robj2slave(cytodata)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation_all,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(cytodata))
save(pvalues_permutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_cytoband_permutationp.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_cytoband_permutationp.RData")
pvalues=compute_cytopvalue_all(cytodata)
mpi.bcast.Robj2slave(pvalues)
res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix=NULL,numit=1000)
save(res2,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_cytoband_qvalues0.25.RData")
rm(res2)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_ascat_TCN_cytoband_qvalues0.25.RData")
cytobandres_aus_ascat_tcn=cbind(cytoband,p=pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_TCN_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
tmp=plotgenome(cytobandres=cytobandres_aus_ascat_tcn,ylab="-log10(p)",main="",qvalues=res2,qcutoff=0.05)
dev.off()

#generate table for some cytobands
LOHs=c("1_p33","1_p31.3","1_p31.1","10_p15.3","10_p15.2","10_p11.1","19_p11","2_p25.3","22_q12.2","3_p26.3","3_q22.1","3_q29","4_q21.21","4_q32.2","5_q35.3","6_q25.3","7_q11.1","7_q31.32","8_q11.1","X_p21.1")
TCNs=c("19_q12","10_q23.1","6_q27","8_q24.3")
generatetable=function()
{
  res=data.frame(matrix(NA,nrow=length(LOHs)+length(TCNs),ncol=4))
  colnames(res)=c("cytoband","freq_sens","freq_resi","pvalue")
  rownames(res)=res$cytoband=c(LOHs,TCNs)
  rownames(res)[(length(LOHs)+1):nrow(res)]=paste0("TCN_",TCNs)
  #LOH
  rm(res1,resist)
  load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_loh_cytoband_0.25.RData") #res1 is the data
  resist$resistance=as.character(resist$resistance)
  resist$sampleID=as.character(resist$sampleID)
  load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
  clinicaltable$residual_disease <- ifelse(clinicaltable$residual_disease=="No Macroscopic disease",1,0)
  clinicaltable$residual_disease[is.na(clinicaltable$residual_disease)] <- 0
  clinicaltable$stage <- ifelse(clinicaltable$stage =="Stage IV",1,0)
  clinicaltable$stage[is.na(clinicaltable$stage)] <- 0
  clinicaltable$grade=as.character(clinicaltable$grade)
  clinicaltable$grade <- ifelse(clinicaltable$grade %in% c("G3","G4"),1,0)
  clinicaltable$grade[is.na(clinicaltable$grade)] <- 0
  clinicaltable$grade=factor(clinicaltable$grade)
  
  resist=merge(resist,clinicaltable,by="sampleID")
  colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
  #training
  train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
  train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
  train_sens=res1[,colnames(res1) %in% train_samples_sens]
  train_resi=res1[,colnames(res1) %in% train_samples_resi]
  cytodata=cbind(train_sens,train_resi)
  pvalues=compute_cytopvalue(cytodata)
  idx=match(rownames(res)[1:length(LOHs)],names(pvalues))
  res$pvalue[1:length(LOHs)]=round(pvalues[idx],3)
  #frequencies
  freq_sens=sapply(1:nrow(train_sens),function(i){
    sum(unlist(train_sens[i,]),na.rm=T)/ncol(train_sens)
  })
  names(freq_sens)=rownames(cytodata)
  idx=match(rownames(res)[1:length(LOHs)],names(freq_sens))
  res$freq_sens[1:length(LOHs)]=round(freq_sens[idx],3)
  
  freq_resi=sapply(1:nrow(train_resi),function(i){
    sum(unlist(train_resi[i,]),na.rm=T)/ncol(train_resi)
  })
  names(freq_resi)=rownames(cytodata)
  idx=match(rownames(res)[1:length(LOHs)],names(freq_resi))
  res$freq_resi[1:length(LOHs)]=round(freq_resi[idx],3)

  #TCNs

  rm(res1,resist)
  load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_tcn_cytoband_0.25.RData") #res1 is the data
  resist$resistance=as.character(resist$resistance)
  resist$sampleID=as.character(resist$sampleID)
  load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
  clinicaltable$residual_disease <- ifelse(clinicaltable$residual_disease=="No Macroscopic disease",1,0)
  clinicaltable$residual_disease[is.na(clinicaltable$residual_disease)] <- 0
  clinicaltable$stage <- ifelse(clinicaltable$stage =="Stage IV",1,0)
  clinicaltable$stage[is.na(clinicaltable$stage)] <- 0
  clinicaltable$grade=as.character(clinicaltable$grade)
  clinicaltable$grade <- ifelse(clinicaltable$grade %in% c("G3","G4"),1,0)
  clinicaltable$grade[is.na(clinicaltable$grade)] <- 0
  clinicaltable$grade=factor(clinicaltable$grade)
  
  resist=merge(resist,clinicaltable,by="sampleID")
  colnames(resist)[which(colnames(resist)=="resistance.x")]="resistance"
  #training
  train_samples_sens=resist$sampleID[resist$resistance=="Sensitive" & resist$before2009==1]
  train_samples_resi=resist$sampleID[resist$resistance=="Resistant" & resist$before2009==1]
  train_sens=res1[,colnames(res1) %in% train_samples_sens]
  train_resi=res1[,colnames(res1) %in% train_samples_resi]
  cytodata=cbind(train_sens,train_resi)
  pvalues=compute_cytopvalue_all(cytodata)
  idx=match(TCNs,names(pvalues))
  res$pvalue[(length(LOHs)+1):nrow(res)]=round(pvalues[idx],3)
  #frequencies
  freq_sens=sapply(1:nrow(train_sens),function(i){
    sum(unlist(train_sens[i,]),na.rm=T)/ncol(train_sens)
  })
  names(freq_sens)=rownames(cytodata)
  idx=match(TCNs,names(freq_sens))
  res$freq_sens[(length(LOHs)+1):nrow(res)]=round(freq_sens[idx],3)
  
  freq_resi=sapply(1:nrow(train_resi),function(i){
    sum(unlist(train_resi[i,]),na.rm=T)/ncol(train_resi)
  })
  names(freq_resi)=rownames(cytodata)
  idx=match(TCNs,names(freq_resi))
  res$freq_resi[(length(LOHs)+1):nrow(res)]=round(freq_resi[idx],3)
  write.table(res,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/lassoselectedcytoband.txt",row.names = F,col.names = T,sep="\t",quote=F)
}

