
##salloc -n 41 -t 24:0:0 mpirun -n 1 Rscript --no-save --no-restore PSCBS_LOH_try.R

rm(list=ls())

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



snp_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno <- snp_anno[,1:4]


library("DNAcopy")
library("PSCBS")



nsample <- length(allsamples)
  
  for (i in 1:nsample){
    
  cat("sample=",i,"\n")  
  mysample=allsamples[i]
  tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[i]],"/",filelist$filename[idxtumor[i]]),
                        skip = 1, sep="\t",header=T,stringsAsFactors = F)
  tumortable$Signal_A[which(!is.na(tumortable$Signal_A) & tumortable$Signal_A<0)] <- NA
  tumortable$Signal_B[which(!is.na(tumortable$Signal_A) & tumortable$Signal_B<0)] <- NA
  
  colnames(tumortable)[1]="Probe.Set.ID"
  
  
  #if it has normal allele files, find heterogeneity probes first
  if (!is.na(idxnormal[i]))
  {
    normaltable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxnormal[i]],"/",filelist$filename[idxnormal[i]]),
                           skip = 1, sep="\t",header=T,stringsAsFactors = F)
    
    normaltable$Signal_A[which(!is.na(normaltable$Signal_A) & normaltable$Signal_A<0)] <- NA
    normaltable$Signal_B[which(!is.na(normaltable$Signal_A) & normaltable$Signal_B<0)] <- NA
    
    
    colnames(normaltable)[1]="Probe.Set.ID"

  tumortable=merge(tumortable,snp_anno,by="Probe.Set.ID")
  normaltable=merge(normaltable,snp_anno,by="Probe.Set.ID")
  
  tumortable=tumortable[order(tumortable$Chromosome,tumortable$Physical.Position),]
  normaltable=normaltable[order(normaltable$Chromosome,normaltable$Physical.Position),]
  
  
  tumortable$Physical.Position=as.integer(tumortable$Physical.Position)
  normaltable$Physical.Position=as.integer(normaltable$Physical.Position)
  
  
  PSdata = data.frame(cbind(tumortable$Chromosome,tumortable$Physical.Position,tumortable$Signal_A+tumortable$Signal_B,tumortable$Signal_B/(tumortable$Signal_A+tumortable$Signal_B), normaltable$Signal_A+normaltable$Signal_B,normaltable$Signal_B/(normaltable$Signal_A+normaltable$Signal_B)))
  names(PSdata) = c("chromosome","x","CT","betaT","CN","betaN")
  PSdata <- PSdata[!is.na(PSdata$x) & PSdata$chromosome != "X" & PSdata$chromosome != "Y",]
  for (i in 3:ncol(PSdata)) PSdata[,i] <- as.numeric(as.character(PSdata[,i]))
  
  PSdata[,1] <- as.integer(as.numeric(as.character(PSdata[,1])))
  PSdata[,2] <- as.integer(as.numeric(as.character(PSdata[,2])))
  
  
  PSdata$betaT[!is.na(PSdata$betaT) & PSdata$betaT==-Inf] <- 0
  PSdata$betaT[!is.na(PSdata$betaT) & PSdata$betaT==+Inf] <- 1
  
  
  PSdata$betaN[!is.na(PSdata$betaN) & PSdata$betaN==-Inf] <- 0
  PSdata$betaN[!is.na(PSdata$betaN) & PSdata$betaN==+Inf] <- 1
  
  PSdata <- dropSegmentationOutliers(PSdata)
  
  gaps<- findLargeGaps(PSdata, minLength = 1e+06)
  
  knownSegments <- gapsToSegments(gaps)
  
  fit <- segmentByPairedPSCBS(PSdata, knownSegments = knownSegments,preserveScale = FALSE, seed = 48879, verbose =F)
  
  deltaAB <- estimateDeltaAB(fit, scale=1)
  fit <- callAB(fit, verbose = -10)
  
  
  deltaLOH <- estimateDeltaLOH(fit, midpoint=1/2)
  fit <- callLOH(fit, delta=deltaLOH)
  
  outseg <- fit$output
  
  outseg <- data.frame(cbind(mysample,outseg))
  names(outseg)[1] <- "TCGA_ID"
  
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.PSCBS.segment.txt")
  write.table(outseg,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  }
}

  


##sum(outseg$dhEnd[!is.na(outseg$lohCall) & outseg$lohCall==TRUE] - outseg$dhStart[!is.na(outseg$lohCall) & outseg$lohCall==TRUE])/10^6
##sum(outseg$tcnEnd[!is.na(outseg$tcnEnd)] - outseg$tcnStart[!is.na(outseg$tcnStart)])/10^6
##outseg$tcnMean





### now try one example using PSCBS without normal sample ###
##i <- which(is.na(idxnormal1))[1]
##mysample=allsamples[i]
##tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[i]],"/",filelist$filename[idxtumor[i]]),
##                      skip = 1, sep="\t",header=T,stringsAsFactors = F)
##tumortable$Signal_A[which(tumortable$Signal_A<0)]=0
##tumortable$Signal_B[which(tumortable$Signal_B<0)]=0
##colnames(tumortable)[1]="Probe.Set.ID"
##tumortable=merge(tumortable,snp_anno,by="Probe.Set.ID")
##tumortable=tumortable[order(tumortable$Chromosome,tumortable$Physical.Position),]
##tumortable$Physical.Position=as.integer(tumortable$Physical.Position)
##PSdata = data.frame(cbind(tumortable$Chromosome,tumortable$Physical.Position,tumortable$Signal_A+tumortable$Signal_B,tumortable$Signal_B/(tumortable$Signal_A+tumortable$Signal_B)))
##names(PSdata) = c("chromosome","x","CT","betaT")
##PSdata <- PSdata[!is.na(PSdata$x) & PSdata$chromosome != "X" & PSdata$chromosome != "Y",]
##for (i in 3:ncol(PSdata)) PSdata[,i] <- as.numeric(as.character(PSdata[,i]))
##PSdata[,1] <- as.integer(as.numeric(as.character(PSdata[,1])))
##PSdata[,2] <- as.integer(as.numeric(as.character(PSdata[,2])))
##PSdata$betaT[!is.na(PSdata$betaT) & PSdata$betaT==-Inf] <- 0
##PSdata$betaT[!is.na(PSdata$betaT) & PSdata$betaT==+Inf] <- 1
##library(PSCBS)
##gaps<- findLargeGaps(PSdata, minLength = 1e+06)
##knownSegments <- gapsToSegments(gaps)
##fit <- segmentByNonPairedPSCBS(PSdata, knownSegments = knownSegments,preserveScale = FALSE, seed = 48879, verbose = -10)















