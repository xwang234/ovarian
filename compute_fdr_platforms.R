#!/usr/bin/env Rscript
library(qvalue)
source("./functions.R")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")
compute_fdrs=function(data=data_copynumber_CGH1M,platform="copynumber_CGH1M")
{
  Sys.time()
  print("compute pvalues...")
  pvalues=compute_pvalues(data=data,platform=platform,runpermutation=T)
  Sys.time()
 
  #read the permutation pvalues genreated in the previous step
  pvalues_permutation=read.table(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,"_permutation1.txt"),header=F)
  Sys.time()
  print("compute fdr...")
  fdrs=sapply_compute_qvalue_permutation(pvalues=pvalues,pvalues_permutation=pvalues_permutation,
                                                    platform=platform)
  sum(fdrs<0.05)
  Sys.time()
  smallfdr=smallfdrtable(fdrs=fdrs,data=data,threshold=0.05)
  result=list(fdrs=fdrs,smallfdr=smallfdr)
  return(result)
}
compute_fdrs(data=data_copynumber_CGH1M,platform="copynumber_CGH1M")
compute_fdrs(data=data_copynumber_1MDUO,platform="copynumber_1MDUO")
platform="copynumber_CGH1M"
platform="copynumber_1MDUO"
platform="copynumber"
platform="mrna"
platform="copynumber_rmborderline"
platform="mrna_rmborderline"
platform="copynumber_extreme"
platform="mrna_extreme"
fdrs=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",platform,"_1000permutation.txt"),header=T)
fdrs=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_4class_",platform,"_1000permutation.txt"),header=T)
# genenames=rownames(fdrs)
# idx=!is.na(rowSums(fdrs))
# pvalues=fdrs[idx,1]
# names(pvalues)=genenames[idx]
# 
# fdrs=fdrs[idx,2]
# names(fdrs)=genenames[idx]
draw_manhattan(fdrs=fdrs,fdrthreshold=fdrthreshold)

fdrs1=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",platform,"_1000permutation.txt"),header=T)
qvalues=qvalue(pvalues,lambda=seq(0.1,0.95,0.01))$qvalues
draw_manhattan(fdrs=fdrs1,fdrthreshold=fdrthreshold)

#for rmborderlines
fdrs1=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_4class_",platform,"_1000permutation.txt"),header=T)
hist(fdrs1[,"pvalues"],main=paste0("Histogram of ",platform),probability = T)
draw_manhattan(fdrs=fdrs1,fdrthreshold=fdrthreshold)


smallfdr=smallfdrtable(fdrs=fdrs1,data=data4_mrna,threshold=0.05)
write.table(smallfdr,file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/",platform,"_result_smallfdr.txt"),row.names = F,col.names = T,sep="\t",quote=F)





