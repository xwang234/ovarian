#!/usr/bin/env Rscript
#salloc -t 0-5 -n 100 mpirun -n 1 R --interactive

njob=100
compute_fwer_permutation=function(pvalue,pvalues_permutation)
{
  numit=nrow(pvalues_permutation)/nrow(pvalues)
  pvalues_permutation_matrix=matrix(pvalues_permutation[,1],ncol=numit)
  res=sapply(1:numit,function(i){
    res1=min(pvalues_permutation_matrix[,i],na.rm=T)<=pvalue
  })
  res1=sum(res)/numit
  return(res1)
}

#parallel computing:
library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
# mpi.bcast.Robj2slave(compute_pvalue_permutation)
mpi.bcast.Robj2slave(compute_fwer_permutation)

mpi_compute_fwer_permutation=function(pvalues,pvalues_permutation,outputprefix)
{
  mpi.bcast.Robj2slave(pvalues)
  mpi.bcast.Robj2slave(pvalues_permutation)
  res1 <- NULL
  nrun <- ceiling(nrow(pvalues)/1000)
  print(paste0("total number of runs: ",nrun))
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(pvalues)
    z=pvalues[cseq,1]
    res=mpi.parSapply(X=z,FUN=compute_fwer_permutation,pvalues_permutation=pvalues_permutation,job.num=njob)
    names(res)=rownames(pvalues)[cseq]
    res1=c(res1,res)
  }
  
  res2=data.frame(pvalues=pvalues,fwer=res1)
  rownames(res2)=names(res1)
  #save result
  filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/fwer_",outputprefix,"_1000permutation.txt")
  write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
  print("done")
}


pvalues_copynumber_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1000permutation.txt",header=T)
pvalues_copynumber=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
outputprefix="copynumber"
mpi_compute_fwer_permutation(pvalues=pvalues_copynumber,pvalues_permutation=pvalues_copynumber_permutation,outputprefix="copynumber")
filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/fwer_",outputprefix,"_1000permutation.txt")
fwer_copynumber=read.table(file=filename,header=T,sep="\t")
pvalues=pvalues_copynumber[,1]
names(pvalues)=rownames(pvalues_copynumber)
draw_manhattan(pvalues=pvalues)
#to find pvalue with fwer=0.05
idx=order(fwer_copynumber[,2])
fwer_copynumber=fwer_copynumber[idx,]
idx=fwer_copynumber[,2]<=0.05
idx=which(idx==T)
threshold=fwer_copynumber[length(idx),1]
abline(h=-log10(threshold),lwd=2,col="red")

pvalues_copynumber_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_filtered_platinum.txt",header=T)
pvalues_copynumber_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_filtered_1000permutation.txt",header=T)
mpi_compute_fwer_permutation(pvalues=pvalues_copynumber_filtered_platinum,pvalues_permutation=pvalues_copynumber_filtered_platinum_permutation,
                               outputprefix="copynumber_filtered")
outputprefix="copynumber_filtered"
filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/fwer_",outputprefix,"_1000permutation.txt")
fwer_copynumber=read.table(file=filename,header=T,sep="\t")
pvalues=pvalues_copynumber_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_filtered_platinum)
draw_manhattan(pvalues=pvalues)
#to find pvalue with fwer=0.05
idx=order(fwer_copynumber[,2])
fwer_copynumber=fwer_copynumber[idx,]
idx=fwer_copynumber[,2]<=0.1
idx=which(idx==T)
threshold=fwer_copynumber[length(idx),1]
abline(h=-log10(threshold),lwd=2,col="red")

pvalues_copynumber_tangent_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_platinum.txt",header=T)
pvalues_copynumber_tangent_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_1000permutation.txt",header=T)
mpi_compute_fwer_permutation(pvalues=pvalues_copynumber_tangent_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_filtered")
outputprefix="copynumber_tangent_filtered"
filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/fwer_",outputprefix,"_1000permutation.txt")
fwer_copynumber=read.table(file=filename,header=T,sep="\t")
pvalues=pvalues_copynumber_tangent_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_tangent_filtered_platinum)
draw_manhattan(pvalues=pvalues)
#to find pvalue with fwer=0.05
idx=order(fwer_copynumber[,2])
fwer_copynumber=fwer_copynumber[idx,]
idx=fwer_copynumber[,2]<=0.1
idx=which(idx==T)
threshold=fwer_copynumber[length(idx),1]
abline(h=-log10(threshold),lwd=2,col="red")

pvalues_copynumber_aocs_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_aocs_filtered_platinum.txt",header=T)
pvalues_copynumber_aocs_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_aocs_filtered_1000permutation.txt",header=T)
mpi_compute_fwer_permutation(pvalues=pvalues_copynumber_aocs_filtered_platinum,pvalues_permutation=pvalues_copynumber_aocs_filtered_platinum_permutation,
                               outputprefix="copynumber_aocs_filtered")
outputprefix="copynumber_aocs_filtered"
filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/fwer_",outputprefix,"_1000permutation.txt")
fwer_copynumber=read.table(file=filename,header=T,sep="\t")
pvalues=pvalues_copynumber_aocs_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_aocs_filtered_platinum)
draw_manhattan(pvalues=pvalues)
#to find pvalue with fwer=0.05
idx=order(fwer_copynumber[,2])
fwer_copynumber=fwer_copynumber[idx,]
idx=fwer_copynumber[,2]<=0.1
idx=which(idx==T)
threshold=fwer_copynumber[length(idx),1]
abline(h=-log10(threshold),lwd=2,col="red")
