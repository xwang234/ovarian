#!/usr/bin/env Rscript
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
#salloc -t 1-1 -n 100 mpirun -n 1 ./compute_pvalues_cor_mrna_copynumber.R

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")

#only include samples having platinum classification
library(gdata)
#####the platinum status was mannually checked and updated.
comparison_reorder2=read.xls("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_table1.xlsx",1,stringsAsFactors=F)
colnames(comparison_reorder2)=as.character(comparison_reorder2[1,])
comparison_reorder2=comparison_reorder2[2:nrow(comparison_reorder2),]
comparison_reorder2$sample=as.character(comparison_reorder2$sample)
sampleswithdef=comparison_reorder2$sample[which(comparison_reorder2$corrected_classification1 %in% c("Resistant","Sensitive"))]
idx=which(colnames(copynumber_comcm) %in% sampleswithdef)
copynumber_comcm=copynumber_comcm[,idx]
mrna_comcm=mrna_comcm[,idx]

selgenes=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genes_smallqvalues",header=F)
cor_mrna_copynumber=data.frame(gene=selgenes[,1],correlation=rep(NA,nrow(selgenes)))
for (i in 1:nrow(selgenes))
{
  gene=selgenes[i,1]
  idx=which(rownames(copynumber_comcm)==gene)
  if (length(idx)>0)
  {
    cor_mrna_copynumber[i,2]=cor(as.numeric(mrna_comcm[idx,]),as.numeric(copynumber_comcm[idx,]),use="complete")
  }
}

computecorrelation=function(copynumber_row,mrna_row)
{
  res=NA
  res=cor(as.numeric(copynumber_row),as.numeric(mrna_row),use="complete")
  return(res)
}

computecorrelation1=function(idx,copynumber_comcm,mrna_comcm)
{
  res=NA
  res=cor(as.numeric(copynumber_comcm[idx,]),as.numeric(mrna_comcm[idx,]),use="complete")
  return(res)
}

mpi_computecorrelation=function(mrna_comcm,copynumber_comcm,filename=NULL,njob=100)
{
  library(Rmpi)
  mpi.spawn.Rslaves(needlog = FALSE)
  mpi.bcast.Robj2slave(copynumber_comcm)
  mpi.bcast.Robj2slave(mrna_comcm)
  mpi.bcast.Robj2slave(computecorrelation1)
  
  numit=1000
  res2 <- NULL
  for (i in 1:numit)
  {
    set.seed(i)
    idx_mrna=sample(ncol(mrna_comcm))
    set.seed(i+1000)
    idx_copynumber=sample(ncol(copynumber_comcm))
    res1=NULL
    nrun <- ceiling(nrow(copynumber_comcm)/1000)
    for (k in 1:nrun)
    {
        if (k < nrun) cseq <- ((k-1)*1000+1):(k*1000)  else  cseq <- ((k-1)*1000+1):nrow(copynumber_comcm)
        z=cseq
        res=mpi.parSapply(X=z,FUN=computecorrelation1,copynumber_comcm=copynumber_comcm[,idx_copynumber],mrna_comcm=mrna_comcm[,idx_mrna],job.num=njob)
        res1=c(res1,res)
    }
    res2=rbind(res2,res1)
    if ( i %% 10 ==0)
    {
      cat(i,"..")
      write.table(res2,file=filename,row.names=F,col.names=T,sep="\t",quote=F)
    }
  }
  write.table(res2,file=filename,row.names=F,col.names=T,sep="\t",quote=F)
  mpi.close.Rslaves()
  mpi.quit()
}

mpi_computecorrelation(mrna_comcm,copynumber_comcm,filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/cor_mrna_copynumber_permutation.txt",njob=100)
