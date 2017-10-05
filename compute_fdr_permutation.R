#!/usr/bin/env Rscript
#salloc -t 0-5 -n 100 mpirun -n 1 R --interactive
#salloc -t 0-5 -n 100 mpirun -n 1 ./compute_fdr_permutation.R
njob=100

# #follow http://statweb.stanford.edu/~tibs/ftp/fdringenomics.pdf Remark B, pg 15
# compute_pvalue_permutation=function(pvalue,pvalues_permutation,numpvalue)
# {
#   #number of iterations
#   numit=1000
#   #compute pvalues using permutation
#   res=sum(pvalue>=pvalues_permutation[,1])/numpvalue/numit
#   return(res)
# }

compute_qvalue_permutation=function(pvalue,pvalues,pvalues_permutation,numit=1000)
{
  #numit=nrow(pvalues_permutation)/nrow(pvalues)
  res=sum(pvalue>=pvalues_permutation,na.rm=T)/numit/sum(pvalue>=pvalues,na.rm=T)
  return(res)
}
#parallel computing:
library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
# mpi.bcast.Robj2slave(compute_pvalue_permutation)
mpi.bcast.Robj2slave(compute_qvalue_permutation)

# mpi_compute_pvalue_permutation=function(pvalues,pvalues_permutation,outputprefix)
# {
#   numpvalues=nrow(pvalues)
#   mpi.bcast.Robj2slave(pvalues)
#   mpi.bcast.Robj2slave(pvalues_permutation)
#   mpi.bcast.Robj2slave(numpvalues)
#   res1 <- NULL
#   nrun <- ceiling(nrow(pvalues)/1000)
#   for (j in 1:nrun){
#     cat(j,"..")
#     if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(pvalues)
#     z=pvalues[cseq,1]
#     res=mpi.parSapply(X=z,FUN=compute_pvalue_permutation,pvalues_permutation=pvalues_permutation,numpvalue=numpvalues,job.num=njob)
#     names(res)=rownames(pvalues)[cseq]
#     res1=c(res1,res)
#   }
#   library(qvalue)
#   #qvalues based on the original pvalues
#   qvalues=qvalue(pvalues[,1],lambda=seq(0.1,0.95,0.01))$qvalues
#   #qvalues based on pvalues from permuatation
#   qvalues1=qvalue(res1,lambda=seq(0.1,0.95,0.01))$qvalues
#   res2=data.frame(pvalues_orig=pvalues,pvalues_permut=res1,qvalues_orig=qvalues,qvalues_permut=qvalues1)
#   rownames(res2)=names(res1)
#   #save result
#   filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,".txt")
#   write.table(res2,file=filename,row.names=T,col.names=T,quote=F)
#   print("done")
# }

mpi_compute_qvalue_permutation=function(pvalues,pvalues_permutation,outputprefix)
{
  mpi.bcast.Robj2slave(pvalues)
  mpi.bcast.Robj2slave(pvalues_permutation)
  res1 <- NULL
  nrun <- ceiling(nrow(pvalues)/1000)
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(pvalues)
    z=pvalues[cseq,1]
    res=mpi.parSapply(X=z,FUN=compute_qvalue_permutation,pvalues=pvalues,pvalues_permutation=pvalues_permutation,job.num=njob)
    names(res)=rownames(pvalues)[cseq]
    res1=c(res1,res)
  }
  
  res2=data.frame(pvalues=pvalues,qvalues_permut=res1)
  rownames(res2)=names(res1)
  #save result
  filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"_1000permutation.q0.01.txt")
  write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
  print("done")
}

sapply_compute_qvalue_permutation=function(pvalues,pvalues_permutation,outputprefix)
{
  res=sapply(pvalues[,1],function(pvalue){
    res1=compute_qvalue_permutation(pvalue,pvalues,pvalues_permutation)
  })
  res2=data.frame(pvalues=pvalues,qvalues=res)
  rownames(res2)=rownames(pvalues)
  filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"1000permutation.txt")
  write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
}
# for data in different platforms:
# pvalues_copynumber_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
# pvalues_copynumber_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1000permutation.txt",header=T)
# mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_platinum,pvalues_permutation=pvalues_copynumber_platinum_permutation,
#                                 outputprefix="copynumber")
# # sapply_compute_qvalue_permutation(pvalues=pvalues_copynumber_platinum,pvalues_permutation=pvalues_copynumber_platinum_permutation,
# #                                  outputprefix="copynumber_platinum")
# pvalues_mrna_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum.txt",header=T)
# pvalues_mrna_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_1000permutation.txt",header=T)
# mpi_compute_qvalue_permutation(pvalues=pvalues_mrna_platinum,pvalues_permutation=pvalues_mrna_platinum_permutation,
#                                 outputprefix="mrna")
# # sapply_compute_qvalue_permutation(pvalues=pvalues_mrna_platinum,pvalues_permutation=pvalues_mrna_platinum_permutation,
# #                                outputprefix="mrna_platinum")

# pvalues_mrna_exon_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_exon_platinum.txt",header=T)
# pvalues_mrna_exon_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_exon_1000permutation.txt",header=T)
# mpi_compute_qvalue_permutation(pvalues=pvalues_mrna_exon_platinum,pvalues_permutation=pvalues_mrna_exon_platinum_permutation,
#                                 outputprefix="mrna_exon")
# # sapply_compute_qvalue_permutation(pvalues=pvalues_mrna_exon_platinum,pvalues_permutation=pvalues_mrna_exon_platinum_permutation,
# #                                outputprefix="mrna_exon_platinum")
# 
# pvalues_mrna_u133_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_u133_platinum.txt",header=T)
# pvalues_mrna_u133_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_u133_1000permutation.txt",header=T)
# mpi_compute_qvalue_permutation(pvalues=pvalues_mrna_u133_platinum,pvalues_permutation=pvalues_mrna_u133_platinum_permutation,
#                                 outputprefix="mrna_u133")
# # sapply_compute_qvalue_permutation(pvalues=pvalues_mrna_u133_platinum,pvalues_permutation=pvalues_mrna_u133_platinum_permutation,
# #                                outputprefix="mrna_u133_platinum")

pvalues_copynumber_CGH1M_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_CGH1M_platinum.txt",header=T)
pvalues_copynumber_CGH1M_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_CGH1M_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_CGH1M_platinum,pvalues_permutation=pvalues_copynumber_CGH1M_platinum_permutation,
                                outputprefix="copynumber_CGH1M")
# sapply_compute_qvalue_permutation(pvalues=pvalues_copynumber_CGH1M_platinum,pvalues_permutation=pvalues_copynumber_CGH1M_platinum_permutation,
#                                outputprefix="copynumber_CGH1M_platinum")

pvalues_copynumber_1MDUO_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1MDUO_platinum.txt",header=T)
pvalues_copynumber_1MDUO_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1MDUO_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_1MDUO_platinum,pvalues_permutation=pvalues_copynumber_1MDUO_platinum_permutation,
                               outputprefix="copynumber_1MDUO")
# sapply_compute_qvalue_permutation(pvalues=pvalues_copynumber_1MDUO_platinum,pvalues_permutation=pvalues_copynumber_1MDUO_platinum_permutation,
#                                outputprefix="copynumber_1MDUO_platinum")

pvalues_4class_copynumber_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_4class_copynumber_platinum.txt",header=T)
pvalues_4class_copynumber_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_4class_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_4class_copynumber_platinum,pvalues_permutation=pvalues_4class_copynumber_platinum_permutation,
                               outputprefix="4class_copynumber")

pvalues_4class_mrna_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_4class_mrna_platinum.txt",header=T)
pvalues_4class_mrna_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_4class_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_4class_mrna_platinum,pvalues_permutation=pvalues_4class_mrna_platinum_permutation,
                               outputprefix="4class_mrna")

pvalues_4class_copynumber_rmborderline_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_4class_copynumber_rmborderline_platinum.txt",header=T)
pvalues_4class_copynumber_rmborderline_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_rmborderline_4class_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_4class_copynumber_rmborderline_platinum,pvalues_permutation=pvalues_4class_copynumber_rmborderline_platinum_permutation,
                               outputprefix="4class_copynumber_rmborderline")

pvalues_4class_mrna_rmborderline_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_4class_mrna_rmborderline_platinum.txt",header=T)
pvalues_4class_mrna_rmborderline_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_rmborderline_4class_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_4class_mrna_rmborderline_platinum,pvalues_permutation=pvalues_4class_mrna_rmborderline_platinum_permutation,
                               outputprefix="4class_mrna_rmborderline")

pvalues_copynumber_extreme_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_extreme_platinum.txt",header=T)
pvalues_copynumber_extreme_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_extreme_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_extreme_platinum,pvalues_permutation=pvalues_copynumber_extreme_platinum_permutation,
                               outputprefix="copynumber_extreme")

pvalues_mrna_extreme_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_extreme_platinum.txt",header=T)
pvalues_mrna_extreme_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_extreme_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_mrna_extreme_platinum,pvalues_permutation=pvalues_mrna_extreme_platinum_permutation,
                               outputprefix="mrna_extreme")

pvalues_copynumber_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_filtered_platinum.txt",header=T)
pvalues_copynumber_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_filtered_1000permutation.txt",header=T)
pvalues_copynumber_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_filtered_platinum,pvalues_permutation=pvalues_copynumber_filtered_platinum_permutation,
                               outputprefix="copynumber_filtered")
pvalues=pvalues_copynumber_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_filtered_platinum)
fdrs=read.table(file="../result/qvalues_copynumber_filtered_1000permutation.txt")
fdrs1=read.table(file="../result/qvalues_copynumber_1000permutation.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)
#[1] 336

pvalues_copynumber_tangent_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_platinum.txt",header=T)
pvalues_copynumber_tangent_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_tangent_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_filtered")
pvalues=pvalues_copynumber_tangent_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_tangent_filtered_platinum)
fdrs=read.table(file="../result/qvalues_copynumber_tangent_filtered_1000permutation.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)

pvalues_copynumber_tangent_gistic_filtered=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_gistic_filtered_platinum.txt",header=T)

pvalues_copynumber_tangent_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_tangent_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_filtered")
pvalues=pvalues_copynumber_tangent_gistic_filtered[,1]
names(pvalues)=rownames(pvalues_copynumber_tangent_gistic_filtered)
draw_manhattan(pvalues=pvalues)
fdrs=read.table(file="../result/qvalues_copynumber_tangent_filtered_1000permutation.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)


pvalues_copynumber_aocs_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_aocs_filtered_platinum.txt",header=T)
pvalues_copynumber_aocs_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_aocs_filtered_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_aocs_filtered_platinum,pvalues_permutation=pvalues_copynumber_aocs_filtered_platinum_permutation,
                               outputprefix="copynumber_aocs_filtered")
pvalues=pvalues_copynumber_aocs_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_aocs_filtered_platinum)
draw_manhattan(pvalues=pvalues)
fdrs=read.table(file="../result/qvalues_copynumber_aocs_filtered_1000permutation.txt")
draw_manhattan(fdrs=fdrs)

pvalues_copynumber_tangent_gistic_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_gistic_filtered_platinum.txt",header=T)
pvalues_copynumber_tangent_gistic_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_gistic_filtered_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_tangent_gistic_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_gistic_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_gistic_filtered")
fdrs=read.table(file="../result/qvalues_copynumber_tangent_gistic_filtered_1000permutation.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)
#[1] 279

#segment from firehose, redo gistic
pvalues_copynumber_gistic_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_gistic_platinum.txt",header=T)
pvalues_copynumber_gistic_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_gistic_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_gistic_platinum,pvalues_permutation=pvalues_copynumber_gistic_platinum_permutation,
                               outputprefix="copynumber_gistic")
fdrs=read.table(file="../result/qvalues_copynumber_gistic_1000permutation.txt")
draw_manhattan(fdrs=fdrs)

#based on firehose segment
pvalues_copynumber_firehosesegment_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_firehosesegment_platinum.txt",header=T)
pvalues_copynumber_firehosesegment_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_firehosesegment_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_firehosesegment_platinum,pvalues_permutation=pvalues_copynumber_firehosesegment_platinum_permutation,
                               outputprefix="copynumber_firehosesegment")
fdrs=read.table(file="../result/qvalues_copynumber_firehosesegment_1000permutation.txt")
draw_manhattan(fdrs=fdrs)

#q<=0.01
pvalues_copynumber_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_filtered_platinum.txt",header=T)
pvalues_copynumber_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_filtered_1000permutation.q0.01.txt",header=T)
pvalues_copynumber_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1000permutation.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_filtered_platinum,pvalues_permutation=pvalues_copynumber_filtered_platinum_permutation,
                               outputprefix="copynumber_filtered")
pvalues=pvalues_copynumber_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_filtered_platinum)
fdrs=read.table(file="../result/qvalues_copynumber_filtered_1000permutation.q0.01.txt")
fdrs1=read.table(file="../result/qvalues_copynumber_1000permutation.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)
#[1] 336

pvalues_copynumber_tangent_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_platinum.txt",header=T)
pvalues_copynumber_tangent_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_1000permutation.q0.01.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_tangent_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_filtered")
pvalues=pvalues_copynumber_tangent_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_tangent_filtered_platinum)
fdrs=read.table(file="../result/qvalues_copynumber_tangent_filtered_1000permutation.q0.01.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)

pvalues_copynumber_tangent_gistic_filtered=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_gistic_filtered_platinum.txt",header=T)
pvalues_copynumber_tangent_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_filtered_1000permutation.q0.01.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_tangent_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_filtered")
pvalues=pvalues_copynumber_tangent_gistic_filtered[,1]
names(pvalues)=rownames(pvalues_copynumber_tangent_gistic_filtered)
draw_manhattan(pvalues=pvalues)
fdrs=read.table(file="../result/qvalues_copynumber_tangent_filtered_1000permutation.q0.01.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)


pvalues_copynumber_aocs_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_aocs_filtered_platinum.txt",header=T)
pvalues_copynumber_aocs_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_aocs_filtered_1000permutation.q0.01.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_aocs_filtered_platinum,pvalues_permutation=pvalues_copynumber_aocs_filtered_platinum_permutation,
                               outputprefix="copynumber_aocs_filtered")
pvalues=pvalues_copynumber_aocs_filtered_platinum[,1]
names(pvalues)=rownames(pvalues_copynumber_aocs_filtered_platinum)
draw_manhattan(pvalues=pvalues)
fdrs=read.table(file="../result/qvalues_copynumber_aocs_filtered_1000permutation.q0.01.txt")
draw_manhattan(fdrs=fdrs)

pvalues_copynumber_tangent_gistic_filtered_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_gistic_filtered_platinum.txt",header=T)
pvalues_copynumber_tangent_gistic_filtered_platinum_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_tangent_gistic_filtered_1000permutation.q0.01.txt",header=T)
mpi_compute_qvalue_permutation(pvalues=pvalues_copynumber_tangent_gistic_filtered_platinum,pvalues_permutation=pvalues_copynumber_tangent_gistic_filtered_platinum_permutation,
                               outputprefix="copynumber_tangent_gistic_filtered")
fdrs=read.table(file="../result/qvalues_copynumber_tangent_gistic_filtered_1000permutation.q0.01.txt")
draw_manhattan(fdrs=fdrs)
sum(fdrs[,2]<=0.05)
#[1] 279

mpi.close.Rslaves()
mpi.quit()