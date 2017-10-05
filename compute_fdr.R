#!/usr/bin/env Rscript
#first calculate p-value of 2df first,then permutation based p-value,then fdr

load("../data/tcga_data.RData")


two_df_pvalues=function(x,y,z)
{
  fit1 <- glm(z~x,family=binomial,x=T,y=T)
  fit2 <- glm(z~y,family=binomial,x=T,y=T)
  score1 <- fit1$x*(fit1$y-fit1$fitted)
  score2 <- fit2$x*(fit2$y-fit2$fitted)
  score <- cbind(score1,score2)
  bread <- matrix(0,4,4)
  bread[1:2,1:2] <- t(fit1$x) %*% (fit1$x*(1-fit1$fitted)*fit1$fitted)
  bread[3:4,3:4] <- t(fit2$x) %*% (fit2$x*(1-fit2$fitted)*fit2$fitted)
  varmat <- solve(bread) %*% (t(score)%*% score) %*% solve(bread)
  bb <- c(fit1$coef[2] ,fit2$coef[2])
  bbmat <- varmat[c(2,4),c(2,4)]
  tt <- drop(bb %*% solve(bbmat) %*% bb)
  pval <- 1-pchisq(tt,2)
}

makecomdata=function(data1,data2)
{
  commongenes=intersect(colnames(data1)[2:ncol(data1)],colnames(data2)[2:ncol(data2)])
  commonsamples=intersect(rownames(data1)[1:nrow(data1)],rownames(data2)[1:nrow(data2)])
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data1)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data1)==mysample)
  })
  #output
  z=data1[idx_row,1]
  data1=data1[idx_row,idx_col]
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data2)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data2)==mysample)
  })
  data2=data2[idx_row,idx_col]
  return(result=list(data1=data1,data2=data2,z=z))
}
# two_df_pvalues_data=function(data1,data2,z)
# {
#   
#   pvalues=sapply(1:ncol(data1),function(i){
#     x=data1[,i]
#     y=data2[,i]
#     res=two_df_pvalues(x=x,y=y,z=z)
#   })
#   names(pvalues)=colnames(data1)
#   ##for some cases when x or y are mostly 0, which causes the problem of p=0
#   idx_pvalue0=which(pvalues==0)
#   if (length(idx_pvalue0)>0) pvalues[idx_pvalue0]=NA
#   
#   return(pvalues)
# }

# calpvalues=function(x1,y1)
# {
#   fit=glm(y1~x1,family="binomial")
#   summaryres=coef(summary(fit))
#   #if z is constant, no pvalue would be found, and just has the result for intercept
#   if (nrow(summaryres)>1)
#   {
#     res=coef(summary(fit))[2,4]
#   }else
#   {
#     #no p-value found
#     res=NA
#   }
#   return(res)
# }
# 
# two_df_pvalues_1col=function(i,data1,data2,z1)
# {
#   x=data1[,i]
#   y=data2[,i]
#   res=two_df_pvalues(x=x,y=y,z=z1)
#   if (res==0) res=NA
#   names(res)=colnames(data1)[i]
#   return(res)
#   return(pvalues)
# }

calallpvalues=function(i,data1,data2,z)
{
  res1=res2=NA
  x=data1[,i]
  y=data2[,i]
  fit=glm(z~x,family="binomial")
  summaryres=coef(summary(fit))
  #if z is constant, no pvalue would be found, and just has the result for intercept
  if (nrow(summaryres)>1)
  {
    res1=coef(summary(fit))[2,4]
  }
  fit=glm(z~y,family="binomial")
  summaryres=coef(summary(fit))
  #if z is constant, no pvalue would be found, and just has the result for intercept
  if (nrow(summaryres)>1)
  {
    res2=coef(summary(fit))[2,4]
  }
  res3=two_df_pvalues(x=x,y=y,z=z)
  if (res3==0) res3=NA
  return(c(res1,res2,res3))
}

data1=data_copynumber
data2=data_copynumber_allele
tmp=makecomdata(data1,data2)
data1=tmp$data1
data2=tmp$data2
z=tmp$z

library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
mpi.bcast.Robj2slave(two_df_pvalues)
# mpi.bcast.Robj2slave(two_df_pvalues_1col)
mpi.bcast.Robj2slave(calpvalues)
mpi.bcast.Robj2slave(data1)
mpi.bcast.Robj2slave(data2)
mpi.bcast.Robj2slave(calallpvalues)

mpi_calpvalues=function(data1,data2,z1,njob=100)
{
  mpi.bcast.Robj2slave(z1)
  
  res=data.frame(matrix(NA,nrow=ncol(data1),ncol=3))
  colnames(res)=c("cn","dh","twodf")
  rownames(res)=colnames(data1)
  
  nrun <- ceiling(ncol(data1)/1000)
  for (j in 1:nrun){
    #cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):ncol(data1)
#     x=data1[,cseq]
#     res1=mpi.parCapply(X=x,FUN=calpvalues,y1=z1,job.num=njob)
#     res$cn[cseq]=res1
#     y=data2[,cseq]
#     res2=mpi.parCapply(X=y,FUN=calpvalues,y1=z1,job.num=njob)
#     res$dh[cseq]=res2
#     res3=mpi.parSapply(X=cseq,FUN=two_df_pvalues_1col,data1=data1,data2=data2,z1=z1,job.num=njob)
#     res$twodf[cseq]=res3
    res1=mpi.parSapply(X=cseq,FUN=calallpvalues,data1=data1,data2=data2,z=z1)
    idx=seq(1,length(res1),3)
    res$cn[cseq]=res1[idx]
    res$dh[cseq]=res1[idx+1]
    res$twodf[cseq]=res1[idx+2]
    
  }
  return(res)
}
#the pvalues
pvalues_tcga=mpi_calpvalues(data1,data2,z1=z,njob=100)
write.table(pvalues_tcga,file="../result/pvalues_tcga_copynumber.txt",col.names=T,row.names=T,sep="\t",quote=F)

#the permutation p values
numit=1000
pvalues_permutation_tcga=c()
for (i in 1:numit)
{
  
  set.seed(i+1000)
  z1=z[sample(1:length(z))]
  tmp=mpi_calpvalues(data1,data2,z1=z1,njob=100)
  rownames(tmp)=paste0(rownames(tmp),".",i)
  pvalues_permutation_tcga=rbind(pvalues_permutation_tcga,tmp)
  if (i %% 10==0)
  {
    write.table(pvalues_permutation_tcga,file="../result/pvalues_permutation_tcga_copynumber.txt",col.names=T,row.names=T,sep="\t",quote=F)
    cat(i,"..")
  }
}

write.table(pvalues_permutation_tcga,file="../result/pvalues_permutation_tcga_copynumber.txt",col.names=T,row.names=T,sep="\t",quote=F)

pvalues_tcga_copynumber=read.table(file="../result/pvalues_tcga_copynumber.txt",header=T,sep='\t')
pvalues_permutation_tcga_copynumber=read.table(file="../result/pvalues_permutation_tcga_copynumber.txt",header=T,sep='\t')
pvalues=data.frame(pvalues=pvalues_tcga_copynumber$cn)
rownames(pvalues)=rownames(pvalues_tcga_copynumber)
fdrs_copynumber_cn=compute_qvalue_permutation(pvalues,pvalues_permutation_tcga_copynumber$cn,numit=100,platform="copynumber_cn")
draw_manhattan(fdrs=fdrs_copynumber_cn,main="TCGA_totalCN")

pvalues=data.frame(pvalues=pvalues_tcga_copynumber$dh)
rownames(pvalues)=rownames(pvalues_tcga_copynumber)
fdrs_copynumber_dh=compute_qvalue_permutation(pvalues,pvalues_permutation_tcga_copynumber$dh,numit=100,platform="copynumber_dh")
draw_manhattan(fdrs=fdrs_copynumber_dh,main="TCGA_BAF")

pvalues=data.frame(pvalues=pvalues_tcga_copynumber$twodf)
rownames(pvalues)=rownames(pvalues_tcga_copynumber)
fdrs_copynumber_2df=compute_qvalue_permutation(pvalues,pvalues_permutation_tcga_copynumber$twodf,numit=100,platform="copynumber_2df")
draw_manhattan(fdrs=fdrs_copynumber_2df,main="TCGA_2df")

pvalues_copynumber=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
pvalues_copynumber_1000permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1000permutation.txt",header=T)
system.time()
fdrs_copynumber=compute_qvalue_permutation(pvalues=pvalues_copynumber,pvalues_permutation=pvalues_copynumber_1000permutation,platform="copynumber")
system.time()
draw_manhattan(fdrs=fdrs_copynumber,fdrthreshold=0.05,main="TCGA")

load("../data/aocs_data.RData")
data1=data_aocs_copynumber
data2=data_aocs_mean_copynumber_allele
tmp=makecomdata(data1,data2)
data1=tmp$data1
data2=tmp$data2
z=tmp$z
pvalues_aocs=mpi_calpvalues(data1,data2,z1=z,njob=100)
write.table(pvalues_aocs,file="../result/pvalues_aocs_copynumber.txt",col.names=T,row.names=T,sep="\t",quote=F)

pvalues_aocs_copynumber=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_aocs_copynumber.txt",header=T)
pvalues=data.frame(pvalues=pvalues_aocs_copynumber$cn)
rownames(pvalues)=rownames(pvalues_aocs_copynumber)
draw_manhattan(pvalues=pvalues,main="ICGC_totalCN")

pvalues=data.frame(pvalues=pvalues_aocs_copynumber$dh)
rownames(pvalues)=rownames(pvalues_aocs_copynumber)
draw_manhattan(pvalues=pvalues,main="ICGC_BAF")

pvalues=data.frame(pvalues=pvalues_aocs_copynumber$twodf)
rownames(pvalues)=rownames(pvalues_aocs_copynumber)
draw_manhattan(pvalues=pvalues,main="ICGC_2df")


pvalues_aocs_copynumber_1000permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_aocs_copynumber_1000permutation.txt",header=T)
fdrs_aocs_copynumber=compute_qvalue_permutation(pvalues=pvalues_aocs_copynumber,pvalues_permutation=pvalues_aocs_copynumber_1000permutation,numit=200,platform="aocs_copynumber")
draw_manhattan(fdrs=fdrs_aocs_copynumber,fdrthreshold=0.05,main="ICGC")