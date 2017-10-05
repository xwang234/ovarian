#!/usr/bin/env Rscript

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25.RData")
TCGA_cytoband_loh=res1

computepermutation=function(jobn)
{
  set.seed(jobn+1000)
  #shuffle the outcome
  resist$resistance=resist$resistance[sample(length(resist$resistance))]
  res=data.frame(matrix(NA,nrow=nrow(TCGA_cytoband_loh),ncol=1))
  colnames(res)="p"
  rownames(res)=rownames(TCGA_cytoband_loh)
  
  for (i in 1:nrow(TCGA_cytoband_loh)){
    temp <- matrix(NA,ncol(TCGA_cytoband_loh),2)
    temp[,1] <- as.vector(colnames(TCGA_cytoband_loh))
    temp[,2] <- as.character(TCGA_cytoband_loh[i,])
    
    temp2 <- data.frame(temp)
    names(temp2) <- c("sampleID","loh")
    temp2$loh <- as.numeric(as.character(temp2$loh))
    
    if (mean(temp2$loh!=0)>0.05) {
      resist1 <- merge(resist,temp2,by="sampleID")
      fit2 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
      res[i,1] <- summary(fit2)$coeff[2,4]
    }
  }
  return(res)
}

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(TCGA_cytoband_loh)
mpi.bcast.Robj2slave(resist)
res=mpi.parSapply(X=1:1000,FUN=computepermutation,job.num=njobs)
pvalues_permutation=matrix(unlist(res),byrow =F, nrow=nrow(TCGA_cytoband_loh))

compute_aqvalue_permutation=function(pvalue,pvalues,pvalues_permutation,numit=1000)
{
  #numit=nrow(pvalues_permutation)/nrow(pvalues)
  qvalue=sum(pvalues_permutation<=pvalue,na.rm=T)/numit/sum(pvalues<=pvalue,na.rm=T)
}
mpi.bcast.Robj2slave(compute_aqvalue_permutation)

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_loh_cytoband_0.25_pvalue.RData")
pvalues=cytobandres$p
names(pvalues)=rownames(TCGA_cytoband_loh)
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
    filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_qvalues_",outputprefix,"_",numit,"permutations.txt")
    write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
  }
  print("done")
  return(qvalues)
}

res2=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=pvalues_permutation,outputprefix="tcga",numit=1000)

