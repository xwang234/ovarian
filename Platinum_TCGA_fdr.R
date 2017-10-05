#!/usr/bin/env Rscript

#salloc -t 1-5 -n 15 mpirun -n 1 R --interactive
#load permutation results
TCGA_6p_permutation=NULL
for (i in 1:10)
{
  load(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_fdr",i,".RData"))
  TCGA_6p_permutation=rbind(TCGA_6p_permutation,res1)
}
#save(TCGA_6p_permutation,"/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_6p_permutation.RData")
# too slow
# compute_qvalue_permutation1=function(pvalues,pvalues_permutation,numit=1000,thename=NULL)
# {
#   #pvalues is observed p-values 
#   #It can be a 1-column dataframe
#   #sort pvalues decreasingly
#   if (is.data.frame(pvalues))
#   {
#     genenames=rownames(pvalues)
#     idx=order(pvalues[,1],decreasing=T)
#     pvalues=pvalues[idx,1]
#     names(pvalues)=genenames[idx]
#   }else #pvalues is a named vector
#   {
#     genenames=names(pvalues)
#     pvalues=pvalues[order(pvalues,decreasing=T)]
#   }
#   
#   idxnoNA=which(!is.na(pvalues))
#   qvalues=rep(NA,length(pvalues))
#   #for (i in idxnoNA)
#   Sys.time()
#   for (i in idxnoNA)
#   {
#     #compute q-value
#     qvalue=sum(pvalues_permutation<=pvalues[i],na.rm=T)/numit/sum(pvalues<=pvalues[i],na.rm=T)
#     #make q values monotone
#     if (i==idxnoNA[1])
#     {
#       #keep the first q value associated with the highest p-value
#       qvalues[i]=qvalue
#     }else
#     {
#       #make sure q-value is non-increasing
#       qvalues[i]=min(qvalues[i-1],qvalue,na.rm=T)
#     }
#     if (i %% 100==0) cat(i,"..")
#   }
#   
#   
#   #use the original order of p-values
#   idx=match(names(pvalues),genenames)
#   oqvalues=qvalues[idx]
#   names(oqvalues)=genenames
#   
#   res2=data.frame(pvalues=pvalues,qvalues=qvalues)
#   rownames(res2)=names(pvalues)
#   if (!is.null(thename))
#   {
#     filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",thename,"_permutation",numit,".txt")
#     write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
#   }
#   return(oqvalues)
# }

compute_aqvalue_permutation=function(pvalue,pvalues,pvalues_permutation,numit=1000)
{
  #numit=nrow(pvalues_permutation)/nrow(pvalues)
  qvalue=sum(pvalues_permutation<=pvalue,na.rm=T)/numit/sum(pvalues<=pvalue,na.rm=T)
}

library("Rmpi")
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(compute_aqvalue_permutation)

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






#load the observed p-values
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_platinum_cnloh.RData")
#compute and store the q-values
for (i in 1:6)
{
  print(i)
  pvalues=outp[,i]
  names(pvalues)=rownames(outp)
  res=mpi_compute_qvalue_permutation(pvalues,pvalues_permutation=TCGA_6p_permutation[,i],outputprefix=paste0("TCGA",i),numit=1000)
  print(min(res,na.rm=T))
}
mpi.close.Rslaves()
mpi.quit()
