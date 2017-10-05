#!/usr/bin/env Rscript
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_gene_copynumber_loh.RData")


compute6pvalues=function(jobn,resist=resist,TCGA_gene_copynumber=TCGA_gene_copynumber,TCGA_gene_loh=TCGA_gene_loh)
{
  set.seed(jobn+1000)
  #shuffle the outcome
  resist$resistance=resist$resistance[sample(length(resist$resistance))]
  outp <- matrix(NA,nrow(TCGA_gene_copynumber),6)
  
  for (i in 1:nrow(TCGA_gene_copynumber)){
    if ((i %% 1000)==0) cat(i,"..")
    temp <- matrix(NA,ncol(TCGA_gene_copynumber),2)
    temp[,1] <- as.vector(names(TCGA_gene_copynumber))
    temp[,2] <- as.character(TCGA_gene_copynumber[i,])
    
    temp1 <- data.frame(temp)
    names(temp1) <- c("sampleID","tcn")
    
    temp1$tcn <- as.numeric(as.character(temp1$tcn))
    
    if (mean(temp1$tcn!=0)>0.05) {
      resist1 <- merge(resist,temp1,by="sampleID")
      fit1 <- glm(I(resistance=="Sensitive")~tcn,family=binomial,data=resist1,x=T,y=T)
      outp[i,1] <- summary(fit1)$coeff[2,4]
    }
    
    if (mean(temp1$tcn>0.5|temp1$tcn<(-0.4))>0.05) {
      resist1 <- merge(resist,temp1,by="sampleID")
      resist1$cn <- 0
      resist1$cn <- ifelse(resist1$tcn>0.5,1,resist1$cn)
      resist1$cn <- ifelse(resist1$tcn<(-0.4),-1,resist1$cn)
      fit3 <- glm(I(resistance=="Sensitive")~cn,family=binomial,data=resist1,x=T,y=T)
      outp[i,4] <- summary(fit3)$coeff[2,4]
    }
    
    temp <- matrix(NA,ncol(TCGA_gene_loh),2)
    temp[,1] <- as.vector(names(TCGA_gene_loh))
    temp[,2] <- as.character(TCGA_gene_loh[i,])
    
    temp2 <- data.frame(temp)
    names(temp2) <- c("sampleID","loh")
    temp2$loh <- as.numeric(as.character(temp2$loh))
    
    if (mean(temp2$loh!=0)>0.05) {
      resist1 <- merge(resist,temp2,by="sampleID")
      fit2 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
      outp[i,2] <- summary(fit2)$coeff[2,4]
    }
    
    if (mean(temp2$loh>0.2)>0.05) {
      resist1 <- merge(resist,temp2,by="sampleID")
      resist1$loh <- 1*(resist1$loh>0.2)
      fit4 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
      outp[i,5] <- summary(fit4)$coeff[2,4]
    }
    
    if (mean(temp2$loh!=0)>0.05  & mean(temp1$tcn!=0)>0.05) {  
      ss1 <- fit1$x*(fit1$y-fit1$fitted)
      ss2 <- fit2$x*(fit2$y-fit2$fitted)
      
      ss <- cbind(ss1,ss2)
      bread <- matrix(0,4,4)
      bread[1:2,1:2] <- t(fit1$x) %*% (fit1$x*(1-fit1$fitted)*fit1$fitted)
      bread[3:4,3:4] <- t(fit2$x) %*% (fit2$x*(1-fit2$fitted)*fit2$fitted)
      
      beef <- t(ss)%*%ss
      
      bcov <- solve(bread) %*% beef %*% solve(bread)
      
      bb <- c(fit1$coef[2],fit2$coef[2])
      outp[i,3] <- 1-pchisq(drop(bb %*% solve(bcov[c(2,4),c(2,4)]) %*% bb),2)
      
    }
    
    if (mean(temp1$tcn>0.5|temp1$tcn<(-0.4))>0.05 & mean(temp2$loh>0.2)>0.05) {  
      ss1 <- fit3$x*(fit3$y-fit3$fitted)
      ss2 <- fit4$x*(fit4$y-fit4$fitted)
      
      ss <- cbind(ss1,ss2)
      bread <- matrix(0,4,4)
      bread[1:2,1:2] <- t(fit3$x) %*% (fit3$x*(1-fit3$fitted)*fit3$fitted)
      bread[3:4,3:4] <- t(fit4$x) %*% (fit4$x*(1-fit4$fitted)*fit4$fitted)
      
      beef <- t(ss)%*%ss
      bcov <- solve(bread) %*% beef %*% solve(bread)
      bb <- c(fit3$coef[2],fit4$coef[2])
      outp[i,6] <- 1-pchisq(drop(bb %*% solve(bcov[c(2,4),c(2,4)]) %*% bb),2)
    }
  }
  return(outp)
}

library("Rmpi")
njobs=mpi.universe.size() - 1
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(resist)
mpi.bcast.Robj2slave(TCGA_gene_copynumber)
mpi.bcast.Robj2slave(TCGA_gene_loh)
mpi.bcast.Robj2slave(compute6pvalues)

# numpermutation=1000
# res=mpi.parSapply(X=1:numpermutation,FUN=compute6pvalues,resist=resist,TCGA_gene_copynumber=TCGA_gene_copynumber,TCGA_gene_loh=TCGA_gene_loh,job.num=numpermutation)
res2=NULL
numpermutation=100
#split into 10 iterations
for (i in 1:10)
{
  cat(i,"..")
  print(Sys.time())
  jobstart=(i-1)*numpermutation+1
  jobend=i*numpermutation
  res=mpi.parSapply(X=jobstart:jobend,FUN=compute6pvalues,resist=resist,TCGA_gene_copynumber=TCGA_gene_copynumber,TCGA_gene_loh=TCGA_gene_loh,job.num=njobs)
  res1=NULL
  for (j in 1:numpermutation)
  {
    start=(j-1)*nrow(TCGA_gene_copynumber)*6+1
    end=j*nrow(TCGA_gene_copynumber)*6
    res1=rbind(res1,matrix(res[start:end],ncol=6))
  }
  res2=rbind(res2,res1)
  Sys.time()
}
TCGA_6p_permutation=res2
save(TCGA_6p_permutation,file="TCGA_6p_permutation.RData")

mpi.close.Rslaves()
mpi.quit()

#the fdr was computed in another program
# compute_qvalue_permutation=function(pvalues,pvalues_permutation,numit=1000,thename=NULL)
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
#   }else #pvalues is a vector
#   {
#     pvalues=pvalues[order(pvalues,decreasing=T)]
#   }
#   
#   idxnoNA=which(!is.na(pvalues))
#   qvalues=rep(NA,length(pvalues))
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
#   }
#   res2=data.frame(pvalues=pvalues,qvalues=qvalues)
#   rownames(res2)=names(pvalues)
#   if (!is.null(thename))
#   {
#     filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",thename,"_permutation",numit,".txt")
#     write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
#   }
#   return(res2)
# }

