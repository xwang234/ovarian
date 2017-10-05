#!/usr/bin/env Rscript

#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
#salloc -t 3-1 -n 100 mpirun -n 1 ./mpi_pvalue_fdr_manhattan.R
#include functions
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
njob=100
#load data_copynumber:
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_4classificationdata_stringent.RData")
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_4classification_otherdata_stringent.RData")
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_copynumber.RData")

#sub-functions section begin---------------------------------------------------------------------
#function to compute fdr. pvalue is an element of pvalues array; pvalue_permutation is the null data;numit is the number of iteration
compute_qvalue_permutation=function(pvalue,pvalues,pvalues_permutation,numit=500)
{
  #numit=nrow(pvalues_permutation)/nrow(pvalues)
  res=sum(pvalue>=pvalues_permutation,na.rm=T)/numit/sum(pvalue>=pvalues,na.rm=T)
  return(res)
}

#compute p-value considering clinicaldata
calpvalues_adjustc=function(z,y,clinicaldata)
{
  tmpdata=cbind(y,z,clinicaldata)
  fit=glm(as.formula(paste0(colnames(tmpdata)[1],"~",paste0(colnames(tmpdata)[2:ncol(tmpdata)],collapse="+"))),data=tmpdata,family="binomial")
  summaryres=coef(summary(fit))
  #if z is constant, no pvalue would be found, and we just have the result for intercept
  idx_z=which(rownames(summaryres)=="z")
  if (length(idx_z)>0)
  {
    res=coef(summary(fit))[idx_z,4]
  }else
  {
    #no p-value found
    res=NA
  }
  return(res)
}


library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
#calpvalues function is in functions.R
mpi.bcast.Robj2slave(calpvalues)
mpi.bcast.Robj2slave(calpvalues_adjustc)
mpi.bcast.Robj2slave(compute_qvalue_permutation)


#parallel computing pvalues
#data is the input x matrix, y is the output y matrix; p-values are stored in filename file; opt=NULL: calculate p-values not considering clinical data
mpi_calpvalues=function(data,y,filename=NULL,opt=NULL,njob=100,clinicaldata=NULL)
{
  res1 <- NULL
  if (is.null(opt)) #if not considering clinical data
  {
    nrun <- ceiling(ncol(data)/1000)
    for (j in 1:nrun){
      #cat(j,"..")
      if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):ncol(data)
      z=data[,cseq]
      if (is.null(opt))
      {
        res=mpi.parCapply(X=z,FUN=calpvalues,y=y,job.num=njob)
      }else
      {
        #clinicaldata is global
        res=mpi.parCapply(X=z,FUN=calpvalues_adjustc,y=y,clinicaldata=clinicaldata,job.num=njob)
      }
      
      res1=c(res1,res)
    }
    res1=res1[1:ncol(data)]
    res2=data.frame(pvalues=res1)
    rownames(res2)=names(res1)
    if (!is.null(filename))
    {
      write.table(res2,file=filename,quote=F,sep="\t",row.names=T,col.names=T)
    }
  }else #if considering clinical data
  {
    nrun <- ceiling(ncol(data)/1000)
    for (j in 1:nrun){
      #cat(j,"..")
      if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):ncol(data)
      
      z=data[,cseq]
      if (is.null(opt))
      {
        res=mpi.parCapply(X=z,FUN=calpvalues,y=y,job.num=njob)
      }else
      {
        res=mpi.parCapply(X=z,FUN=calpvalues_adjustc,y=y,clinicaldata=clinicaldata,job.num=njob)
      }
      res1=c(res1,res)
    }
    res1=res1[1:ncol(data)]
    res2=data.frame(pvalues=res1)
    rownames(res2)=names(res1)
    if (!is.null(filename))
    {
      write.table(res2,file=filename,quote=F,sep="\t",row.names=T,col.names=T)
    }
  }
 
  return(res2)
}


#parallel computing for fdr
mpi_compute_qvalue_permutation=function(pvalues,pvalues_permutation,platform,numit=500)
{
  mpi.bcast.Robj2slave(pvalues)
  mpi.bcast.Robj2slave(pvalues_permutation)
  res1 <- NULL
  nrun <- ceiling(nrow(pvalues)/1000)
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(pvalues)
    z=pvalues[cseq,1]
    res=mpi.parSapply(X=z,FUN=compute_qvalue_permutation,pvalues=pvalues,pvalues_permutation=pvalues_permutation,numit=numit,job.num=njob)
    names(res)=rownames(pvalues)[cseq]
    res1=c(res1,res)
  }
  
  res2=data.frame(pvalues=pvalues,fdrs=res1)
  rownames(res2)=names(res1)
  #save fdr result
  filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",platform,"_1000permutation.txt")
  write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
  print("done")
  return(res2)
}
#sub-functions section end---------------------------------------------------------------------

#main parallel function begin------------------------------------------------------------------------------
#alldata contains x and y, y in its first column; opt=NULL: calculate p-values not considering clinical data; numit: number of permutation iterations.
#platform is a string, used to specify the filenames of results
mpi_pvalue_fdr_manhattan=function(alldata,platform="copynumber",opt=NULL,numit=500)
{
  idx_y=which(colnames(alldata)=="platinumclass")
  #outcome
  y=alldata[,idx_y]
  y=as.factor(y)
  #gene matrix
  x=alldata[,(idx_y+1):ncol(alldata)]
  colnames(x)=gsub("-","_",colnames(x))
  mpi.bcast.Robj2slave(x)
  #columns of clinicals if applicable. 3 columns for residual disease, 1 column for age
  ncol_clinicals=2:5
  if (!is.null(opt)) #for the case considering clinical data
  {
    clinicaldata=alldata[,ncol_clinicals]
    mpi.bcast.Robj2slave(clinicaldata)
    #platform=paste0(platform,"_clinical1")
  }

  #to compute permutations
  #permutation from iteration 2:
  permutationp=NULL
  permutationp_filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,"_1000permutation.txt")

  for (i in 1:(numit+1))
  {
    y1=y
    if (i>1)
    {
      #permutation on outcome
      set.seed(i+1000)
      y1=y[sample(length(y))]
      names(y1)=names(y)
    }
    mpi.bcast.Robj2slave(y1)
    if (i==1)
    {
      #when i==1, calculate p values without permutation, and store the pvalues
      pvalue_filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,".txt")
      if (is.null(opt)) #if not considering clinical data
      {
        pvalues=mpi_calpvalues(x,y1,filename=pvalue_filename,njob=njob)
      }else
      {
        pvalues=mpi_calpvalues(data=x,y=y1,filename=pvalue_filename,opt="includeclinical",njob=njob,clinicaldata=clinicaldata)
      }
      
    }else
    {
      #when i>1, calculate p values with permutation
      if (is.null(opt))
      {
        tmp=mpi_calpvalues(x,y1,njob=njob) #if not considering clinical data
      }else
      {
        tmp=mpi_calpvalues(data=x,y=y1,opt="includeclinical",njob=njob,clinicaldata=clinicaldata)
      }
      
      #order the p values
      tmp1=tmp[order(tmp[,1]),1]
      tmp2=data.frame(pvalues=tmp1)
      rownames(tmp2)=paste0(rownames(tmp)[order(tmp[,1])],".",i-1)
      permutationp=rbind(permutationp,tmp2)
    }
    if (i %% 50==0)
    {
      cat(i,"..")
      write.table(permutationp,file=permutationp_filename)
    }
  }
  write.table(permutationp,file=permutationp_filename)
  #compute fdrs
  fdrs=mpi_compute_qvalue_permutation(pvalues=pvalues,pvalues_permutation=permutationp,platform=platform,numit=numit)
  
  #the filename for manhattan plot
  tiff(filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/manhattan_",platform,".tiff"))
  #draw_manhattan is in functions.R
  draw_manhattan(fdrs=fdrs,maxy=5,fdrthreshold=0.05)
  dev.off()
  #smallfdrtable generates the table with genes having small fdr. it is in functions.R
  smallfdr=smallfdrtable(fdrs=fdrs,data=data_copynumber,threshold=0.05)
  write.table(smallfdr,file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/smallfdr_",platform,".txt"),row.names=F,col.names=T,sep="\t",quote=F)
}
#main parallel function end------------------------------------------------------------------------------


#All aboves are functions and data
#Starting run code from here:
opt=NULL #without considering clinical data
#opt="includeclinical" #for the case wich clinical
if (is.null(opt))
{
  print("train data...")
  alldata=data_copynumber$data1
  alldata=formwholedataframe(alldata)
  mpi_pvalue_fdr_manhattan(alldata=alldata,platform="copynumber_train",opt=NULL,numit=100)
  
  print("test data...")
  alldata=data_copynumber$data2
  #formwholedataframe transform the data into a dataframe containing x and y. it is in functions.R
  alldata=formwholedataframe(alldata)
  mpi_pvalue_fdr_manhattan(alldata=alldata,platform="copynumber_test",opt=NULL,numit=100)
  
  print("all data")
  alldata=rbind(data_copynumber$data1,data_copynumber$data2)
  alldata=formwholedataframe(alldata)
  mpi_pvalue_fdr_manhattan(alldata=alldata,platform="copynumber_all",opt=NULL,numit=100)
  print("randomization...")
  #do randomization
  alldata=rbind(data_copynumber$data1,data_copynumber$data2)
  alldata=formwholedataframe(alldata)
  for (i in 1:3)
  {
    print(i)
    idxtrain=rep(T,nrow(alldata))
    set.seed(i+1000)
    idxtrain[sample(1:nrow(alldata),ceiling(nrow(alldata)/2))]=F
    alldata1=alldata[idxtrain,]
    mpi_pvalue_fdr_manhattan(alldata=alldata1,platform=paste0("copynumber_train",i),opt=NULL,numit=100)
    alldata2=alldata[!idxtrain,]
    mpi_pvalue_fdr_manhattan(alldata=alldata2,platform=paste0("copynumber_test",i),opt=NULL,numit=100)
  }
}else
{
  print("train...")
  alldata=rbind(data_copynumber$data1,data_copynumber$data2)
  numtrain=nrow(formwholedataframe(data_copynumber$data1))
  alldata=formwholedataframe(alldata)
  alldata$residual_disease_largest_nodule=as.character(alldata$residual_disease_largest_nodule)
  alldata$residual_disease_largest_nodule=gsub("11_20mm","10mm_",alldata$residual_disease_largest_nodule)
  alldata$residual_disease_largest_nodule=gsub("20mm_","10mm_",alldata$residual_disease_largest_nodule)
  idxna=is.na(alldata$residual_disease_largest_nodule)
  alldata$residual_disease_largest_nodule[idxna]="missing"
  alldata$residual_disease_largest_nodule=as.factor(alldata$residual_disease_largest_nodule)
  ncol_clinical=which(colnames(alldata)=="residual_disease_largest_nodule")
  xfactors=model.matrix(as.formula(paste0("~",colnames(alldata)[ncol_clinical])),data=alldata)[,-1]
  xfactors=cbind(sample=rownames(xfactors),xfactors)
  alldata=alldata[,-ncol_clinical]
  alldata=cbind(sample=rownames(alldata),alldata)
  alldata=merge(xfactors,alldata,by="sample")
  print("train...")
  mpi_pvalue_fdr_manhattan(alldata=alldata[1:numtrain,],platform="copynumber_clinical_train",opt=opt,numit=100)
  
  print("test...")
  mpi_pvalue_fdr_manhattan(alldata=alldata[(numtrain+1):nrow(alldata),],platform="copynumber_clinical_test",opt=opt,numit=100)
  
  print("randomization...")
  #do randomization

  for (i in 1:3)
  {
    print(i)
    idxtrain=rep(T,nrow(alldata))
    set.seed(i+1000)
    idxtrain[sample(1:nrow(alldata),ceiling(nrow(alldata)/2))]=F
    alldata1=alldata[idxtrain,]
    mpi_pvalue_fdr_manhattan(alldata=alldata1,platform=paste0("copynumber_clinical_train",i),opt=opt,numit=100)
    alldata2=alldata[!idxtrain,]
    mpi_pvalue_fdr_manhattan(alldata=alldata2,platform=paste0("copynumber_clinical_test",i),opt=opt,numit=100)
  }
}


mpi.close.Rslaves()
#program ends here
mpi.quit()