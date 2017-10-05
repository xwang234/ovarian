#!/usr/bin/env Rscript
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")
clinical_classification=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/clinical_vairables.txt",header=T,sep="\t")
train=rep(F,nrow(clinical_classification))
train[1:287]=T
clinical_classification=cbind(clinical_classification,train)

formtrain_test_data=function(data=copynumber,colnum=11,includecategories=c("1-10 mm","11-20 mm",">20 mm","No Macroscopic disease"))
{
  idxrows=clinical_classification[,colnum] %in% includecategories
  trainsamples=as.character(clinical_classification$sample[idxrows & clinical_classification$train])
  trainsamples=intersect(trainsamples,colnames(data))
  testsamples=as.character(clinical_classification$sample[idxrows & (! clinical_classification$train)])
  testsamples=intersect(testsamples,colnames(data))
  x=y=x1=y1=NULL
  for (i in 1:length(trainsamples))
  {
    x=cbind(x,data[,colnames(data)==trainsamples[i]])
    y=c(y,as.character(clinical_classification[which(clinical_classification$sample==trainsamples[i]),colnum]))
  }
  R0idx=y == "No Macroscopic disease"
  y[R0idx]=0
  y[!R0idx]=1
  rowmeanx=rowMeans(x,na.rm=T)
  for (i in 1:nrow(x))
  {
    x[i,is.na(x[i,])]=rowmeanx[i]
  }
  rownames(x)=rownames(data)
  colnames(x)=trainsamples
  for (i in 1:length(testsamples))
  {
    x1=cbind(x1,data[,colnames(data)==testsamples[i]])
    y1=c(y1,as.character(clinical_classification[which(clinical_classification$sample==testsamples[i]),colnum]))
  }
  R0idx=y1 == "No Macroscopic disease"
  y1[R0idx]=0
  y1[!R0idx]=1
  rowmeanx1=rowMeans(x1,na.rm=T)
  for (i in 1:nrow(x1))
  {
    x1[i,is.na(x1[i,])]=rowmeanx1[i]
  }
  rownames(x1)=rownames(data)
  colnames(x1)=testsamples
  result=list(x=t(x),y=y,x1=t(x1),y1=y1)
}

glmnet_classifiction=function(data=copynumber,main="copynumber")
{
  train_testdata=formtrain_test_data(data)
  trainresult=selectlambda(train_testdata$x,train_testdata$y,numlambda=100,ntime=20,opt="binomial",main=main,nfolds=10)
  #lambda.sel=trainresult$lambda.1se
  lambda.sel=trainresult$lambda.min
  
  fit=glmnet(as.matrix(train_testdata$x),train_testdata$y,family="binomial",nlambda = 100)
  tmp=predict(fit,type="coefficients",s=lambda.sel)
  coeff=tmp@x
  names(coeff)=c("intercept",colnames(train_testdata$x)[tmp@i[tmp@i>0]])
  selgenes=names(coeff)[! names(coeff) %in% c("intercept")]
  
  pfit_train = predict(fit,as.matrix(train_testdata$x),s=lambda.sel,type="response")
  pfit_test = predict(fit,as.matrix(train_testdata$x1),s=lambda.sel,type="response")
  plotroc2(pfit_train,train_testdata$y,pfit_test,train_testdata$y1,main=paste0("ROC, ",length(selgenes), " genes selected"))
  
  data1=data_copynumber_commc$data1
  data2=data_mrna_commc$data1
  if (length(selgenes)>0)
  {
    corr=sapply(selgenes,function(gene){
      idx=which(colnames(data1)==gene)
      res=NA
      names(res)=gene
      if (length(idx)>0)
      {
        res=cor(data1[,idx],data2[,idx],use="complete")
      }
      return(res)
    })
    
    idx=which(corr>0.1)
    selgenes1=selgenes[idx]
    idx=which(colnames(train_testdata$x) %in% selgenes1)
    trainresult1=selectlambda(train_testdata$x[,idx],train_testdata$y,numlambda=100,ntime=20,opt="binomial",main=main,nfolds=10)
    lambda.sel1=trainresult1$lambda.1se
    
    fit1=glmnet(as.matrix(train_testdata$x[,idx]),train_testdata$y,family="binomial",nlambda = 100)
    tmp=predict(fit1,type="coefficients",s=lambda.sel1)
    coeff1=tmp@x
    names(coeff1)=c("intercept",colnames(train_testdata$x)[tmp@i[tmp@i>0]])
    selgenes1=names(coeff1)[! names(coeff1) %in% c("intercept")]
    
    pfit_train1 = predict(fit1,as.matrix(train_testdata$x[,idx]),s=lambda.sel1,type="response")
    pfit_test1 = predict(fit1,as.matrix(train_testdata$x1[,idx]),s=lambda.sel1,type="response")
    plotroc2(pfit_train1,train_testdata$y,pfit_test1,train_testdata$y1,main=paste0("ROC, ",length(selgenes1), " genes selected"))
  }
  
}

sapply_compute_pvalue=function(data=copynumber)
{
  train_testdata=formtrain_test_data(data)
  x=rbind(train_testdata$x,train_testdata$x1)
  y=c(train_testdata$y,train_testdata$y1)
  y=factor(y)
  res=sapply(1:ncol(x),function(i){
    tmp=glm(y~x[,i],family="binomial")
    summaryres=coef(summary(tmp))
    res1=NA
    if (nrow(summaryres)>1)
    {
      res1=summaryres[2,4]
    }
    names(res1)=colnames(x)[i]
    return(res1)
  })
  return(res)
}
Sys.time()
p_copynumber=sapply_compute_pvalue(data=copynumber)
Sys.time()
p_mrna=sapply_compute_pvalue(data=mrna)
Sys.time()
hist(p_copynumber,probability = T)
hist(p_mrna,probability = T)
