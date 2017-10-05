#!/usr/bin/env Rscript

library(glmnet)
library(MASS)
library(ROCR)
library(pROC)
#load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
selectlambda=function(x,y,numlambda=200,ntime=100,opt="binomial",penalty=NULL,standardize=TRUE,main="",nfolds=10)
{
  #ntime: number of times of cvfit
  lamb <- matrix(NA,numlambda,ntime)
  err  <- matrix(NA,numlambda,ntime)
  esd  <- matrix(NA,numlambda,ntime)
  
  if (is.null(penalty))
  {
    penalty=rep(1,ncol(x))
  }
  for (i in 1:ntime){
    cat(i,"..")
    if (opt=="binomial")
    {
      set.seed(i)
      cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=nfolds,penalty.factor=penalty,standardize=standardize)
    }
    if (opt=="numeric")
    {
      set.seed(i)
      cvfit <- cv.glmnet(as.matrix(x),y,nlambda=numlambda,nfolds=10,penalty.factor=penalty)
    }
    nn <- length(cvfit$lambda)
    lamb[1:nn,i] <- cvfit$lambda
    err[1:nn,i] <-  cvfit$cvm
    esd[1:nn,i] <-  cvfit$cvsd
  }
  
  #some small lambdas may not be used all the time
  allnarow=apply(err,1,function(x)
  {
    res=F
    if (sum(is.na(x))==length(x))
      res=T
    return(res)
  })
  lamb=lamb[!allnarow,]
  err=err[!allnarow,]
  esd=esd[!allnarow,]
  lambseq <- rep(NA,nrow(lamb))
  #mean cv error for a lambda 
  errseq <- rep(NA,nrow(lamb))
  for (i in 1:nrow(lamb)) {
    
    lambseq[i] <- min(lamb[i,!is.na(lamb[i,])])
    errseq[i] <- mean(err[i,!is.na(err[i,])])
  }
  
  idxmin=which.min(errseq)
  sdmin=sd(err[idxmin,])/sqrt(ntime-1)
  theerror=errseq[idxmin]+sdmin
  idxsel=which(errseq>=theerror)
  idxsel=idxsel[idxsel<idxmin]
  if (length(idxsel)>0)
  {
    lambda.1se=lambseq[idxsel[length(idxsel)]] #the last one
  }else
  {
    #the error curve is decreasing monotically, pick the largest lambda and no variable is selected
    lambda.1se=lambseq[1]
  }
  plot(lambseq,errseq,main=main)
  abline(v=lambda.1se,col="red")
  
  #the selected lambda
  lambda.min=lambseq[which.min(errseq)]
  abline(v=lambda.min,col="blue")
  result=list(lamb=lamb,err=err,esd=esd,lambseq=lambseq,errseq=errseq,lambda.min=lambda.min,lambda.1se=lambda.1se)
}

#plot 2 roc functions, use predicted values
plotroc2=function(pfit_train,y,pfit_test,y1,plotflag=1,main="")
{
  rocdata=function(pfit,y)
  {
    yy=as.numeric(y[!is.na(pfit)])
    pfit=pfit[!is.na(pfit)]
    xx <- cbind(pfit,yy)
    if (max(xx[,2])>=2)
      xx[,2] <- ifelse(xx[,2]>=2,1,0)
    xx <- xx[order(xx[,1]),]
    pred <- prediction(xx[,1],xx[,2])
    auc.tmp <- performance(pred,"auc")
    if (auc.tmp@y.values<0.5)
    {
      xx[,1]=1-xx[,1]
      pred <- prediction(xx[,1],xx[,2])
      auc.tmp <- performance(pred,"auc")
    }
    auc <- format(as.numeric(auc.tmp@y.values),digits=3)
    result=list(pred=pred,auc=auc)
  }
  res_train=rocdata(pfit_train,y)
  res_test=rocdata(pfit_test,y1)
  
  roc.perf_train = performance(res_train$pred, measure = "tpr", x.measure = "fpr")
  roc.perf_test = performance(res_test$pred, measure = "tpr", x.measure = "fpr")
  plot(roc.perf_train,col=3,lwd=4,main=main)
  points(unlist(roc.perf_test@x.values),unlist(roc.perf_test@y.values),col="red",type="l",lwd=4)
  text(0.5,0.8,paste0("Training AUC=",res_train$auc),cex=1)
  text(0.5,0.25,paste0("Validation AUC=",res_test$auc),cex=1)
}

#update cilinical items of "clinical stage" and "residual disease"
updateclinicalitem=function(data=data_mrna)
{
  if (is.data.frame(data)==F)
  {
    alldata=NULL
    for (i in 1:length(data))
    {
      alldata=rbind(alldata,data[[i]])
    }
  }else
  {
    alldata=data
  }
  
  #process clinical stage, merge A,B,C
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IA","Stage IB","Stage IC"))
  alldata[tmp,"clinical_stage"]="Stage1"
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IIA","Stage IIB","Stage IIC"))
  alldata[tmp,"clinical_stage"]="Stage2"
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IIIA","Stage IIIB","Stage IIIC"))
  alldata[tmp,"clinical_stage"]="Stage3"
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IV"))
  alldata[tmp,"clinical_stage"]="Stage4"
  table(alldata[,"clinical_stage"])
  #Stage1 Stage2 Stage3 Stage4 
  #8     18    300     54 
  #process residual names
  alldata[,"residual_disease_largest_nodule"]=as.character(alldata[,"residual_disease_largest_nodule"])
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "No Macroscopic disease")
  alldata[tmp,"residual_disease_largest_nodule"]="No_Macroscopic_disease"
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "1-10 mm")
  alldata[tmp,"residual_disease_largest_nodule"]="1_10mm"
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "11-20 mm")
  alldata[tmp,"residual_disease_largest_nodule"]="11_20mm"
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% ">20 mm")
  alldata[tmp,"residual_disease_largest_nodule"]="20mm_"
  alldata[,"residual_disease_largest_nodule"]=as.factor(alldata[,"residual_disease_largest_nodule"])
  return(alldata)
}

#creat data considering only from a set of genes
dataselgenes=function(data,selgenes,num_allclinical=13)
{
  if (is.null(selgenes))
  {
    data1=data
  }else
  {
    
    tmp=names(data)
    #for data structure as data1,data2,data3
    if (names(data)[1]=="data1")
    {
      idxgenes=which(colnames(data$data1) %in% selgenes)
      data1=list(data1=cbind(data$data1[,1:num_allclinical],data$data1[,idxgenes]),
                 data2=cbind(data$data2[,1:num_allclinical],data$data2[,idxgenes]))
    }else
    {
      #for randomized trainingdata and testingdata
      idxgenes=which(colnames(data$trainingdata) %in% selgenes)
      data1=list(trainingdata=cbind(data$trainingdata[,1:num_allclinical],data$trainingdata[,idxgenes]),
                 testingdata=cbind(data$testingdata[,1:num_allclinical],data$testingdata[,idxgenes]))
    }
    
  }
  return(data1)
}

classify_plat4=function(data=data_copynumber,num_clinical=4,prefix="copynumber",selgenes=NULL)
{
  #if genes are specified (as not NULL), use data only from selgenes.
  data=dataselgenes(data,selgenes)
  #update stage and residual
  alldata=updateclinicalitem(data)
  #add indicator of train/test as the last column
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(data$data1)]=TRUE
  alldata=cbind(alldata,train=train)
  
  #clinicals to remove before analysis
  remove_clinical=c("race","initial_pathologic_dx_year","vital_status","death_months_to","treatment_outcome_first_course","progression_free_survival","uselastcontact")
  alldata=alldata[,!colnames(alldata) %in% remove_clinical]
  #number of clinical variables to keep (age,grade,stage,residual)
  num_clinical=4
  
  #some genes have name with "-" character
  tmp=gsub("-","_",colnames(alldata))
  colnames(alldata)=tmp
  #to process gene names which contain "|" character
  tmp=which(grepl("?",colnames(alldata),fixed=T)==T)
  if (length(tmp)>0)
  {
    alldata=alldata[,-tmp]
  }
  tmp=sapply(colnames(alldata),function(x){
    if (grepl("|",x,fixed=T))
    {
      res=unlist(strsplit(x,"|",fixed=T))[1]
    }else
    {
      res=x
    }
    return(res)
  })
  colnames(alldata)=tmp
  
  #remove constant columns,# 2 extra variables: platinumclass and drug_interval_computed
  tmp1=sapply((num_clinical+3):ncol(alldata),function(x){
    res=FALSE
    if (var(alldata[,x])>0)
    {
      res=TRUE
    }
    return(res)
  })
  alldata=cbind(alldata[,1:(num_clinical+2)],alldata[(num_clinical+3):ncol(alldata)][,tmp1])
  
  #remove drug interval
  alldata=alldata[,-which(colnames(alldata)=="drug_interval_computed")]
  
  #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
  alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
  #only 1 sample in "G4", GX -missing data
  alldata=alldata[!alldata[,"tumor_grade"] %in% c("GX","G4"),]
  alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
  alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  
  #without clinicals
  #form the design matrix for clinical data
  xfactors=model.matrix(as.formula(paste0("~",paste0(colnames(alldata)[1:num_clinical],collapse="+"))),data=alldata)[,-1]
  #some samples may be removed due to missing data in clinical, more than 30 removed for missing in residual
  tmp=sapply(rownames(xfactors),function(x){
    res=which(rownames(alldata)==x)
  })
  alldata1=alldata[tmp,]
  #combine clinical categorical data and other numeric data
  alldata2=cbind(xfactors,alldata1[,(num_clinical+1):ncol(alldata1)])
  colnames(alldata2)=gsub(" ","_",colnames(alldata2))
  
  #the column indices of train indicator and outcome
  idx_train=which(colnames(alldata2)=="train")
  idx_y=which(colnames(alldata2)=="platinumclass")
#   
  cor_flag=F
  if (cor_flag==T)
  {
    #select top 000 genes based on correlation with output
    corrs=sapply((idx_y+1):(ncol(alldata2)-1),function(i){
      res=cor(alldata2[which(alldata2[,idx_train]==1),i],as.numeric(alldata2[which(alldata2[,idx_train]==1),idx_y]),use="complete")
      names(res)=colnames(alldata2)[i]
      return(res)
    })
    corrs_order=corrs[order(abs(corrs),decreasing=T)]
    topgenes=names(corrs_order)[1:1000]
    tmp=which(colnames(alldata2) %in% topgenes)
    alldata2=cbind(alldata2[,1:idx_y],alldata2[,tmp],train=alldata2[,idx_train])
    idx_train=which(colnames(alldata2)=="train")
    idx_y=which(colnames(alldata2)=="platinumclass")
    
  }
  
  #training data
  ally=alldata2[,idx_y]
  x=alldata2[which(alldata2[,idx_train]==1),]
  y=ally[which(alldata2[,idx_train]==1)]
  x=x[,-c(idx_y,idx_train)]
  
  #testing data
  x1=alldata2[which(alldata2[,idx_train]==0),]
  y1=ally[which(alldata2[,idx_train]==0)]
  x1=x1[,-c(idx_y,idx_train)]
  
  #standardize numeric data
  x[,idx_y:ncol(x)]=scale(x[,idx_y:ncol(x)])
  x1[,idx_y:ncol(x1)]=scale(x1[,idx_y:ncol(x1)])

  trainingsamples=rownames(x)
  testingsamples=rownames(x1)
  
  penalty=rep(1,ncol(x))
  penalty[c(7)]=0
  #trainresult=selectlambda(x,y,numlambda=200,ntime=20,opt="binomial",penalty=rep(1,ncol(x)),standardize=T,main=prefix,nfolds=10)
  trainresult=selectlambda(x,y,numlambda=100,ntime=20,opt="binomial",penalty=penalty,standardize=F,main=prefix,nfolds=10)
  #lambda.sel=trainresult$lambda.1se
  lambda.sel=trainresult$lambda.min
  #[1] [1] 0.0206985
  
  #fit=glmnet(as.matrix(x),y,family="binomial",standardize=T,nlambda = 200)
  fit=glmnet(as.matrix(x),y,family="binomial",nlambda = 100,penalty.factor=penalty)
  tmp=predict(fit,type="coefficients",s=lambda.sel)
  coeff=tmp@x
  names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  fcon=file("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_glmnetselectedgenes.txt","w")
  for (k in 2:length(coeff))
  {
    writeLines(names(coeff)[k],fcon)
  }
  close(fcon)
  save(alldata2,x,y,x1,y1,coeff,trainresult,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_glmnetclassresult.RData")
  
  pfit_train = predict(fit,as.matrix(x),s=lambda.sel,type="response")
  pfit_test = predict(fit,as.matrix(x1),s=lambda.sel,type="response")
  plotroc2(pfit_train,y,pfit_test,y1,main="ROC")
  
#   #inclue a few clinicals
#   penalty=rep(1,ncol(x))
#   penalty[which(colnames(x) %in% c("age","residual_disease_largest_noduleNo_Macroscopic_disease"))]=0
#   trainresult=selectlambda(x,y,numlambda=200,ntime=20,opt="binomial",penalty=penalty,standardize=T,main=prefix,nfolds=10)
#   lambda.sel=trainresult$lambda.min
#   fit=glmnet(as.matrix(x),y,family="binomial",standardize=T,nlambda = 200,penalty.factor=penalty)
#   tmp=predict(fit,type="coefficients",s=lambda.sel)
#   coeff=tmp@x
#   names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
#   #coeff
#   pfit_train = predict(fit,as.matrix(x),s=lambda.sel,type="response")
#   pfit_test = predict(fit,as.matrix(x1),s=lambda.sel,type="response")
#   plotroc2(pfit_train,y,pfit_test,y1,main="ROC")
  
  return(list(coeff=coeff,x=x,y=y,x1=x1,y1=y1,pfit_train=pfit_train,pfit_test=pfit_test))
}

#run the following lines to get results
#include all the genes
copynumber_res=classify_plat4(data=data_copynumber,num_clinical=4,prefix="copynumber",selgenes=NULL)

#include top 1000 genes with smallest pval
pvaluegenes_copynumber=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_train.txt",header=T)
pvaluegenes_copynumber_order=data.frame(pvalues=pvaluegenes_copynumber[order(pvaluegenes_copynumber[,1]),1])
rownames(pvaluegenes_copynumber_order)=rownames(pvaluegenes_copynumber)[order(pvaluegenes_copynumber[,1])]
copynumber_top1000_res=classify_plat4(data=data_copynumber,num_clinical=4,prefix="copynumber_top1000",selgenes=rownames(pvaluegenes_copynumber_order)[1:1000])

#include all the genes
mrna_res=classify_plat4(data=data_mrna,num_clinical=4,prefix="mrna",selgenes=NULL)

#include top 1000 genes with smallest pval
pvaluegenes_mrna=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum_train.txt",header=T)
pvaluegenes_mrna_order=data.frame(pvalues=pvaluegenes_mrna[order(pvaluegenes_mrna[,1]),1])
rownames(pvaluegenes_mrna_order)=rownames(pvaluegenes_mrna)[order(pvaluegenes_mrna[,1])]
mrna_top1000_res=classify_plat4(data=data_mrna,num_clinical=4,prefix="mrna_top1000",selgenes=rownames(pvaluegenes_mrna_order)[1:1000])



#use 1500 genes
vargenes=as.character(read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/1500genes.txt")[,1])
mrna_var1500_res=classify_plat4(data=data_mrna,num_clinical=4,prefix="mrna_1500genes",selgenes=vargenes)



