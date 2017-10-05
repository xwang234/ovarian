#!/usr/bin/env Rscript

#load the training/testing data
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/mRNA_platinum_classification.RData")

library(glmnet)
library(MASS)
library(ROCR)
library(caret)
library(pROC)
library(caTools)
#select lambda
selectlambda=function(x,y,numlambda=200,ntime=100,opt="binomial",penalty=NULL,standardize=TRUE,main="")
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
      #       cvfit1=cv.glmnet(as.matrix(allx[,-c(idx_y,idx_train)]),allx[,idx_y],nlambda=numlambda,nfolds=10,penalty.factor=penalty)
      #       cvfit1=cv.glmnet(as.matrix(allx[,-c(idx_y,idx_train)]),allx[,13],family="binomial",nlambda=numlambda,nfolds=10,penalty.factor=penalty)
      #       plot(cvfit1)
      cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=10,penalty.factor=penalty,standardize=standardize)
      #cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=length(y),penalty.factor=penalty)
    }
    if (opt=="numeric")
    {
      cvfit <- cv.glmnet(as.matrix(x),y,nlambda=numlambda,nfolds=10,penalty.factor=penalty)
    }
    #cvfit <- cv.glmnet(as.matrix(x),y,nlambda=numlambda,nfolds=10)
    #cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=length(y))
    #plot(cvfit)
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
  sdmin=sd(err[idxmin,])
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


getlambda1sebyloo=function(x,y,numlambda=200,opt="binomial",opt1="nopenalty",idx_y)
{
  penalty=rep(1,ncol(x))
  if (opt1=="penalty")
  {
    penalty[1:idx_y]=0
  }
  if (opt=="binomial")
  {
    cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=length(y),penalty.factor=penalty)
  }
  if (opt=="numeric")
  {
    cvfit <- cv.glmnet(as.matrix(x),y,nlambda=numlambda,nfolds=length(y),penalty.factor=penalty)
  }
  lambda.sel=cvfit$lambda.1se
}

plotroc=function(fit1,testingdata1,opt="glm",plotflag=1,main="")
{
  testingdata1=as.data.frame(testingdata1)
  pfit=predict(fit1,newdata=testingdata1,type="response")
  if (opt %in% c("lda","qda"))
  {
    pfit=pfit$posterior[,2]
  }
  yy=testingdata1$y[!is.na(pfit)]
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
  if (plotflag==1)
  {
    roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
    plot(roc.perf,col=3,lwd=4,main=main)
    text(0.5,0.5,paste0("AUC=",auc),cex=1.5)
  }
  return(auc)
}

#use glmnet model
plotroc1=function(fit,x,y,lambda.sel=0.05,plotflag=1,main="")
{
  #   pfit = predict(fit1,newdata=x1[,selected],type="response")
  #   pfit=predict(fit1,newdata=testingdata1,type="response")
  #   ldapredict=predict(ldafit,testingdata1)
  #   pfit=ldapredict$posterior[,2]
  #   ldapredict1=predict(ldafit1,testingdata1)
  #   pfit=ldapredict1$posterior[,2]
  #   pfit=predict(glmfit,newdata=testingdata1,type="response")
  #   pfit=predict(glmfit1,type="response")
  pfit = predict(fit,as.matrix(x),s=lambda.sel,type="response")
  #   test=cbind(pfit,x[,selected])
  #   tmp=coeff
  #   tmp1=cbind(intercept=rep(1,length(pfit)),x[,selected])
  #   i=86
  #   exp(sum(tmp*tmp1[i,]))
  #   exp(sum(tmp*tmp1[i,]))/(1+exp(sum(tmp*tmp1[i,])))
  #   as.numeric(tmp)
  #   as.numeric(tmp1[i,])
  
  yy=y[!is.na(pfit)]
  pfit=pfit[!is.na(pfit)]
  xx <- cbind(pfit,yy)
  #xx <- cbind(pfit,trainingdata1$y)
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
  if (plotflag==1)
  {
    roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
    plot(roc.perf,col=3,lwd=4,main=main)
    text(0.5,0.5,paste0("AUC=",auc),cex=1.5)
  }
  return(auc)
}

plotroc2=function(fit1,testingdata1,plotflag=1,main="")
{
  testingdata1=as.data.frame(testingdata1)
  x1=testingdata1[,2:ncol(testingdata1)]
  pfit=predict(fit1,x1)
  yy=testingdata1$y[!is.na(pfit)]
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
  if (plotflag==1)
  {
    roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
    plot(roc.perf,col=3,lwd=4,main=main)
    text(0.5,0.5,paste0("AUC=",auc),cex=1.5)
  }
  return(auc)
}

plotroc3=function(pfit,y1,plotflag=1,main="")
{
  yy=as.numeric(y1[!is.na(pfit)])
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
  if (plotflag==1)
  {
    roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
    plot(roc.perf,col=3,lwd=4,main=main)
    text(0.5,0.5,paste0("AUC=",auc),cex=1.5)
  }
  return(auc)
}

classify_plat=function(data=data_u133)
{
  #   #standardize using all the data
  #   alldata=rbind(data$data1,data$data2,data$data3)
  #   alldata=alldata[,-c(1,2,3,4,5,6)]
  #   alldata=scale(alldata)
  #   data$data1[,7:ncol(data$data1)]=alldata[1:nrow(data$data1),]
  #   data$data2[,7:ncol(data$data2)]=alldata[(nrow(data$data1)+1):(nrow(data$data1)+nrow(data$data2)),]
  #   data$data3[,7:ncol(data$data3)]=alldata[(nrow(data$data1)+nrow(data$data2)+1):(nrow(alldata)),]
  
  trainingdata=data$data1
  #trainingdata=rbind(data$data1,data$data2)
  #not use the clinical variables
  x=trainingdata[,-c(1,2,3,4,5,6)]
  y=trainingdata$platinumclass
  #   y=unname(y)
  #   y=as.numeric(y)
  #   #impute missing data in x
  #   for (i in 1:ncol(x))
  #   {
  #     x[,i] <- ifelse(is.na(x[,i]),mean(x[!is.na(x[,i]),i]),x[,i])
  #   }
  trainresult=selectlambda(x,y,numlambda=200,ntime=10)
  lambda.sel=trainresult$lambda.sel
  fit=glmnet(as.matrix(x),y,family="binomial")
  plotroc1(fit,x,y,lambda.sel)
  
  #use glm model
  selected <- which(as.matrix(coef(fit,s=lambda.sel))[,1]!=0)
  #remove intercept
  if (names(selected)[1]=="(Intercept)")
  {
    selected=selected[-1]-1
  }else
  {
    selected=selected-1
  }
  
  #build a logistic regression model
  trainingdata1=cbind(y=y,trainingdata)
  fm=as.formula(paste0("y~",paste(names(selected),collapse="+")))
  glmfit=glm(fm,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit,trainingdata1)
  
  #use variables with small p values
  tmp=summary(glmfit)
  idx=which(tmp$coefficients[,4]<0.05)
  if (names(idx)[1] == "(Intercept)")
    idx=idx[-1]
  fm1=as.formula(paste0("y~",paste(names(idx),collapse="+")))
  glmfit1=glm(fm1,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit1,trainingdata1)
  #use lda
  ldafit=lda(fm,data=trainingdata1)
  plotroc(fit1=ldafit,trainingdata1,opt="lda")
  ldafit1=lda(fm1,data=trainingdata1)
  plotroc(fit1=ldafit1,trainingdata1,opt="lda")
  #use qda
  qdafit=qda(fm,data=trainingdata1)
  plotroc(fit1=qdafit,trainingdata1,opt="qda")
  qdafit1=qda(fm1,data=trainingdata1)
  plotroc(fit1=qdafit1,trainingdata1,opt="lda")
  
  #for testing data
  #testingdata=data$data2
  #testingdata=data$data3
  testingdata=rbind(data$data2,data$data3)
  x1=testingdata[,-c(1,2,3,4,5,6)]
  x1=scale(x1)
  y1=testingdata$platinumclass
  
  #   for (i in 1:ncol(x1))
  #   {
  #     x1[,i] <- ifelse(is.na(x1[,i]),mean(x1[!is.na(x[,i]),i]),x1[,i])
  #   }
  lambdas=seq(0.05,0.1,0.0025)
  for (lambda in lambdas)
  {
    plotroc1(fit,x1,y1,lambda)
  }
  plotroc1(fit,x1,y1,lambda.sel)
  testingdata1=cbind(y=y1,x1)
  plotroc(fit1=glmfit,testingdata1)
  plotroc(fit1=glmfit1,testingdata1)
  plotroc(fit1=ldafit,testingdata1,opt="lda")
  plotroc(fit1=ldafit1,testingdata1,opt="lda")
  plotroc(fit1=qdafit,testingdata1,opt="qda")
  plotroc(fit1=qdafit1,testingdata1,opt="qda")
}



interpl_factor=function(factorvector)
{
  tmp=table(factorvector)
  tmp1=tmp[order(tmp,decreasing=T)]
  tmp2=which(is.na(factorvector))
  if (length(tmp2)>0)
  {
    factorvector[tmp2]=names(tmp1)[1]
  }
  return(factorvector)
}

#include clinical variables
classify_plat1=function(data=data_u133)
{
  #for predicted copy number data+ residual copy number
  tmp=ncol(data_residualscopynumber$data1)
  colnames(data_residualscopynumber$data1)[9:tmp]=colnames(data_residualscopynumber$data2)[9:tmp]=
    colnames(data_residualscopynumber$data3)[9:tmp]=paste0("res_",colnames(data_residualscopynumber$data1[9:tmp]))
  data=data_fittedcopynumber
  data$data1=cbind(data$data1,data_residualscopynumber$data1[,9:tmp])
  data$data2=cbind(data$data2,data_residualscopynumber$data2[,9:tmp])
  data$data3=cbind(data$data3,data_residualscopynumber$data3[,9:tmp])
  
  #number of clinical variables included
  num_clinical=5
  alldata=rbind(data$data1,data$data2,data$data3)
  #remove "uselastcontact"
  tmp=which(colnames(alldata)=="uselastcontact")
  alldata=alldata[,-tmp] 
  
  #some genes have name with "-" character
  tmp=gsub("-","_",colnames(alldata))
  colnames(alldata)=tmp
  #to process colnames in methylation data which contain |
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
  
  #remove constant columns
  tmp1=sapply((num_clinical+3):ncol(alldata),function(x){
    res=FALSE
    if (var(alldata[,x])>0)
    {
      res=TRUE
    }
    return(res)
  })
  
  alldata=cbind(alldata[,1:(num_clinical+2)],alldata[(num_clinical+3):ncol(alldata)][,tmp1])
  #scale those numeric data
  #alldata[,c(1,(num_clinical+3):ncol(alldata))]=scale(alldata[,c(1,(num_clinical+3):ncol(alldata))])
  #alldata[,c((num_clinical+3):ncol(alldata))]=scale(alldata[,c((num_clinical+3):ncol(alldata))])
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(data$data1)]=TRUE
  alldata=cbind(alldata,train=train)
  #mannualy correct those mismatches
  mannualclass=list(samples=c("TCGA-13-0803","TCGA-29-1774","TCGA-24-2029","TCGA-61-2113,TCGA-61-2009"),class=c("Sensitive","Resistant","Resistant","Resistant","Resistant"))
  for (i in 1:length(mannualclass[[1]]))
  {
    tmp=mannualclass[[1]][i]
    tmp1=which(rownames(alldata) %in% tmp)
    tmp2=mannualclass[[2]][i]
    alldata$platinumclass[tmp1]=tmp2
  }
  #remove samples having computed interval close to the boundary
  tmp=which(alldata$drug_interval_computed>5.9 & alldata$drug_interval_computed<6.1)
  if (length(tmp)>0)
  {
    alldata=alldata[-tmp,]
  }
  
  #remove drug interval
  computed_interval=alldata[,num_clinical+2]
  alldata=alldata[,-(num_clinical+2)]
  #   #check clinical data
  #   table(alldata[,2],useNA="ifany")
  #   #AMERICAN INDIAN OR ALASKA NATIVE                            ASIAN        BLACK OR AFRICAN AMERICAN                            WHITE 
  #   #2                               14                               17                              330 
  # #  <NA> 
  # #    15
  #   table(alldata[,3],useNA="ifany")
  #   #G2   G3   G4   GX   G1 <NA> 
  #   #46  321    1    7    2    1 
  #   table(alldata[,4],useNA="ifany")
  #   #Stage IIA  Stage IIC Stage IIIA Stage IIIB Stage IIIC   Stage IV  Stage IIB   Stage IA   Stage IB   Stage IC       <NA> 
  # #  2         15          6         16        275         54          1          1          1          6          1 
  #   table(alldata[,5],useNA="ifany")
  #   #>20 mm                1-10 mm               11-20 mm No Macroscopic disease                   <NA> 
  #   #  75                    164                     21                     85                     33 
  #   mod1=as.formula(paste0("platinumclass~",paste0(colnames(alldata)[1:(num_clinical)],collapse="+")))
  #   test=glm(mod1,family="binomial",data=alldata)
  #   summary(test)
  # #   Coefficients:
  # #     Estimate Std. Error z value Pr(>|z|)  
  # #   (Intercept)                                             34.34835 4755.23646   0.007   0.9942  
  # #   age                                                     -0.02309    0.01204  -1.919   0.0550 .
  # #   raceASIAN                                              -14.97427 2638.35395  -0.006   0.9955  
  # #   raceBLACK OR AFRICAN AMERICAN                          -16.10638 2638.35380  -0.006   0.9951  
  # #   raceWHITE                                              -16.29546 2638.35375  -0.006   0.9951  
  # #   tumor_gradeG3                                           -0.24777    0.41002  -0.604   0.5456  
  # #   tumor_gradeG4                                          -18.40400 3956.18035  -0.005   0.9963  
  # #   tumor_gradeGX                                           16.67351 1703.21012   0.010   0.9922  
  # #   tumor_gradeG1                                           16.04374 2786.00784   0.006   0.9954  
  # #   clinical_stageStage IIC                                -14.50047 3956.18048  -0.004   0.9971  
  # #   clinical_stageStage IIIA                                 0.54970 4298.82242   0.000   0.9999  
  # #   clinical_stageStage IIIB                               -15.07398 3956.18042  -0.004   0.9970  
  # #   clinical_stageStage IIIC                               -15.76184 3956.18035  -0.004   0.9968  
  # #   clinical_stageStage IV                                 -15.80852 3956.18036  -0.004   0.9968  
  # #   clinical_stageStage IIB                                  0.08613 5594.88388   0.000   1.0000  
  # #   clinical_stageStage IA                                  -0.14479 5594.88388   0.000   1.0000  
  # #   clinical_stageStage IC                                   0.08966 4415.69452   0.000   1.0000  
  # #   residual_disease_largest_nodule1-10 mm                  -0.27543    0.32306  -0.853   0.3939  
  # #   residual_disease_largest_nodule11-20 mm                 -0.21999    0.55714  -0.395   0.6930  
  # #   residual_disease_largest_noduleNo Macroscopic disease    1.01414    0.45401   2.234   0.0255 *
  #race is not related, and we have 15 NAs, so remove it
  alldata=alldata[,colnames(alldata) !="race"]
  num_clinical=4
  #there are few patines in gradeG1,combine G1 and G2
  alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
  tmp=which(alldata[,"tumor_grade"] %in% c("G1","G2"))
  alldata[tmp,"tumor_grade"]="G1or2"
  #there are few patines in gradeG4,combine G3 and G4
  tmp=which(alldata[,"tumor_grade"] %in% c("G3","G4"))
  alldata[tmp,"tumor_grade"]="G3or4"
  #GX: Grade cannot be assessed (undetermined grade)
  tmp=which(alldata[,"tumor_grade"] %in% c("GX"))
  alldata[tmp,"tumor_grade"]=NA
  alldata[,"tumor_grade"]=as.factor(alldata[,"tumor_grade"])
  
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
  #8     18    297     54 
  #there are not many patients in sage1, combine stages1 and 2
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage1","Stage2"))
  alldata[tmp,"clinical_stage"]="Stage1or2"
  alldata[,"clinical_stage"]=as.factor(alldata[,"clinical_stage"])
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
  #interplate categorical clinical data
  interplate_flag=F
  if (interplate_flag==T)
  {
    alldata[,2]=interpl_factor(alldata[,2])
    alldata[,3]=interpl_factor(alldata[,3])
    alldata[,4]=interpl_factor(alldata[,4])
    
  }
  
  #form the design matrix for clinical data
  xfactors=model.matrix(as.formula(paste0("~",paste0(colnames(alldata)[1:num_clinical],collapse="+"))),data=alldata)[,-1]
  #some samples may be removed due to missing data in clinical
  tmp=sapply(rownames(xfactors),function(x){
    res=which(rownames(alldata)==x)
  })
  alldata1=alldata[tmp,]
  #combine clinical categorical data and other numeric data
  alldata2=cbind(xfactors,alldata1[,(num_clinical+1):ncol(alldata1)])
  colnames(alldata2)=gsub(" ","_",colnames(alldata2))
  
  idx_train=which(colnames(alldata2)=="train")
  idx_y=which(colnames(alldata2)=="platinumclass")
  #multiple genes may have the save data (in copynumber),just keep one of them
  correlationMatrix <- cor(alldata2[,-c(idx_y,idx_train)])
  highlyCorrelated <- findCorrelation(correlationMatrix, names=T,cutoff=0.9999)
  alldata2=alldata2[,-which(colnames(alldata2) %in% highlyCorrelated)]
  rm(correlationMatrix)
  idx_train=which(colnames(alldata2)=="train")
  #   mod1=as.formula(paste0("platinumclass~",paste0(colnames(alldata2)[1:ncol(xfactors)],collapse="+")))
  #   test=glm(mod1,family="binomial",data=alldata2[1:250,1:15])
  #   summary(test)
  #form training/testing data
  
  ally=alldata2[,idx_y]
  x=alldata2[which(alldata2[,idx_train]==1),]
  y=ally[which(alldata2[,idx_train]==1)]
  x=x[,-c(idx_y,idx_train)]
  trainingsamples=rownames(x)
  #testing data
  x1=alldata2[which(alldata2[,idx_train]==0),]
  y1=ally[which(alldata2[,idx_train]==0)]
  x1=x1[,-c(idx_y,idx_train)]
  testingsamples=rownames(x1)
  
  #cmopute correlation between mrna and class, filter genes using correlation first
  corrs=sapply(1:ncol(x),function(i){
    res=cor(x[,i],as.numeric(y),use="complete")
  })
  #1 cor is 0, due to constant 
  corrs[is.na(corrs)]=0
  names(corrs)=colnames(x)
  corrs1=corrs[order(abs(corrs),decreasing=T)]
  #include clinical data
  selvars=unique(c(colnames(x)[1:(idx_y-1)],names(corrs1)[1:1000]))
  #selvars=unique(c(colnames(x)[1:(idx_y-1)],names(corrs1)[1:2000]))
  length(selvars)
  tmp=sapply(selvars,function(i){
    res=which(colnames(x)==i)
  })
  x=x[,tmp]
  #standardize numeric data
  x[,idx_y:ncol(x)]=scale(x[,idx_y:ncol(x)])
  x1=x1[,tmp]
  x1[,idx_y:ncol(x1)]=scale(x1[,idx_y:ncol(x1)])
  #plot(colMeans(x))
  #   points(colMeans(x1),col="red")
  
  #train
  trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",opt1="penalty",idx_y=idx_y-1,standardize=FALSE)
  lambda.sel=trainresult$lambda.1se
  #lambda.sel=trainresult$lambda.min
  #lambda.sel=getlambda1sebyloo(x,y,numlambda=200,opt="binomial",opt1="penalty",idx_y)
  
  penalty=rep(1,ncol(x))
  penalty[1:(idx_y-1)]=0
  
  #not standardize categorical data
  #fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty)
  fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty,standardize=FALSE)
  tmp=predict(fit,type="coefficients",s=lambda.sel)
  coeff=tmp@x
  names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  coeff
  plotroc1(fit,x,y,lambda.sel,main="glmnet,train")
  plotroc1(fit,x1,y1,lambda.sel,main="glmnet,test") #0.733
  
  #   fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty,pmax=45,standardize=FALSE)
  #   tmp=fit$lambda
  #   lambda.sel=tmp[length(tmp)]
  #   plotroc1(fit,x,y,lambda.sel,main="glmnet,train")
  #   plotroc1(fit,x1,y1,lambda.sel,main="glmnet,test")
  
  
  #use glm model
  selected <- which(as.matrix(coef(fit,s=lambda.sel))[,1]!=0)
  #remove intercept
  if (names(selected)[1]=="(Intercept)")
  {
    selected=selected[-1]-1
  }else
  {
    selected=selected-1
  }
  
  #   #use pca and glmnet,not work
  #   xsel=x[,selected]
  #   x1sel=x1[,selected]
  #   x1sel1=x1sel[,idx_y:ncol(x1sel)]
  #   xsel1=xsel[,idx_y:ncol(xsel)]
  #   tmp=prcomp(xsel1)
  #   pcax=tmp$x
  #   pcax1=as.matrix(x1sel1) %*% tmp$rotation
  #   pcax=cbind(x[,1:(idx_y-1)],pcax)
  #   pcax1=cbind(x1[,1:(idx_y-1)],pcax1)
  #   penalty=rep(1,ncol(x))
  #   penalty[1:(idx_y-1)]=0
  #   
  #   trainresult=selectlambda(pcax,y,numlambda=200,ntime=10,opt="binomial",opt1="penalty",idx_y=idx_y-1,standardize=F)
  #   lambda.sel=trainresult$lambda.1se
  #   fit=glmnet(as.matrix(pcax),y,family="binomial",penalty.factor=penalty,standardize=T)
  #   selected1 <- which(as.matrix(coef(fit,s=lambda.sel))[,1]!=0)
  #   #remove intercept
  #   if (names(selected1)[1]=="(Intercept)")
  #   {
  #     selected1=selected1[-1]-1
  #   }else
  #   {
  #     selected1=selected1-1
  #   }
  #   plotroc1(fit,pcax,y,lambda.sel,main="glmnet,train,pca")
  #   plotroc1(fit,pcax1,y1,lambda.sel,main="glmnet,test,pca") #0.71
  #   
  #   
  
  #mrna_exon_selected=names(selected)[(idx_y+1):length(selected)]
  #mrna_u133_selected=names(selected)[(idx_y+1):length(selected)]
  #mrna_4502a_selected=names(selected)[(idx_y+1):length(selected)]
  #copynumber_snp6_selected=names(selected)[(idx_y+1):length(selected)]
  #copynumber_CGH1M_selected=names(selected)[(idx_y+1):length(selected)]
  #methylation_selected=names(selected)[(idx_y+1):length(selected)]
  
  
  #build a logistic regression model
  trainingdata1=as.data.frame(cbind(y=y,x))
  trainingdata1[,1]=as.factor(trainingdata1[,1])
  fm=as.formula(paste0("y~",paste(names(selected),collapse="+")))
  glmfit=glm(fm,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit,trainingdata1,main="glmfit,train")
  #only include clinicals
  fm_clinical=as.formula(paste0("y~",paste0(colnames(trainingdata1)[2:idx_y],collapse="+")))
  #fm_clinical=as.formula(paste0("y~",paste(names(selected)[1:11],collapse="+")))
  glmfit_clinical=glm(fm_clinical,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit_clinical,trainingdata1,main="glmfit,clinical,train")
  
  testingdata1=as.data.frame(cbind(y=y1,x1))
  plotroc(fit1=glmfit,testingdata1,main="glmfit,test")
  plotroc(fit1=glmfit_clinical,testingdata1,main="glmfit,clinical,test")
  
  library(pls)
  trainingdata1=as.data.frame(cbind(y=as.numeric(y)-1,x))
  pcr.fit=pcr(fm,data=trainingdata1,validation="CV")
  summary(pcr.fit)
  validationplot(pcr.fit,val.type="MSEP")
  plotroc2(fit1=pcr.fit,testingdata1,ncomp=20,plotflag=1,main="")
  #   #use variables with small p values
  #   tmp=summary(glmfit)
  #   idx=which(tmp$coefficients[,4]<0.05)
  #   if (names(idx)[1] == "(Intercept)")
  #     idx=idx[-1]
  #   fm1=as.formula(paste0("y~",paste(names(idx),collapse="+")))
  #   glmfit1=glm(fm1,data=trainingdata1,family="binomial")
  #   plotroc(fit1=glmfit1,trainingdata1)
  #use lda
  ldafit=lda(fm,data=trainingdata1)
  plotroc(fit1=ldafit,trainingdata1,opt="lda",main="lda,train")
  lda_clinical_fit=lda(fm_clinical,data=trainingdata1)
  plotroc(fit1=lda_clinical_fit,trainingdata1,opt="lda",main="lda,clinical,train")
  #   ldafit1=lda(fm1,data=trainingdata1)
  #   plotroc(fit1=ldafit1,trainingdata1,opt="lda")
  #use qda
  #   qdafit=qda(fm,data=trainingdata1)
  #   plotroc(fit1=qdafit,trainingdata1,opt="qda")
  #   qdafit1=qda(fm1,data=trainingdata1)
  #   plotroc(fit1=qdafit1,trainingdata1,opt="lda")
  
  #for testing data
  #   for (i in 1:ncol(x1))
  #   {
  #     x1[,i] <- ifelse(is.na(x1[,i]),mean(x1[!is.na(x[,i]),i]),x1[,i])
  #   }
  
  plotroc1(fit,x1,y1,lambda.sel,main="glmnet,test")
  testingdata1=as.data.frame(cbind(y=y1,x1))
  testingdata1[,1]=as.factor(testingdata1[,1])
  plotroc(fit1=glmfit,testingdata1,main="glmfit,test")
  #   #variables selected by training, parameters were refined by testing data
  #   glmfit_test=glm(fm,data=testingdata1,family="binomial")
  #   plotroc(fit1=glmfit_test,testingdata1)
  #check if the 
  plotroc(fit1=glmfit_clinical,testingdata1,main="glmfit,clinical,test")
  glmfit_clinical1=glm(fm_clinical,data=testingdata1,family="binomial")
  plotroc(fit1=glmfit_clinical1,testingdata1,main="glmfit,clinical,test1")
  
  # plotroc(fit1=glmfit1,testingdata1)
  plotroc(fit1=ldafit,testingdata1,opt="lda")
  plotroc(fit1=lda_clinical_fit,testingdata1,opt="lda",main="lda,clinical,test")
  # plotroc(fit1=ldafit1,testingdata1,opt="lda")
  #   plotroc(fit1=qdafit,testingdata1,opt="qda")
  #   plotroc(fit1=qdafit1,testingdata1,opt="qda")
  
  #use other feature selection method
  # define the control using a random forest selection function
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)
  control <- rfeControl(functions=rfFuncs, method="boot", number=10)
  results <- rfe(x, y, sizes=c(1:50), rfeControl=control)
  results <- rfe(x, y, sizes=c(1:10), rfeControl=control)
  plot(results)
  tmp=predictors(results)
  plotroc2(fit1=results$fit,trainingdata1,main="rf,train")
  plotroc2(fit1=results$fit,testingdata1,main="rf,test")
  
  set.seed(1)
  control=trainControl(method = "repeatedCV", number = 10, repeats = 5, 
                       returnResamp = "all", classProbs = TRUE, summaryFunction = twoClassSummary)
  model <- train(x, y, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(.alpha =1, .lambda = seq(0.01, 0.05, by = 0.005)), trControl = control)
  
  plot(model, metric = "ROC")
  tmp <- predict(model, newdata = x, type = "prob")
  plotroc3(tmp$Sensitive,y,main="caret,glmnet,train")
  tmp <- predict(model, newdata = x1, type = "prob")
  #colAUC(test, y1)
  plotroc3(tmp$Sensitive,y1,main="caret,glmnet,test") #0.732
  
  RFE <- rfe(x, y, sizes = seq(30, 70, by = 3), metric = "ROC", maximize = TRUE, 
             rfeControl = control, method = "glmnet", tuneGrid = expand.grid(.alpha = 1, .lambda = c(0.01, 0.02)), trControl = control)
  
  NewVars <- RFE$optVariables
  plot(RFE)
  
  #use pls
  plsFit <- train(y~., data=trainingdata1[,selected], method = "gpls", tuneLength = 20,trControl = control,metric = "ROC")
  gpfit <- predict(plsFit, newdata = trainingdata1[,selected],type = "prob")
  plotroc3(pfit$Sensitive,y,main="caret,pls,glmnet,train")
  pfit <- predict(plsFit, newdata = testingdata1[,selected],type = "prob")
  plotroc3(pfit$Sensitive,y1,main="caret,pls,glmnet,test") #0.68
  pfit <- predict(plsFit, newdata = testingdata1[,selected])
  confusionMatrix(pfit,y1)
  
}

#use specific training/testing samples and platform
classify_plat2=function(data=data_u133,trainingsamples,testingsamples)
{
  #standardize using all the data
  alldata=rbind(data$data1,data$data2,data$data3)
  #some genes have name with "-" character
  tmp=gsub("-","_",colnames(alldata))
  colnames(alldata)=tmp
  #scale those numeric data
  alldata[,c(1,8:ncol(alldata))]=scale(alldata[,c(1,8:ncol(alldata))])
  #   train=rep(FALSE,nrow(alldata))
  #   train[1:nrow(data$data1)]=TRUE
  #   alldata=cbind(alldata,train=train)
  #remove drug interval
  alldata=alldata[,-6]
  #process clinical stage, merge A,B,C
  alldata[,3]=as.character(alldata[,3])
  tmp=which(alldata[,3] %in% c("Stage IA","Stage IB","Stage IC"))
  alldata[tmp,3]="Stage1"
  tmp=which(alldata[,3] %in% c("Stage IIA","Stage IIB","Stage IIC"))
  alldata[tmp,3]="Stage2"
  tmp=which(alldata[,3] %in% c("Stage IIIA","Stage IIIB","Stage IIIC"))
  alldata[tmp,3]="Stage3"
  tmp=which(alldata[,3] %in% c("Stage IV"))
  alldata[tmp,3]="Stage4"
  alldata[,3]=as.factor(alldata[,3])
  #process residual names
  alldata[,4]=as.character(alldata[,4])
  tmp=which(alldata[,4] %in% "No Macroscopic disease")
  alldata[tmp,4]="No_Macroscopic_disease"
  tmp=which(alldata[,4] %in% "1-10 mm")
  alldata[tmp,4]="1_10mm"
  tmp=which(alldata[,4] %in% "11-20 mm")
  alldata[tmp,4]="11_20mm"
  tmp=which(alldata[,4] %in% ">20 mm")
  alldata[tmp,4]="20mm_"
  alldata[,4]=as.factor(alldata[,4])
  #interplate categorical clinical data
  alldata[,2]=interpl_factor(alldata[,2])
  alldata[,3]=interpl_factor(alldata[,3])
  alldata[,4]=interpl_factor(alldata[,4])
  
  #form the design matrix for clinical data
  tmp=model.matrix(~.-1,data=alldata[,1:6])
  #remove variables with low variation
  tmp1=sapply(2:ncol(tmp),function(x){
    res=FALSE
    tmp2=min(table(tmp[,x]))
    if (tmp2<nrow(tmp)/100)
    {
      res=TRUE
    }
    return(res)
  })
  tmp=cbind(age=tmp[,1],tmp[,2:ncol(tmp)][,!tmp1])
  
  #lose of samples due to missing values
  tmp1=sapply(rownames(tmp),function(x){
    res=which(rownames(alldata)==x)
  })
  alldata=alldata[tmp1,]
  allx=cbind(tmp,alldata[,7:ncol(alldata)])
  #allx=model.matrix(~.-1,data=alldata)
  #idx_train=which(colnames(allx)=="train")
  #allx[,idx_train]=as.logical(allx[,idx_train])
  idx_y=which(colnames(allx)=="platinumclassSensitive")
  #training data
  
  x=allx[which(rownames(allx) %in% trainingsamples),]
  y=x[,idx_y]
  x=x[,-c(idx_y)]
  #testing data
  x1=allx[which(rownames(allx) %in% testingsamples),]
  y1=x1[,idx_y]
  x1=x1[,-c(idx_y)]
  
  #cmopute correlation between mrna and class, filter genes using correlation first
  corrs=sapply(1:ncol(x),function(i){
    res=cor(x[,i],y,use="complete")
  })
  #1 cor is 0, due to constant 
  corrs[is.na(corrs)]=0
  names(corrs)=colnames(x)
  corrs1=corrs[order(abs(corrs),decreasing=T)]
  #include clinical data
  selvars=unique(c(colnames(x)[1:idx_y],names(corrs1)[1:1000]))
  length(selvars)
  tmp=sapply(selvars,function(i){
    res=which(colnames(x)==i)
  })
  x=x[,tmp]
  x[,(idx_y+1):ncol(x)]=scale(x[,(idx_y+1):ncol(x)])
  x1=x1[,tmp]
  x1[,(idx_y+1):ncol(x1)]=scale(x1[,(idx_y+1):ncol(x1)])
  
  plot(colMeans(x))
  points(colMeans(x1),col="red")
  
  trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",opt1="penalty",idx_y)
  lambda.sel=trainresult$lambda.sel
  
  penalty=rep(1,ncol(x))
  penalty[1:idx_y]=0
  
  fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty)
  plotroc1(fit,x,y,lambda.sel)
  
  #use glm model
  selected <- which(as.matrix(coef(fit,s=lambda.sel))[,1]!=0)
  #remove intercept
  if (names(selected)[1]=="(Intercept)")
  {
    selected=selected[-1]-1
  }else
  {
    selected=selected-1
  }
  
  #build a logistic regression model
  trainingdata1=as.data.frame(cbind(y=y,x))
  fm=as.formula(paste0("y~",paste(names(selected),collapse="+")))
  glmfit=glm(fm,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit,trainingdata1)
  #only include clinicals
  fm_clinical=as.formula(paste0("y~",paste(names(selected)[1:12],collapse="+")))
  glmfit_clinical=glm(fm_clinical,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit_clinical,trainingdata1)
  
  #use variables with small p values
  tmp=summary(glmfit)
  idx=which(tmp$coefficients[,4]<0.05)
  if (names(idx)[1] == "(Intercept)")
    idx=idx[-1]
  fm1=as.formula(paste0("y~",paste(names(idx),collapse="+")))
  glmfit1=glm(fm1,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit1,trainingdata1)
  #use lda
  ldafit=lda(fm,data=trainingdata1)
  plotroc(fit1=ldafit,trainingdata1,opt="lda")
  ldafit1=lda(fm1,data=trainingdata1)
  plotroc(fit1=ldafit1,trainingdata1,opt="lda")
  #use qda
  #   qdafit=qda(fm,data=trainingdata1)
  #   plotroc(fit1=qdafit,trainingdata1,opt="qda")
  #   qdafit1=qda(fm1,data=trainingdata1)
  #   plotroc(fit1=qdafit1,trainingdata1,opt="lda")
  
  #for testing data
  #   for (i in 1:ncol(x1))
  #   {
  #     x1[,i] <- ifelse(is.na(x1[,i]),mean(x1[!is.na(x[,i]),i]),x1[,i])
  #   }
  
  plotroc1(fit,x1,y1,lambda.sel)
  testingdata1=as.data.frame(cbind(y=y1,x1))
  plotroc(fit1=glmfit,testingdata1)
  #variables selected by training, parameters were refined by testing data
  glmfit_test=glm(fm,data=testingdata1,family="binomial")
  plotroc(fit1=glmfit_test,testingdata1)
  #check if the 
  plotroc(fit1=glmfit_clinical,testingdata1)
  plotroc(fit1=glmfit1,testingdata1)
  plotroc(fit1=ldafit,testingdata1,opt="lda")
  plotroc(fit1=ldafit1,testingdata1,opt="lda")
  #   plotroc(fit1=qdafit,testingdata1,opt="qda")
  #   plotroc(fit1=qdafit1,testingdata1,opt="qda")
}
#update cilinical items
updateclinicalitem=function(data=data_exon)
{
  if (is.list(data)==T)
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
dataselgenes=function(data,selgenes,num_allclinical=12)
{
  if (is.null(selgenes))
  {
    data1=data
  }else
  {
    
    tmp=names(data)
    if (names(data)[1]=="data1")
    {
      idxgenes=which(colnames(data$data1) %in% selgenes)
      data1=list(data1=cbind(data$data1[,1:num_allclinical],data$data1[,idxgenes]),
                 data2=cbind(data$data2[,1:num_allclinical],data$data2[,idxgenes]),
                 data3=cbind(data$data3[,1:num_allclinical],data$data3[,idxgenes]))
    }else
    {
      idxgenes=which(colnames(data$trainingdata) %in% selgenes)
      data1=list(trainingdata=cbind(data$trainingdata[,1:num_allclinical],data$trainingdata[,idxgenes]),
                 testingdata=cbind(data$testingdata[,1:num_allclinical],data$testingdata[,idxgenes]))
    }
    
  }
  return(data1)
}
classify_plat3=function(data=data_copynumber_commc,useclinical=T,num_clinical=4,prefix="copynumber,useclinical,",selgenes=NULL)
{
  num_clinical=4
  data=dataselgenes(data,selgenes)
  alldata=updateclinicalitem(data)
  
  #clinicals to remove before analysis
  remove_clinical=c("race","initial_pathologic_dx_year","vital_status","death_months_to","progression_free_survival","uselastcontact")
  alldata=alldata[,!colnames(alldata) %in% remove_clinical]
  
  #some genes have name with "-" character
  tmp=gsub("-","_",colnames(alldata))
  colnames(alldata)=tmp
  #to process colnames in methylation data which contain |
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
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(data$data1)]=TRUE
  alldata=cbind(alldata,train=train)
  usecorrection=F
  if (usecorrection==T)
  {
    #mannualy correct those mismatches
    mannualclass=list(samples=c("TCGA-13-0803","TCGA-29-1774","TCGA-24-2029","TCGA-61-2113,TCGA-61-2009"),class=c("Sensitive","Resistant","Resistant","Resistant","Resistant"))
    for (i in 1:length(mannualclass[[1]]))
    {
      tmp=mannualclass[[1]][i]
      tmp1=which(rownames(alldata) %in% tmp)
      tmp2=mannualclass[[2]][i]
      alldata$platinumclass[tmp1]=tmp2
    }
    #remove samples having computed interval close to the boundary
    tmp=which(alldata$drug_interval_computed>5.9 & alldata$drug_interval_computed<6.1)
    if (length(tmp)>0)
    {
      alldata=alldata[-tmp,]
    }
  }
  
  #remove drug interval
  computed_interval=alldata[,"drug_interval_computed"]
  alldata=alldata[,-which(colnames(alldata)=="drug_interval_computed")]
  
  #interplate categorical clinical data
  interplate_flag=F
  if (interplate_flag==T)
  {
    alldata[,2]=interpl_factor(alldata[,2])
    alldata[,3]=interpl_factor(alldata[,3])
    alldata[,4]=interpl_factor(alldata[,4])
    
  }
  #remove grade 1(G1),grad4 (G4) in tumor_grade and Stage1 in clinical_stage
  alldata=alldata[alldata[,"tumor_grade"]!="G1",]
  alldata=alldata[alldata[,"tumor_grade"]!="G4",]
  alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
  alldata=alldata[alldata[,"clinical_stage"]!="Stage1",]
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  
  if (useclinical==T)
  {
    #form the design matrix for clinical data
    xfactors=model.matrix(as.formula(paste0("~",paste0(colnames(alldata)[1:num_clinical],collapse="+"))),data=alldata)[,-1]
    #some samples may be removed due to missing data in clinical
    tmp=sapply(rownames(xfactors),function(x){
      res=which(rownames(alldata)==x)
    })
    alldata1=alldata[tmp,]
    #combine clinical categorical data and other numeric data
    alldata2=cbind(xfactors,alldata1[,(num_clinical+1):ncol(alldata1)])
  }else
  {
    alldata2=alldata[,(num_clinical+1):ncol(alldata)]
  }
  colnames(alldata2)=gsub(" ","_",colnames(alldata2))
  
  idx_train=which(colnames(alldata2)=="train")
  idx_y=which(colnames(alldata2)=="platinumclass")
  #check clinical variables
  test=glm(as.formula(paste0("platinumclass~",paste0(colnames(alldata2)[1:(idx_y-1)],collapse="+"))),data=alldata2[1:255,],family="binomial")
  #summary(test)
  #     #multiple genes may have the save data (in copynumber),just keep one of them
  #     correlationMatrix <- cor(alldata2[,-c(idx_y,idx_train)])
  #     highlyCorrelated <- findCorrelation(correlationMatrix, names=T,cutoff=0.999)
  #     if (length(highlyCorrelated)>0)
  #     {
  #       alldata2=alldata2[,-which(colnames(alldata2) %in% highlyCorrelated)]
  #     }
  #     rm(correlationMatrix)
  #     idx_train=which(colnames(alldata2)=="train")
  #   mod1=as.formula(paste0("platinumclass~",paste0(colnames(alldata2)[1:ncol(xfactors)],collapse="+")))
  #   test=glm(mod1,family="binomial",data=alldata2[1:250,1:15])
  #   summary(test)
  #form training/testing data
  
  ally=alldata2[,idx_y]
  x=alldata2[which(alldata2[,idx_train]==1),]
  y=ally[which(alldata2[,idx_train]==1)]
  x=x[,-c(idx_y,idx_train)]
  
  #testing data
  x1=alldata2[which(alldata2[,idx_train]==0),]
  y1=ally[which(alldata2[,idx_train]==0)]
  x1=x1[,-c(idx_y,idx_train)]
  
  #save(x,x1,y,y1,file="../data/test_mrna.RData")
  #save(x,x1,y,y1,file="../data/test_copynumber.RData")
  
  #standardize numeric data
  x[,idx_y:ncol(x)]=scale(x[,idx_y:ncol(x)])
  x1[,idx_y:ncol(x1)]=scale(x1[,idx_y:ncol(x1)])
  dim(x)
  trainingsamples=rownames(x)
  testingsamples=rownames(x1)
  
  #train
  if (useclinical==T)
  {
    penalty=rep(1,ncol(x))
    penalty[1:(idx_y-1)]=0
    trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=T,main=prefix)
    #     penalty=rep(1,ncol(x))
    #     penalty[(idx_y-1)]=0
    #     trainresult6=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=T) #only include x[,7]
  }else
  {
    trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=TRUE,main=prefix)
  }
  if (trainresult$lambda.1se==trainresult$lambseq[1])
  {
    print("use min lambda")
    lambda.sel=trainresult$lambda.min
  }else
  {
    print("use 1se lambda")
    lambda.sel=trainresult$lambda.1se
  }
  #lambda.sel=getlambda1sebyloo(x,y,numlambda=200,opt="binomial",opt1="penalty",idx_y)
  
  fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty,standardize=T,nlambda = 200)
  tmp=predict(fit,type="coefficients",s=lambda.sel)
  coeff=tmp@x
  names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  #coeff
  plotroc1(fit,x,y,lambda.sel,main=paste0(prefix," glmnet,train"))
  plotroc1(fit,x1,y1,lambda.sel,main=paste0(prefix, "glmnet,test")) #0.733
  
  trainresult1=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=rep(1,ncol(x)),standardize=T,main=prefix)
  if (trainresult1$lambda.1se==trainresult1$lambseq[1])
  {
    print("use min lambda")
    lambda.sel1=trainresult1$lambda.min
  }else
  {
    print("use 1se lambda")
    lambda.sel1=trainresult1$lambda.1se
  }
  fit1=glmnet(as.matrix(x),y,family="binomial",standardize=T,nlambda = 200)
  tmp=predict(fit1,type="coefficients",s=lambda.sel1)
  coeff1=tmp@x
  names(coeff1)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  #coeff
  plotroc1(fit1,x,y,lambda.sel1,main=paste0(prefix," no penalty, glmnet,train"))
  plotroc1(fit1,x1,y1,lambda.sel1,main=paste0(prefix, " no penalty glmnet,test")) #0.733
  
  control=trainControl(method = "repeatedcv", number = 10, repeats = 10, 
                       returnResamp = "all", classProbs = TRUE, summaryFunction = twoClassSummary)
  set.seed(1)
  model <- train(x, y, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(.alpha =1, .lambda = seq(0.01, 0.1, by = 0.005)), trControl = control)
  
  print(plot(model, metric = "ROC"))
  allauc <- predict(model, newdata = x, type = "prob")
  tmp <- predict(model, newdata = x, type = "prob")
  plotroc3(tmp$Sensitive,y,main=paste0(prefix,"caret,glmnet,train"))
  tmp <- predict(model, newdata = x1, type = "prob")
  #colAUC(test, y1)
  plotroc3(tmp$Sensitive,y1,main=paste0(prefix,"caret,glmnet,test")) #0.732
  coeff2=predictors(model)
  test=varImp(model)
  
  if (useclinical==T)
  {
    trainingdata1=as.data.frame(cbind(y=y,x))
    #only include clinicals
    fm_clinical=as.formula(paste0("y~",paste0(colnames(trainingdata1)[2:idx_y],collapse="+")))
    #fm_clinical=as.formula(paste0("y~",paste(names(selected)[1:11],collapse="+")))
    glmfit_clinical=glm(fm_clinical,data=trainingdata1,family="binomial")
    plotroc(fit1=glmfit_clinical,trainingdata1,main=paste0(prefix, "glmfit,clinical,train"))
    
    testingdata1=as.data.frame(cbind(y=y1,x1))
    plotroc(fit1=glmfit_clinical,testingdata1,main=paste0(prefix,"glmfit,clinical,test"))
  }
  
  return(result=list(penaltycoef=coeff,nopenaltycoef=coeff1,roccoef=coeff2,allauc=allauc))
}

sel1=classify_plat3(data=data_mrna_commc,useclinical=T,num_clinical=4,prefix="mrna,useclinical,allgenes,")
sel2=classify_plat3(data=data_copynumber_commc,useclinical=T,num_clinical=4,prefix="copynumber,useclinical,allgenes,")
# sel3=classify_plat3(data=data_mrna_commc,useclinical=F,num_clinical=4,prefix="mrna,noclinical,allgenes,")
# sel4=classify_plat3(data=data_copynumber_commc,useclinical=T,num_clinical=4,prefix="copynumber,noclinical,allgenes,")

# only consider subset of genes
# top copynumber genes, with high expression
data=data_copynumber_commc
num_allclinical=12
traindata1=rbind(data$data1,data$data2)
avgmaggenes=colMeans(abs(traindata1[,(num_allclinical+1):ncol(traindata1)]))
avgmaggenes=avgmaggenes[order(avgmaggenes,decreasing=T)]
selgenes=names(avgmaggenes)[1:1000]
sel5=classify_plat3(data=data_copynumber_commc,useclinical=T,num_clinical=4,prefix="copynumber,useclinical,topexpression,",selgenes=selgenes)

data=data_mrna_commc
num_allclinical=12
traindata1=rbind(data$data1,data$data2)
avgmaggenes=colMeans(abs(traindata1[,(num_allclinical+1):ncol(traindata1)]))
avgmaggenes=avgmaggenes[order(avgmaggenes,decreasing=T)]
selgenes=names(avgmaggenes)[1:1000]
sel8=classify_plat3(data=data_mrna_commc,useclinical=T,num_clinical=4,prefix="mrna,useclinical,topexpression,",selgenes=selgenes)


#top copynumber genes have correlation with mrna
data=data_copynumber_commc
num_allclinical=12
traindata1=rbind(data$data1,data$data2)
data=data_mrna_commc
traindata2=rbind(data$data1,data$data2)
corrs=sapply((num_allclinical+1):ncol(traindata1),function(x){
  res=cor(traindata1[,x],traindata2[,x],use="complete")
})
names(corrs)=colnames(traindata1)[(num_allclinical+1):ncol(traindata1)]
corrs=corrs[order(abs(corrs),decreasing=T)]
selgenes=names(corrs)[1:1000]
sel6=classify_plat3(data=data_copynumber_commc,useclinical=T,num_clinical=4,prefix="copynumber,useclinical,topcorrelation,",selgenes=selgenes)
sel7=classify_plat3(data=data_mrna_commc,useclinical=T,num_clinical=4,prefix="mrna,useclinical,topcorrelation,",selgenes=selgenes)


classify_plat_rand=function(data=data_exon_rand,useclinical=T,num_clinical=4,prefix="copynumber,useclinical,",selgenes=NULL)
{
  num_clinical=4
  data=dataselgenes(data,selgenes)
  alldata=updateclinicalitem(data)
  #clinicals to remove before analysis
  remove_clinical=c("race","initial_pathologic_dx_year","vital_status","death_days_to","death_months_to","progression_free_survival","uselastcontact")
  alldata=alldata[,!colnames(alldata) %in% remove_clinical]
  
  #some genes have name with "-" character
  tmp=gsub("-","_",colnames(alldata))
  colnames(alldata)=tmp
  #to process colnames in methylation data which contain |
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
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(data$trainingdata)]=TRUE
  alldata=cbind(alldata,train=train)
  usecorrection=F
  if (usecorrection==T)
  {
    #mannualy correct those mismatches
    mannualclass=list(samples=c("TCGA-13-0803","TCGA-29-1774","TCGA-24-2029","TCGA-61-2113,TCGA-61-2009"),class=c("Sensitive","Resistant","Resistant","Resistant","Resistant"))
    for (i in 1:length(mannualclass[[1]]))
    {
      tmp=mannualclass[[1]][i]
      tmp1=which(rownames(alldata) %in% tmp)
      tmp2=mannualclass[[2]][i]
      alldata$platinumclass[tmp1]=tmp2
    }
    #remove samples having computed interval close to the boundary
    tmp=which(alldata$drug_interval_computed>5.9 & alldata$drug_interval_computed<6.1)
    if (length(tmp)>0)
    {
      alldata=alldata[-tmp,]
    }
  }
  #remove drug interval
  computed_interval=alldata[,"drug_interval_computed"]
  alldata=alldata[,-which(colnames(alldata)=="drug_interval_computed")]
  
  #interplate categorical clinical data
  interplate_flag=F
  if (interplate_flag==T)
  {
    alldata[,2]=interpl_factor(alldata[,2])
    alldata[,3]=interpl_factor(alldata[,3])
    alldata[,4]=interpl_factor(alldata[,4])
    
  }
  #remove grade 1(G1),grad4 (G4) in tumor_grade and Stage1 in clinical_stage
  alldata=alldata[alldata[,"tumor_grade"]!="G1",]
  alldata=alldata[alldata[,"tumor_grade"]!="G4",]
  alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
  alldata=alldata[alldata[,"clinical_stage"]!="Stage1",]
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  
  if (useclinical==T)
  {
    #form the design matrix for clinical data
    xfactors=model.matrix(as.formula(paste0("~",paste0(colnames(alldata)[1:num_clinical],collapse="+"))),data=alldata)[,-1]
    #some samples may be removed due to missing data in clinical
    tmp=sapply(rownames(xfactors),function(x){
      res=which(rownames(alldata)==x)
    })
    alldata1=alldata[tmp,]
    #combine clinical categorical data and other numeric data
    alldata2=cbind(xfactors,alldata1[,(num_clinical+1):ncol(alldata1)])
  }else
  {
    alldata2=alldata[,(num_clinical+1):ncol(alldata)]
  }
  colnames(alldata2)=gsub(" ","_",colnames(alldata2))
  
  idx_train=which(colnames(alldata2)=="train")
  idx_y=which(colnames(alldata2)=="platinumclass")
  #check clinical variables
  test=glm(as.formula(paste0("platinumclass~",paste0(colnames(alldata2)[1:(idx_y-1)],collapse="+"))),data=alldata2[1:255,],family="binomial")
  #summary(test)
  #   #multiple genes may have the save data (in copynumber),just keep one of them
  #   correlationMatrix <- cor(alldata2[,-c(idx_y,idx_train)])
  #   highlyCorrelated <- findCorrelation(correlationMatrix, names=T,cutoff=0.999)
  #   if (length(highlyCorrelated)>0)
  #   {
  #     alldata2=alldata2[,-which(colnames(alldata2) %in% highlyCorrelated)]
  #   }
  #   rm(correlationMatrix)
  #   idx_train=which(colnames(alldata2)=="train")
  #   mod1=as.formula(paste0("platinumclass~",paste0(colnames(alldata2)[1:ncol(xfactors)],collapse="+")))
  #   test=glm(mod1,family="binomial",data=alldata2[1:250,1:15])
  #   summary(test)
  #form training/testing data
  
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
  dim(x)
  trainingsamples=rownames(x)
  testingsamples=rownames(x1)
  
  #train
  if (useclinical==T)
  {
    penalty=rep(1,ncol(x))
    penalty[1:(idx_y-1)]=0
    trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=T,main=prefix)
    #     penalty=rep(1,ncol(x))
    #     penalty[(idx_y-1)]=0
    #     trainresult6=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=T) #only include x[,7]
  }else
  {
    trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=TRUE)
  }
  if (trainresult$lambda.1se==trainresult$lambseq[1])
  {
    print("use min lambda")
    lambda.sel=trainresult$lambda.min
  }else
  {
    print("use 1se lambda")
    lambda.sel=trainresult$lambda.1se
  }
  #lambda.sel=getlambda1sebyloo(x,y,numlambda=200,opt="binomial",opt1="penalty",idx_y)
  
  fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty,standardize=T,nlambda = 200)
  tmp=predict(fit,type="coefficients",s=lambda.sel)
  coeff=tmp@x
  names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  #coeff
  plotroc1(fit,x,y,lambda.sel,main=paste0(prefix," glmnet,train"))
  plotroc1(fit,x1,y1,lambda.sel,main=paste0(prefix, "glmnet,test")) #0.733
  
  trainresult1=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",penalty=rep(1,ncol(x)),standardize=T,main=prefix)
  if (trainresult1$lambda.1se==trainresult1$lambseq[1])
  {
    print("use min lambda")
    lambda.sel1=trainresult1$lambda.min
  }else
  {
    print("use 1se lambda")
    lambda.sel1=trainresult1$lambda.1se
  }
  
  fit1=glmnet(as.matrix(x),y,family="binomial",standardize=T,nlambda = 200)
  tmp=predict(fit1,type="coefficients",s=lambda.sel1)
  coeff1=tmp@x
  names(coeff1)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  #coeff
  plotroc1(fit1,x,y,lambda.sel1,main=paste0(prefix," no penalty, glmnet,train"))
  plotroc1(fit1,x1,y1,lambda.sel1,main=paste0(prefix, " no penalty glmnet,test")) #0.733
  allauc=computeauc(x,y,fit=fit1,lambdas=fit1$lambda)
  
  selected1=which(colnames(x1) %in% names(coeff1[2:length(coeff1)]))
  trainingdata1=as.data.frame(cbind(y=y,x[,selected1]))
  glmfit_train=glm(as.formula(y~.),data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit_train,trainingdata1,main=paste0(prefix,"glmfit_train,train"))
  testingdata1=as.data.frame(cbind(y=y1,x1[,selected1]))
  plotroc(fit1=glmfit_train,testingdata1,main=paste0(prefix,"glmfit_train,test"))
  
  
  glmfit_test=glm(as.formula(y~.),data=testingdata1,family="binomial")
  plotroc(fit1=glmfit_test,testingdata1,main=paste0(prefix,"glmfit,clinical,test,learnedcoeff"))
  
  control=trainControl(method = "repeatedcv", number = 10, repeats = 10, 
                       returnResamp = "all", classProbs = TRUE, summaryFunction = twoClassSummary)
  set.seed(1)
  model <- train(x, y, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(.alpha =1, .lambda = seq(0.01, 0.1, by = 0.005)), trControl = control)
  
  print(plot(model, metric = "ROC"))
  tmp <- predict(model, newdata = x, type = "prob")
  plotroc3(tmp$Sensitive,y,main=paste0(prefix,"caret,glmnet,train"))
  tmp <- predict(model, newdata = x1, type = "prob")
  #colAUC(test, y1)
  plotroc3(tmp$Sensitive,y1,main=paste0(prefix,"caret,glmnet,test")) #0.732
  coeff2=predictors(model)
  test=varImp(model)
  
  if (useclinical==T)
  {
    trainingdata1=as.data.frame(cbind(y=y,x))
    #only include clinicals
    fm_clinical=as.formula(paste0("y~",paste0(colnames(trainingdata1)[2:idx_y],collapse="+")))
    #fm_clinical=as.formula(paste0("y~",paste(names(selected)[1:11],collapse="+")))
    glmfit_clinical=glm(fm_clinical,data=trainingdata1,family="binomial")
    plotroc(fit1=glmfit_clinical,trainingdata1,main=paste0(prefix, "glmfit,clinical,train"))
    testingdata1=as.data.frame(cbind(y=y1,x1))
    plotroc(fit1=glmfit_clinical,testingdata1,main=paste0(prefix,"glmfit,clinical,test"))
    tmp=which(colnames(trainingdata1) %in% c("age","residual_disease_largest_noduleNo_Macroscopic_disease"))
    fm_clinical1=as.formula(paste0("y~",paste0(colnames(trainingdata1)[tmp],collapse="+")))
    glmfit_clinical1=glm(fm_clinical1,data=trainingdata1,family="binomial")
    plotroc(fit1=glmfit_clinical1,trainingdata1,main=paste0(prefix, "glmfit,clinical1,train"))
    plotroc(fit1=glmfit_clinical1,testingdata1,main=paste0(prefix,"glmfit,clinical1,test"))
  }
  return(result=list(penaltycoef=coeff,nopenaltycoef=coeff1,roccoef=coeff2,allauc=allauc))
}

sel1_rand=classify_plat_rand(data=data_mrna_commc_rand,useclinical=T,num_clinical=4,prefix="mrna,rand,useclinical,allgenes,")
sel2_rand=classify_plat_rand(data=data_copynumber_commc_rand,useclinical=T,num_clinical=4,prefix="copynumber,rand,useclinical,allgenes,")
#sel3_rand=classify_plat_rand(data=data_mrna_commc_rand,useclinical=F,num_clinical=4,prefix="mrna,rand,noclinical,allgenes,")
#sel4_rand=classify_plat_rand(data=data_copynumber_commc_rand,useclinical=T,num_clinical=4,prefix="copynumber,rand,noclinical,allgenes,")

# only consider subset of genes
# top copynumber genes, with high expression
data=data_copynumber_commc_rand
num_allclinical=12
traindata1=data$trainingdata
avgmaggenes=colMeans(abs(traindata1[,(num_allclinical+1):ncol(traindata1)]))
avgmaggenes=avgmaggenes[order(avgmaggenes,decreasing=T)]
selgenes=names(avgmaggenes)[1:1000]
sel5_rand=classify_plat_rand(data=data_copynumber_commc_rand,useclinical=T,num_clinical=4,prefix="copynumber,rand,useclinical,topexpression,",selgenes=selgenes)

#top copynumber genes have correlation with mrna
data=data_copynumber_commc_rand
num_allclinical=12
traindata1=data$trainingdata
data=data_mrna_commc_rand
traindata2=data$trainingdata
corrs=sapply((num_allclinical+1):ncol(traindata1),function(x){
  res=cor(traindata1[,x],traindata2[,x],use="complete")
})
names(corrs)=colnames(traindata1)[(num_allclinical+1):ncol(traindata1)]
corrs=corrs[order(abs(corrs),decreasing=T)]
selgenes=names(corrs)[1:1000]
sel6_rand=classify_plat_rand(data=data_copynumber_commc_rand,useclinical=T,num_clinical=4,prefix="copynumber,rand,useclinical,topcorrelation,",selgenes=selgenes)
sel7_rand=classify_plat_rand(data=data_mrna_commc_rand,useclinical=T,num_clinical=4,prefix="mrna,rand,useclinical,topcorrelation,",selgenes=selgenes)

combine2platforms=function(data1,data2)
{
  comsamples=rownames(data1)[rownames(data1) %in% rownames(data2)]
  tmp=sapply(1:length(comsamples),function(x){
    res=which(rownames(data1)==comsamples[x])
  })
  data1=data1[tmp,]
  tmp=sapply(1:length(comsamples),function(x){
    res=which(rownames(data2)==comsamples[x])
  })
  data2=data2[tmp,]
  data=cbind(data1,data2)
}
#use selected genes from all platforms
classify_plat4=function(wholedata=list(data1=data_exon,data2=data_4502a,data3=data_copynumber_fh,data4=data_copynumber_CGH1M,data5=data_methylation_fh),
                        selgenes=list(genes1=mrna_exon_selected,genes2=mrna_4502a_selected,genes3=copynumber_snp6_selected,genes4=copynumber_CGH1M_selected,genes5=methylation_selected),
                        prefs=c("mrnaexon","mrna4502a","copynumbersnp6","copynumberCGH1M","methylation"))
{
  #get all the data
  data=wholedata[[1]]
  traindatas=cbind(data$data1[,1:7],data$data1[,which(colnames(data$data1) %in% selgenes[[1]])])
  colnames(traindatas)[8:ncol(traindatas)]=paste0(prefs[1],"_",colnames(traindatas)[8:ncol(traindatas)])
  testdata=rbind(cbind(data$data2[,1:7],data$data2[,which(colnames(data$data2) %in% selgenes[[1]])]),cbind(data$data3[,1:7],data$data3[,which(colnames(data$data3) %in% selgenes[[1]])]))
  colnames(testdata)[8:ncol(testdata)]=paste0(prefs[1],"_",colnames(testdata)[8:ncol(testdata)])
  if (length(wholedata)>1)
  {
    for (i in 2:length(wholedata))
    {
      data=wholedata[[i]]
      traindata1=data$data1[,which(colnames(data$data1) %in% selgenes[[i]])]
      colnames(traindata1)=paste0(prefs[i],"_",colnames(traindata1))
      traindatas=combine2platforms(traindatas,traindata1)
      
      testdata1=rbind(data$data2[,which(colnames(data$data2) %in% selgenes[[i]])],data$data3[,which(colnames(data$data3) %in% selgenes[[i]])])
      colnames(testdata1)=paste0(prefs[i],"_",colnames(testdata1))
      testdata=combine2platforms(testdata,testdata1)
    }
  }
  #standardize using all the data
  alldata=rbind(traindatas,testdata)
  #some genes have name with "-" character
  tmp=gsub("-","_",colnames(alldata),fixed=T)
  colnames(alldata)=tmp
  tmp=gsub("|","_",colnames(alldata),fixed=T)
  colnames(alldata)=tmp
  #remove constant columns
  tmp1=sapply(8:ncol(alldata),function(x){
    res=FALSE
    if (var(alldata[,x])>0)
    {
      res=TRUE
    }
    return(res)
  })
  alldata=cbind(alldata[,1:7],alldata[8:ncol(alldata)][,tmp1])
  #scale those numeric data
  alldata[,c(1,8:ncol(alldata))]=scale(alldata[,c(1,8:ncol(alldata))])
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(traindatas)]=TRUE
  alldata=cbind(alldata,train=train)
  #remove drug interval
  alldata=alldata[,-6]
  #process clinical stage, merge A,B,C
  alldata[,3]=as.character(alldata[,3])
  tmp=which(alldata[,3] %in% c("Stage IA","Stage IB","Stage IC"))
  alldata[tmp,3]="Stage1"
  tmp=which(alldata[,3] %in% c("Stage IIA","Stage IIB","Stage IIC"))
  alldata[tmp,3]="Stage2"
  tmp=which(alldata[,3] %in% c("Stage IIIA","Stage IIIB","Stage IIIC"))
  alldata[tmp,3]="Stage3"
  tmp=which(alldata[,3] %in% c("Stage IV"))
  alldata[tmp,3]="Stage4"
  alldata[,3]=as.factor(alldata[,3])
  #process residual names
  alldata[,4]=as.character(alldata[,4])
  tmp=which(alldata[,4] %in% "No Macroscopic disease")
  alldata[tmp,4]="No_Macroscopic_disease"
  tmp=which(alldata[,4] %in% "1-10 mm")
  alldata[tmp,4]="1_10mm"
  tmp=which(alldata[,4] %in% "11-20 mm")
  alldata[tmp,4]="11_20mm"
  tmp=which(alldata[,4] %in% ">20 mm")
  alldata[tmp,4]="20mm_"
  alldata[,4]=as.factor(alldata[,4])
  #interplate categorical clinical data
  alldata[,2]=interpl_factor(alldata[,2])
  alldata[,3]=interpl_factor(alldata[,3])
  alldata[,4]=interpl_factor(alldata[,4])
  
  #form the design matrix for clinical data
  tmp=model.matrix(~.-1,data=alldata[,1:6])
  #remove variables with low variation
  tmp1=sapply(2:ncol(tmp),function(x){
    res=FALSE
    tmp2=min(table(tmp[,x]))
    if (tmp2<nrow(tmp)/100)
    {
      res=TRUE
    }
    return(res)
  })
  tmp=cbind(age=tmp[,1],tmp[,2:ncol(tmp)][,!tmp1])
  
  #lose of samples due to missing values
  tmp1=sapply(rownames(tmp),function(x){
    res=which(rownames(alldata)==x)
  })
  alldata=alldata[tmp1,]
  allx=cbind(tmp,alldata[,7:ncol(alldata)])
  #allx=model.matrix(~.-1,data=alldata)
  idx_train=which(colnames(allx)=="train")
  #allx[,idx_train]=as.logical(allx[,idx_train])
  idx_y=which(colnames(allx)=="platinumclassSensitive")
  #training data
  x=allx[which(allx[,idx_train]==1),]
  y=x[,idx_y]
  x=x[,-c(idx_y,idx_train)]
  trainingsamples=rownames(x)
  #testing data
  x1=allx[which(allx[,idx_train]==0),]
  y1=x1[,idx_y]
  x1=x1[,-c(idx_y,idx_train)]
  testingsamples=rownames(x1)
  
  #   #cmopute correlation between mrna and class, filter genes using correlation first
  #   corrs=sapply(1:ncol(x),function(i){
  #     res=cor(x[,i],y,use="complete")
  #   })
  #   #1 cor is 0, due to constant 
  #   corrs[is.na(corrs)]=0
  #   names(corrs)=colnames(x)
  #   corrs1=corrs[order(abs(corrs),decreasing=T)]
  #   #include clinical data
  #   selvars=unique(c(colnames(x)[1:idx_y],names(corrs1)[1:1000]))
  #   length(selvars)
  #   tmp=sapply(selvars,function(i){
  #     res=which(colnames(x)==i)
  #   })
  #   x=x[,tmp]
  x[,(idx_y+1):ncol(x)]=scale(x[,(idx_y+1):ncol(x)])
  #  x1=x1[,tmp]
  x1[,(idx_y+1):ncol(x1)]=scale(x1[,(idx_y+1):ncol(x1)])
  
  #plot(colMeans(x))
  #points(colMeans(x1),col="red")
  
  trainresult=selectlambda(x,y,numlambda=200,ntime=10,opt="binomial",opt1="penalty",idx_y)
  lambda.sel=trainresult$lambda.sel
  #lambda.sel=getlambda1sebyloo(x,y,numlambda=200,opt="binomial",opt1="penalty",idx_y)
  
  penalty=rep(1,ncol(x))
  penalty[1:idx_y]=0
  
  fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty)
  plotroc1(fit,x,y,lambda.sel,main="glmnet,train")
  
  #use glm model
  selected <- which(as.matrix(coef(fit,s=lambda.sel))[,1]!=0)
  #remove intercept
  if (names(selected)[1]=="(Intercept)")
  {
    selected=selected[-1]-1
  }else
  {
    selected=selected-1
  }
  #mrna_exon_selected=names(selected)[(idx_y+1):length(selected)]
  #mrna_4502a_selected=names(selected)[(idx_y+1):length(selected)]
  #copynumber_snp6_selected=names(selected)[(idx_y+1):length(selected)]
  #copynumber_CGH1M_selected=names(selected)[(idx_y+1):length(selected)]
  #methylation_selected=names(selected)[(idx_y+1):length(selected)]
  
  
  #build a logistic regression model
  trainingdata1=as.data.frame(cbind(y=y,x))
  fm=as.formula(paste0("y~",paste(names(selected),collapse="+")))
  glmfit=glm(fm,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit,trainingdata1,main="glmfit, train")
  #only include clinicals
  fm_clinical=as.formula(paste0("y~",paste(names(selected)[1:12],collapse="+")))
  glmfit_clinical=glm(fm_clinical,data=trainingdata1,family="binomial")
  plotroc(fit1=glmfit_clinical,trainingdata1,main="glmfit,clinical,train")
  
  #   #use variables with small p values
  #   tmp=summary(glmfit)
  #   idx=which(tmp$coefficients[,4]<0.05)
  #   if (names(idx)[1] == "(Intercept)")
  #     idx=idx[-1]
  #   fm1=as.formula(paste0("y~",paste(names(idx),collapse="+")))
  #   glmfit1=glm(fm1,data=trainingdata1,family="binomial")
  #   plotroc(fit1=glmfit1,trainingdata1)
  #   #use lda
  #   ldafit=lda(fm,data=trainingdata1)
  #   plotroc(fit1=ldafit,trainingdata1,opt="lda")
  #   ldafit1=lda(fm1,data=trainingdata1)
  #   plotroc(fit1=ldafit1,trainingdata1,opt="lda")
  #   #use qda
  #   #   qdafit=qda(fm,data=trainingdata1)
  #   #   plotroc(fit1=qdafit,trainingdata1,opt="qda")
  #   #   qdafit1=qda(fm1,data=trainingdata1)
  #   #   plotroc(fit1=qdafit1,trainingdata1,opt="lda")
  #   
  #   #for testing data
  #   #   for (i in 1:ncol(x1))
  #   #   {
  #   #     x1[,i] <- ifelse(is.na(x1[,i]),mean(x1[!is.na(x[,i]),i]),x1[,i])
  #   #   }
  
  plotroc1(fit,x1,y1,lambda.sel,main="glmnet,test")
  testingdata1=as.data.frame(cbind(y=y1,x1))
  plotroc(fit1=glmfit,testingdata1,main="glmfit,test")
  #variables selected by training, parameters were refined by testing data
  glmfit_test=glm(fm,data=testingdata1,family="binomial")
  plotroc(fit1=glmfit_test,testingdata1,main="glmfit,clinical,test")
  #   #check if the 
  #   plotroc(fit1=glmfit_clinical,testingdata1)
  #   plotroc(fit1=glmfit1,testingdata1)
  #   plotroc(fit1=ldafit,testingdata1,opt="lda")
  #   plotroc(fit1=ldafit1,testingdata1,opt="lda")
  #   #   plotroc(fit1=qdafit,testingdata1,opt="qda")
  #   #   plotroc(fit1=qdafit1,testingdata1,opt="qda")
}

classify_plat3(wholedata=list(data1=data_exon,data3=data_copynumber_fh,data5=data_methylation_fh),
               selgenes=list(genes1=mrna_exon_selected,genes3=copynumber_snp6_selected,genes5=methylation_selected),
               prefs=c("mrnaexon","copynumbersnp6","methylation"))

classify_plat3(wholedata=list(data1=data_exon,data3=data_copynumber_fh),
               selgenes=list(genes1=mrna_exon_selected,genes3=copynumber_snp6_selected),
               prefs=c("mrnaexon","copynumbersnp6"))

computeauc=function(x,y,fit=fit,lambdas)
{
  res=data.frame(lambda=lambdas,selnum=rep(0,length(lambdas)),auc1=rep(NA,length(lambdas)),auc2=rep(NA,length(lambdas)))
  for (i in 1:length(lambdas))
  {
    lambda=lambdas[i]
    selected <- which(as.matrix(coef(fit,s=lambda))[,1]!=0)
    #remove intercept
    if (names(selected)[1]=="(Intercept)")
    {
      selected=selected[-1]-1
    }else
    {
      selected=selected-1
    }
    if (length(selected)>0)
    {
      res[i,2]=length(selected)
      #build a logistic regression model
      trainingdata1=as.data.frame(cbind(y=y,x))
      fm=as.formula(paste0("y~",paste(names(selected),collapse="+")))
      glmfit=glm(fm,data=trainingdata1,family="binomial")
      res[i,3]=plotroc(fit1=glmfit,trainingdata1,plotflag=0)
      testingdata1=as.data.frame(cbind(y=y1,x1))
      res[i,4]=plotroc(fit1=glmfit,testingdata1,plotflag=0)
    }
    
  }
  return(res)
}