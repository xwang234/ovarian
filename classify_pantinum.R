#!/usr/bin/env Rscript

#load the training/testing data
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/mRNA_platinum_classification.RData")

library(glmnet)
library(MASS)
library(ROCR)
library(caret)
library(pROC)
library(caTools)
#select lambda using crossvalidation
traindata=function(x,y,numlambda=200,ntime=100,opt="binomial",penalty=NULL,standardize=TRUE,main="",nfolds=10)
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
      #       cvfit1=cv.glmnet(as.matrix(allx[,-c(idx_y,idx_train)]),allx[,idx_y],nlambda=numlambda,nfolds=10,penalty.factor=penalty)
      #       cvfit1=cv.glmnet(as.matrix(allx[,-c(idx_y,idx_train)]),allx[,13],family="binomial",nlambda=numlambda,nfolds=10,penalty.factor=penalty)
      #       plot(cvfit1)
      cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=nfolds,penalty.factor=penalty,standardize=standardize)
      #cvfit <- cv.glmnet(as.matrix(x),y,family="binomial",nlambda=numlambda,nfolds=length(y),penalty.factor=penalty)
    }
    if (opt=="numeric")
    {
      set.seed(i)
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

#plot roc functions, fit1:logistic regression model, linear discriminat analysis model
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

#plot roc functions, fit: glmnet model
plotroc1=function(fit,x,y,lambda.sel=0.05,plotflag=1,main="")
{
  pfit = predict(fit,as.matrix(x),s=lambda.sel,type="response")
  yy=y[!is.na(pfit)]
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

#plot roc functions,fit1: principle component regression model, other model
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

#plot roc functions, no model, use predicted values instead
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

#compute AUC using all the lambdas, used for glmnet model
computeauc=function(x=x,y=y,x1=x1,y1=y1,fit=fit,lambdas)
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


#interplate missing factor values,not used
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


#update cilinical items of "clinical stage" and "residual disease"
updateclinicalitem=function(data=data_mrna)
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

classify_plat3=function(data=data_mrna,useclinical=T,num_clinical=4,prefix="mrna,useclinical,",selgenes=NULL)
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
    
    #some genes have name with "-" character, change to "_"
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
    
    usecorrection=F
    if (usecorrection==T)
    {
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
    #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
    alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
    alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
    alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
    alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
    
    if (useclinical==T)
    {
      #form the design matrix for clinical data
      xfactors=model.matrix(as.formula(paste0("~",paste0(colnames(alldata)[1:num_clinical],collapse="+"))),data=alldata)[,-1]
      #some samples may be removed due to missing data in clinical, more than 30 removed for missing in residual
      tmp=sapply(rownames(xfactors),function(x){
        res=which(rownames(alldata)==x)
      })
      alldata1=alldata[tmp,]
      #combine clinical categorical data and other numeric data
      alldata2=cbind(xfactors,alldata1[,(num_clinical+1):ncol(alldata1)])
    }else
    {
      #without clinicals
      alldata2=alldata[,(num_clinical+1):ncol(alldata)]
    }
    colnames(alldata2)=gsub(" ","_",colnames(alldata2))
    
    #the column indices of train indicator and outcome
    idx_train=which(colnames(alldata2)=="train")
    idx_y=which(colnames(alldata2)=="platinumclass")
    #check clinical variables
    test=glm(as.formula(paste0("platinumclass~",paste0(colnames(alldata2)[1:(idx_y-1)],collapse="+"))),data=alldata2[alldata2[,idx_train]==T,],family="binomial")
    #summary(test)
#     #multiple genes may have the same data (in copynumber),just keep one of them
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
  
  #train, keep all the clinicals setting penalty=0
  if (useclinical==T)
  {
    penalty=rep(1,ncol(x))
    penalty[1:(idx_y-1)]=0
    trainresult=traindata(x,y,numlambda=200,ntime=10,opt="binomial",penalty=penalty,standardize=T,main=prefix,nfolds=5)
    if (trainresult$lambda.1se==trainresult$lambseq[1])
    {
      print("use min lambda")
      lambda.sel=trainresult$lambda.min
    }else
    {
      print("use 1se lambda")
      lambda.sel=trainresult$lambda.1se
    }
    fit=glmnet(as.matrix(x),y,family="binomial",penalty.factor=penalty,standardize=T,nlambda = 200)
    tmp=predict(fit,type="coefficients",s=lambda.sel)
    coeff=tmp@x
    names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
    #coeff
    plotroc1(fit,x,y,lambda.sel,main=paste0(prefix," glmnet,train"))
    plotroc1(fit,x1,y1,lambda.sel,main=paste0(prefix, "glmnet,test"))
    #try without forcing keep clinicals
    trainresult1=traindata(x,y,numlambda=200,ntime=10,opt="binomial",penalty=rep(1,ncol(x)),standardize=T,main=prefix)
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
    plotroc1(fit1,x1,y1,lambda.sel1,main=paste0(prefix, " no penalty glmnet,test"))
    allauc=computeauc(x,y,x1,y1,fit=fit1,lambdas=fit1$lambda)
    #the model only considering clinicals
    trainingdata1=as.data.frame(cbind(y=y,x))
    #only include clinicals
    fm_clinical=as.formula(paste0("y~",paste0(colnames(trainingdata1)[2:idx_y],collapse="+")))
    #fm_clinical=as.formula(paste0("y~",paste(names(selected)[1:11],collapse="+")))
    glmfit_clinical=glm(fm_clinical,data=trainingdata1,family="binomial")
    plotroc(fit1=glmfit_clinical,trainingdata1,main=paste0(prefix, "glmfit,clinical,train"))
    testingdata1=as.data.frame(cbind(y=y1,x1))
    plotroc(fit1=glmfit_clinical,testingdata1,main=paste0(prefix,"glmfit,clinical,test"))
    
  }else
  {
    coeff=NULL
    trainresult1=traindata(x,y,numlambda=200,ntime=10,opt="binomial",penalty=rep(1,ncol(x)),standardize=T,main=prefix)
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
    allauc=computeauc(x,y,x1,y1,fit=fit1,lambdas=fit1$lambda)
  }
  
#use caret
  control=trainControl(method = "repeatedcv", number = 10, repeats = 10, 
                       returnResamp = "all", classProbs = TRUE, summaryFunction = twoClassSummary)
  set.seed(1)
  model <- train(x, y, method = "glmnet", metric = "ROC", tuneGrid = expand.grid(.alpha =1, .lambda = seq(0.01, 0.1, by = 0.005)), trControl = control)
  
  print(plot(model, metric = "ROC"))
  allauc1 <- predict(model, newdata = x, type = "prob")
  tmp <- predict(model, newdata = x, type = "prob")
  plotroc3(tmp$Sensitive,y,main=paste0(prefix,"caret,glmnet,train"))
  tmp <- predict(model, newdata = x1, type = "prob")
  #colAUC(test, y1)
  plotroc3(tmp$Sensitive,y1,main=paste0(prefix,"caret,glmnet,test")) #0.732
  coeff2=predictors(model)
  test=varImp(model)
  
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

