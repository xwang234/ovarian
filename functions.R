#!/usr/bin/env Rscript
library(glmnet)
library(gap) #Manhattan plot
library(GenomicRanges) #find overlaps
library(ROCR)
library(pROC)
library(gdata) #read.xls
#alias genes used to map genes with alias symbol names
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genesymbolposition.RData")

#functions used to draw ROC------------------------------------------------------------------------------------------------------

#creat data considering only from a set of genes, considering alias gene names
#if selgenes came outside of the data, use searchliasgene=T, but time consumming.
dataselgenes=function(data,selgenes,num_allclinical=13,searchaliasgene=F)
{
  #if not specify selgenes, use all the data
  if (is.null(selgenes))
  {
    data1=data
  }else
  {
      if (searchaliasgene==T)
      {
        map_datagenes=mapaliasgenes(colnames(data$data1))
        map_datagenes[1:num_allclinical]=colnames(data$data1)[1:num_allclinical]
        map_selgenes=mapaliasgenes(selgenes)
        idxgenes=which(map_datagenes %in% map_selgenes)
      }else
      {
        idxgenes=which(colnames(data$data1) %in% selgenes)
      }
      data1=list(data1=cbind(data$data1[,1:num_allclinical],data$data1[,idxgenes]),
      data2=cbind(data$data2[,1:num_allclinical],data$data2[,idxgenes]))
  }
  return(data1)
}

#update cilinical items of "clinical stage" and "residual disease", combine training data testing data into a dataframe
formwholedataframe=function(data=data_mrna)
{
  if (class(data)=="list")
  {
    alldata=NULL
    for (i in 1:2)
    {
      alldata=rbind(alldata,data[[i]])
    }
  }else
  {
    alldata=data
  }
  
  #process clinical stage, merge A,B,C
  if ("clinical_stage" %in% colnames(alldata))
  {
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
  }
  if ("residual_disease_largest_nodule" %in% colnames(alldata))
  {
    alldata[,"residual_disease_largest_nodule"]=as.character(alldata[,"residual_disease_largest_nodule"])
    tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "No Macroscopic disease")
    alldata[tmp,"residual_disease_largest_nodule"]="No_Macroscopic_disease"
    tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "1-10 mm")
    alldata[tmp,"residual_disease_largest_nodule"]="1_10mm"
    tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "11-20 mm")
    alldata[tmp,"residual_disease_largest_nodule"]="11_20mm"
    tmp=which(alldata[,"residual_disease_largest_nodule"] %in% ">20 mm")
    alldata[tmp,"residual_disease_largest_nodule"]="20mm_"
#     idxna=is.na(alldata$residual_disease_largest_nodule)
#     alldata$residual_disease_largest_nodule[idxna]="missing"
    alldata[,"residual_disease_largest_nodule"]=as.factor(alldata[,"residual_disease_largest_nodule"])
  }
  if ("tumor_grade" %in% colnames(alldata))
  {
    #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
    alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
    alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  }
  
  #clinicals to remove before analysis, they are not useful for the model. only age and residual_disease are kept
  remove_clinicals=c("tumor_grade","clinical_stage","drug_interval_computed","race","initial_pathologic_dx_year","vital_status","death_months_to","treatment_outcome_first_course","progression_free_survival","uselastcontact")
  alldata=alldata[,!colnames(alldata) %in% remove_clinicals]
  return(alldata)
}

#select glmnet lambda
library(glmnet)
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
    yy=as.numeric(as.factor(y[!is.na(pfit)]))
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
  text(0.4,0.8,paste0("Training AUC=",res_train$auc),cex=1)
  text(0.4,0.25,paste0("Validation AUC=",res_test$auc),cex=1)
}

#calculate AUC
calauc=function(pfit_train,y)
{
  rocdata=function(pfit,y)
  {
    yy=as.numeric(as.factor(y[!is.na(pfit)]))
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
    auc <- format(as.numeric(auc.tmp@y.values),digits=5)
    result=list(pred=pred,auc=auc)
  }
  res_train=rocdata(pfit_train,y)
  res=as.numeric(res_train$auc)
}


#draw ROC using glm model based on selected variables
drawroc_onselectedgenes=function(data=data_copynumber_hrgenes,selgenes=selgenes)
{
  selgenes=intersect(selgenes,colnames(data$data1))
  idx=which(colnames(data$data1) %in% selgenes)
  traindata=data$data1[,c(1:13,idx)]
  testdata=data$data2[,c(1:13,idx)]
  traindata[,"platinumclass"]=as.factor(traindata[,"platinumclass"])
  testdata[,"platinumclass"]=as.factor(testdata[,"platinumclass"])
  colnames(traindata)=colnames(testdata)=gsub("|","___",colnames(traindata),fixed=T)
  selgenes=gsub("|","___",selgenes,fixed=T)
  fm=glm(as.formula(paste0("platinumclass~",paste0(selgenes,collapse="+"))),data=traindata,family="binomial")
  #remove collinear covariates:
  if (length(fm$coefficients)>fm$rank)
  {
    warning("some collinear covariates were removed")
    tmp=summary(fm)
    selgenes1=selgenes[selgenes %in% rownames(tmp$coefficients)]
    fm=glm(as.formula(paste0("platinumclass~",paste0(selgenes1,collapse="+"))),data=traindata,family="binomial")
  }
  print(summary(fm))
  pfit_train=predict(fm,type="response")
  pfit_test=predict(fm,newdata=testdata,type="response")
  plotroc2(pfit_train,traindata$platinumclass,pfit_test,testdata$platinumclass,main="ROC, without clinical")
#   traindata=formwholedataframe(traindata)
#   testdata=formwholedataframe(testdata)
#   fm=glm(as.formula(paste0("platinumclass~",paste0(selgenes,collapse="+"),"+residual_disease_largest_nodule")),data=traindata,family="binomial")
#   if (length(fm$coefficients)>fm$rank)
#   {
#     warning("some collinear covariates were removed")
#     tmp=summary(fm)
#     selgenes1=selgenes[selgenes %in% rownames(tmp$coefficients)]
#     fm=glm(as.formula(paste0("platinumclass~",paste0(selgenes1,collapse="+"),"+residual_disease_largest_nodule")),data=traindata,family="binomial")
#   }
#   pfit_train=predict(fm,type="response")
#   pfit_test=predict(fm,newdata=testdata,type="response")
#   plotroc2(pfit_train,traindata$platinumclass,pfit_test,testdata$platinumclass,main="ROC, with clinical")
  colnames(traindata)=colnames(testdata)=gsub("___","|",colnames(traindata),fixed=T)
  return(result=list(traindata=traindata,testdata=testdata))
}

#main function to draw ROC curve:
classify_platinum=function(data=data_copynumber,platform="copynumber",selgenes=NULL,includeclinical=T)
{
  #if genes are specified (as not NULL), use data only from selgenes.
  data=dataselgenes(data,selgenes)
  alldata=rbind(data$data1,data$data2)
  #add indicator of train/test as the last column
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(data$data1)]=TRUE
  alldata=cbind(alldata,train=train)
  
  #update stage and residual, combine training and testing into one dataframe
  alldata=formwholedataframe(alldata)
  
#   #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
#   alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
#   alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  
#   #clinicals to remove before analysis, they are not useful for the model. only age and residual_disease are kept
#   remove_clinicals=c("tumor_grade","clinical_stage","drug_interval_computed","race","initial_pathologic_dx_year","vital_status","death_months_to","treatment_outcome_first_course","progression_free_survival","uselastcontact")
#   alldata=alldata[,!colnames(alldata) %in% remove_clinicals]
  
  #number of clinical variables to keep (age,residual)
  num_clinical=2
  
  #some genes have name with "-" character,replace it to "_". Or it would cause problem in modeling
  colnames(alldata)=gsub("-","_",colnames(alldata))
  
  #with clinicals or not
  if (includeclinical==T)
  {
    #form the design matrix for clinical data
    xfactors=model.matrix(as.formula(paste0("~",paste0(colnames(alldata)[1:num_clinical],collapse="+"))),data=alldata)[,-1]
    #some samples may be removed due to missing data in clinical, more than 30 removed for missing in residual
    idx=sapply(rownames(xfactors),function(x){
      res=which(rownames(alldata)==x)
    })
    alldata1=alldata[idx,]
    #combine clinical categorical data and other numeric data
    alldata2=cbind(xfactors,alldata1[,(num_clinical+1):ncol(alldata1)])
  }else #without clinical data
  {
    alldata2=alldata[,(num_clinical+1):ncol(alldata)]
  }
  
  #the column indices of train indicator and outcome
  idx_train=which(colnames(alldata2)=="train")
  #the column of output
  idx_y=which(colnames(alldata2)=="platinumclass")
  
  #form training data
  ally=alldata2[,idx_y]
  x=alldata2[which(alldata2[,idx_train]==1),]
  y=ally[which(alldata2[,idx_train]==1)]
  x=x[,-c(idx_y,idx_train)]
  
  #form testing data
  x1=alldata2[which(alldata2[,idx_train]==0),]
  y1=ally[which(alldata2[,idx_train]==0)]
  x1=x1[,-c(idx_y,idx_train)]
  
  #remove constant genes in training set
  idx=sapply(idx_y:ncol(x),function(i){
    res=T
    if (var(x[,i],na.rm=T)==0)
    {
      res=F
    }
    return(res)
  })
  idx=c(rep(T,(idx_y-1)),idx)
  x=x[,idx]
  x1=x1[,idx]
  
  penalty=rep(1,ncol(x))
  if ("residual_disease_largest_noduleNo_Macroscopic_disease" %in% colnames(x))
  {
    penalty[which(colnames(x)=="residual_disease_largest_noduleNo_Macroscopic_disease")]=0
  }
  
 #select lambda
  trainresult=selectlambda(x,y,numlambda=100,ntime=20,opt="binomial",penalty=penalty,main=platform,nfolds=10)
  #lambda.sel=trainresult$lambda.1se
  lambda.sel=trainresult$lambda.min
  
  fit=glmnet(as.matrix(x),y,family="binomial",nlambda = 100,penalty.factor=penalty)
  tmp=predict(fit,type="coefficients",s=lambda.sel)
  coeff=tmp@x
  names(coeff)=c("intercept",colnames(x)[tmp@i[tmp@i>0]])
  selectedgenes=names(coeff)[! names(coeff) %in% c("intercept","age","residual_disease_largest_noduleNo_Macroscopic_disease")]
  
  pfit_train = predict(fit,as.matrix(x),s=lambda.sel,type="response")
  pfit_test = predict(fit,as.matrix(x1),s=lambda.sel,type="response")
  
  plotroc2(pfit_train,y,pfit_test,y1,main=paste0("ROC, ",length(selectedgenes), " genes selected"))
  return(list(coeff=coeff,x=x,y=y,x1=x1,y1=y1,pfit_train=pfit_train,pfit_test=pfit_test))
}

#function to compute correlation between mRNA and CNA
compute_corr_mrna_cna=function(data1=data_copynumber,data2=data_mrna)
{
  alldata1=formwholedataframe(data1)
  #add indicator of train/test as the last column
  train=rep(FALSE,nrow(alldata1))
  train[1:nrow(data1$data1)]=TRUE
  alldata1=cbind(alldata1,train=train)
#   #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
#   alldata1=alldata1[is.na(alldata1[,"tumor_grade"]) | (! is.na(alldata1[,"tumor_grade"]) & alldata1[,"tumor_grade"]!="G1"),]
#   alldata1=alldata1[is.na(alldata1[,"clinical_stage"]) | (! is.na(alldata1[,"clinical_stage"]) & alldata1[,"clinical_stage"]!="Stage1"),]
  data1=alldata1[alldata1[,"train"]==T,4:ncol(alldata1)] #only use train data
    
  alldata2=formwholedataframe(data2)
  #add indicator of train/test as the last column
  train=rep(FALSE,nrow(alldata2))
  train[1:nrow(data2$data1)]=TRUE
  alldata2=cbind(alldata2,train=train)
#   #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
#   alldata2=alldata2[is.na(alldata2[,"tumor_grade"]) | (! is.na(alldata2[,"tumor_grade"]) & alldata2[,"tumor_grade"]!="G1"),]
#   alldata2=alldata2[is.na(alldata2[,"clinical_stage"]) | (! is.na(alldata2[,"clinical_stage"]) & alldata2[,"clinical_stage"]!="Stage1"),]
  data2=alldata2[alldata2[,"train"]==T,4:ncol(alldata2)]
  comsamples=intersect(rownames(data1),rownames(data2))
  comgenes=intersect(colnames(data1),colnames(data2))
  idxrows=sapply(1:length(comsamples),function(i){
    res=which(rownames(data1)==comsamples[i])
  })
  data1=data1[idxrows,]
  idxrows=sapply(1:length(comsamples),function(i){
    res=which(rownames(data2)==comsamples[i])
  })
  data2=data2[idxrows,]
  idxcols=sapply(1:length(comgenes),function(i){
    res=which(colnames(data1)==comgenes[i])
  })
  data1=data1[,idxcols]
  idxcols=sapply(1:length(comgenes),function(i){
    res=which(colnames(data2)==comgenes[i])
  })
  data2=data2[,idxcols]
  res=sapply(1:ncol(data1),function(i){
    res1=cor(data1[,i],data2[,i],use="complete")
  })
  names(res)=colnames(data1)
  return(res)
}


#functions to compute p-values-------------------------------------------------------------------------------
#compute p-value of a gene z --> outcome y
calpvalues=function(z,y)
{
  fit=glm(y~z,family="binomial")
  summaryres=coef(summary(fit))
  #if z is constant, no pvalue would be found, and just has the result for intercept
  if (nrow(summaryres)>1)
  {
    res=coef(summary(fit))[2,4]
  }else
  {
    #no p-value found
    res=NA
  }
  return(res)
}
#parallel computing, data is the X matrix
mpi_calpvalues=function(data,y,filename=NULL,njob=100)
{
  
  res1=rep(NA,ncol(data))
  nrun <- ceiling(ncol(data)/1000)
  print(paste0("total runs: ",nrun))
  for (j in 1:nrun){
    if (j %% 10==0) cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):ncol(data)
    z=data[,cseq]
    res=mpi.parCapply(X=z,FUN=calpvalues,y=y,job.num=njob)
    res1[cseq]=res
  }
  #res1=res1[1:ncol(data)]
  res2=data.frame(pvalues=res1)
  rownames(res2)=colnames(data)
  if (!is.null(filename))
  {
    write.table(res2,file=filename,quote=F,sep="\t",row.names=T,col.names=T)
  }
  return(res2)
}

#for only platinumclass added data (1st column is platinumclass)
compute_pvalues1=function(data)
{
  pvalues=sapply(2:ncol(data),function(i){
    fm=glm(as.factor(data[,"platinumclass"])~as.numeric(as.character(data[,i])),family = "binomial")
    res=NA
    if (nrow(coef(summary(fm)))>1)
    {
      res=coef(summary(fm))[2,4]
    }
    names(res)=colnames(data)[i-1]
    return(res)
  })
  return(pvalues)
}

#main function
#return the pvalue as a numeric array, the pvalues generated in permutation are stored in a dataframe and is saved as a text file.
#the mpi verion is compute_pvalues_genes_platinumclass.R
compute_pvalues=function(data=data_copynumber,platform="copynumber",runpermutation=F,usetrain=F)
{
  #update stage and residual
  alldata=formwholedataframe(data)
  #indicator of if a sample is for training
  train=rep(T,nrow(alldata))
  train[(nrow(data$data1)+1):nrow(alldata)]=F
  alldata=cbind(alldata,train=train)
#   #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
#   alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
#   alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  #number of training, used for the case if we only consider training samples
  numtrain=sum(alldata[,"train"]==T)
  #if only works on train data
  if (usetrain==T)
  {
    alldata=alldata[1:numtrain,]
    platform=paste0(platform,"_train")
  }
  #remove train column
  alldata=alldata[,colnames(alldata)!="train"]
  idx_y=which(colnames(alldata)=="platinumclass")
  #outcome
  y=alldata[,idx_y]
  #gene matrix
  data=alldata[,(idx_y+1):ncol(alldata)]
  colnames(data)=gsub("-","_",colnames(data))
  #loop to comupte p-values
  pvalues=sapply(1:ncol(data),function(i){
    res=calpvalues(data[,i],y)
  })
  names(pvalues)=colnames(data)
  
  #to run permutation,time consumming, may take 2.5 hours
  if (runpermutation==T)
  {
    print("start counting pvalues of permutations")
    permutationp=NULL
    #number of permutation
    numit=100
    for (i in 1:numit)
    {
      cat(i,"..")
      set.seed(i+1001)
      #shuffle the output
      y1=y[sample(length(y))]
      names(y1)=names(y)
      #pvalues of pumutation data
      ppvalues=sapply(1:ncol(data),function(i){
        res=calpvalues(data[,i],y1)
      })
      #order the p values
      ppvalues=ppvalues[order(ppvalues)]
      #it is a dataframe, each column stores a result from a permuation
      permutationp=cbind(permutationp,ppvalues)
    }
    filename2=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,"_permutation100.txt")
    write.table(permutationp,file=filename2,sep="\t",quote=F,row.names=F,col.names=F)
  }
  return(pvalues)
}

#function to compute fdr--------------------------------------------------------------------------------------
#the mpi version is compute_fdr_permutation.R
# sapply_compute_qvalue_permutation=function(pvalues,pvalues_permutation,numit=100,platform)
# {
#   res=sapply(pvalues,function(pvalue){
#     res=sum(pvalue>=pvalues_permutation,na.rm=T)/numit/sum(pvalue>=pvalues,na.rm=T)
#   })
#   res2=data.frame(pvalues=pvalues,qvalues=res)
#   rownames(res2)=names(pvalues)
#   filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",platform,"_permutation100.txt")
#   write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")
#   names(res)=names(pvalues)
#   return(res)
# }

#monotone q values
compute_qvalue_permutation=function(pvalues,pvalues_permutation,numit=1000,platform)
{
  #sort pvalues decreasingly
  genenames=rownames(pvalues)
  idx=order(pvalues[,1],decreasing=T)
  pvalues=pvalues[idx,1]
  names(pvalues)=genenames[idx]
  idxnoNA=which(!is.na(pvalues))
  qvalues=rep(NA,length(pvalues))
  for (i in idxnoNA)
  {
    qvalue=sum(pvalues_permutation<=pvalues[i],na.rm=T)/numit/sum(pvalues<=pvalues[i],na.rm=T)
    if (i==idxnoNA[1])
    {
      #the first q value
      qvalues[i]=qvalue
    }else
    {
      #make q values monotone
      qvalues[i]=min(qvalues[i-1],qvalue,na.rm=T)
    }
  }
  res2=data.frame(pvalues=pvalues,qvalues=qvalues)
  rownames(res2)=names(pvalues)
  filename=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",platform,"_permutation",numit,".txt")
  write.table(res2,file=filename,row.names=T,col.names=T,quote=F,sep="\t")

  return(res2)
}
#functions to draw Manhattan plot-----------------------------------------------------------------------------
#sort table based on chromosome
sortgenetable=function(genetable)
{
  genetable$chr=as.character(genetable$chr)
  genetable$chr=gsub("chr","",genetable$chr)
  genetable$chr=gsub("23","X",genetable$chr)
  genetable$chr=gsub("24","Y",genetable$chr)
  if (class(genetable$start[1])=="factor")
  {
    genetable$start=as.numeric(as.character(genetable$start))
  }
  if (class(genetable$start[1])=="character")
  {
    genetable$start=as.numeric(genetable$start)
  }
 
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genetable)))
  for (chr in chrs)
  {
    tmptable=genetable[which(genetable$chr==chr),]
    tmptable=tmptable[order(tmptable$start),]
    res=rbind(res,tmptable)
  }
  return(res)
}

#the function mhtplot has a problem in placing chomosomes on x axis, the following is the corrected one:
mhtplot1=function (data, control = mht.control(), hcontrol = hmht.control(),  ...) 
{
  for (p in c("grid")) {
    if (length(grep(paste("^package:", p, "$", sep = ""), 
                    search())) == 0) {
      if (!require(p, quietly = TRUE, character.only = TRUE)) 
        warning(paste("mhtplot needs package `", p, "' to be fully functional; please install", 
                      sep = ""))
    }
  }
  
  data2 <- data[!apply(is.na(data), 1, any), ]
  n2 <- dim(data2[1])
  chr <- data2[, 1]
  pos <- newpos <- data2[, 2]
  p <- data2[, 3]
  tablechr <- table(chr)
  #this part has problem, the orders is chr1,chr10,chr2,...; should change to chr1,chr2,...,chr10,... 
  #allchr <- as.vector(tablechr)
  uniqchr=unique(chr)
  allchr=NULL
  for (i in 1:length(uniqchr))
  {
    allchr=c(allchr,sum(chr==uniqchr[i]))
  }
  n.chr <- length(allchr)
  type <- control$type
  usepos <- control$usepos
  logscale <- control$logscale
  base <- control$base
  cutoffs <- control$cutoffs
  colors <- control$colors
  labels <- control$labels
  srt <- control$srt
  gap <- control$gap
  pcex <- control$cex
  yline <- control$yline
  xline <- control$xline
  colorlist <- colors()
  if (is.null(colors)) 
    colors <- sample(colorlist, n.chr)
  tablechr <- unique(chr)
  if (is.null(labels)) 
    labels <- tablechr
  if (is.null(gap)) 
    gap <- 0
  if (!is.null(hcontrol$data)) {
    hdata <- hcontrol$data
    hdata2 <- hdata[!apply(is.na(hdata), 1, any), ]
    hchr <- hdata2[, 1]
    hpos <- hnewpos <- hdata2[, 2]
    hp <- hdata2[, 3]
    hname <- hdata2[, 4]
    hcolors <- hcontrol$colors
    hyoffs <- hcontrol$yoffs
    hboxed <- hcontrol$boxed
    hcex <- hcontrol$cex
    htablechr <- unique(hchr)
    hn.chr <- length(htablechr)
    hlabels <- unique(hname)
    htablechrname <- unique(data.frame(hchr, hname))
    if (is.null(hcolors)) 
      hcolors <- rep("red", length(hlabels))
    else hcolors <- hcontrol$colors
  }
  CMindex <- cumsum(allchr)
  for (i in 1:n.chr) {
    u <- CMindex[i]
    l <- CMindex[i] - allchr[i] + 1
    chr <- l:u
    if (usepos) 
      d <- diff(pos[chr])
    else d <- rep(1, allchr[i] - 1)
    newpos[chr] <- c(gap, d)
  }
  if (usepos)
    CM <- cumsum(as.numeric(newpos))+pos[1]
  else CM <- cumsum(as.numeric(newpos))
  args <- list(...)
  if ("ylim" %in% names(args)) 
    dp <- seq(args$ylim[1], args$ylim[2], length = sum(allchr))
  else dp <- seq(min(p), max(p), length = sum(allchr))
  if (logscale) 
    y <- -log(dp, base)
  else y <- dp
  y1 <- min(y)
  par(xaxt = "n", yaxt = "n")
  xy <- xy.coords(CM, y)
  plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE, ...)
  axis(1)
  axis(2)
  par(xaxt = "s", yaxt = "s")
  for (i in 1:n.chr) {
    u <- CMindex[i]
    l <- CMindex[i] - allchr[i] + 1
    chr <- l:u
    #cat("Plotting points ", l, "-", u, "\n")
    if (logscale) 
      y <- -log(p[chr], base)
    else y <- p[chr]
    col.chr <- colors[i]
    if (type == "l") 
      lines(CM[chr], y, col = col.chr, cex = pcex, ...)
    else points(CM[chr], y, col = col.chr, cex = pcex, ...)
    #text(ifelse(i == 1, CM[1], CM[l]), y1, pos = 1, offset = 1, 
    #     labels[i], srt = srt, ...)
    text(1/2*(CM[u]+CM[l]), y1, pos = 1, offset = 1, 
         labels[i], srt = srt, ...)
  }
  j <- 1
  for (i in 1:n.chr) {
    u <- CMindex[i]
    l <- CMindex[i] - allchr[i] + 1
    chr <- l:u
    if (logscale) 
      y <- -log(p[chr], base)
    else y <- p[chr]
    col.chr <- colors[i]
    if (!is.null(hcontrol$data)) {
      chrs <- htablechrname[tablechr[i] == htablechrname[, 
                                                         1], ]
      if (dim(chrs)[1] > 0) {
        hchrs <- as.character(chrs[, 2])
        for (k in 1:length(hchrs)) {
          hregion <- hpos[hchr == chrs[k, 1] & hname == 
                            hchrs[k]]
          hl <- chr[pos[chr] == hregion[1]]
          hu <- chr[pos[chr] == hregion[length(hregion)]]
          cat("  ... highlighting", hl, "-", hu, hchrs[k], 
              "\n")
          l1 <- hl - l + 1
          l2 <- hu - l + 1
          col.chr[l1:l2] <- hcolors[j]
          if (hboxed) {
            tg <- grid::textGrob(hchrs[k])
            rg <- grid::rectGrob(x = CM[chr][l1], y = max(y[l1:l2]) + 
                                   hyoffs, width = 1.1 * grid::grobWidth(tg), 
                                 height = 1.3 * grid::grobHeight(tg), gp = grid::gpar(col = "black", 
                                                                                      lwd = 2.5))
            boxedText <- grid::gTree(children = grid::gList(tg, 
                                                            rg))
            grid::grid.draw(boxedText)
          }
          else text(CM[chr][l1], max(y[l1:l2] + hyoffs), 
                    hchrs[k], cex = hcex)
          points(CM[l + (l1:l2)], y[l1:l2], col = col.chr[l1:l2], 
                 cex = pcex, ...)
          j <- j + 1
        }
      }
    }
  }
  if (!is.null(cutoffs)) 
    segments(0, cutoffs, n2 + gap * n.chr, cutoffs)
  if ("ylab" %in% names(args)) 
    mtext(args$ylab, 2, line = yline, las = 0)
  else mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
                                    sep = ""), "Observed value"), 2, line = yline, las = 0)
  if ("xlab" %in% names(args)) 
    xlabel <- args$xlab
  else xlabel <- ifelse(is.null(names(chr)), "Chromosome", 
                        names(chr))
  mtext(xlabel, 1, line = xline, las = 0)
  return(res=list(CM=CM,CMindex=CMindex)) # x coordinate and number of points in each chr
}

#data have 3 columns: chr, pos, pvalue
mhtplot2=function (data, control = mht.control(),  ...) 
{
  for (p in c("grid")) {
    if (length(grep(paste("^package:", p, "$", sep = ""), 
                    search())) == 0) {
      if (!require(p, quietly = TRUE, character.only = TRUE)) 
        warning(paste("mhtplot needs package `", p, "' to be fully functional; please install", 
                      sep = ""))
    }
  }
  #length of each chromosome
  chrs=c(1:22,"X","Y")
  dictfile="/fh/fast/dai_j/CancerGenomics/Tools/database/reference/compact/ucsc.hg19.compact.dict"
  chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
  chrlen=chrlen$V3
  chrlen=gsub("LN:","",chrlen,fixed=T)
  chrlen=as.numeric(chrlen)
  names(chrlen)=chrs
  
  data2 <- data[!apply(is.na(data), 1, any), ]
  idxnoNA <- !apply(is.na(data), 1, any)
  #the index of nonNA values. NAs will not be plotted
  idxnoNA <- which(idxnoNA==T)                  
  n2 <- dim(data2[1])
  chr <- as.character(data2[, 1])
  uniqchr=unique(chr)
  chrs1=chrs[chrs %in% uniqchr]
  idx=match(chrs1,uniqchr)
  uniqchr=uniqchr[idx]
  #if only draw chrom with data
  chrlen1=chrlen[names(chrlen) %in% chr]
  chrlen2=c(0,cumsum(chrlen1)) #start point of each chr
  names(chrlen2)=c(names(chrlen1),"nextchr")
  pos <- newpos <- data2[, 2]
  for (i in 1:length(uniqchr))
  {
    idx=which(data2[,1]==uniqchr[i])
    startp=chrlen2[which(names(chrlen2)==uniqchr[i])]
    # print(uniqchr[i])
    # print(range(idx))
    # print(startp)
    #the overall cooridinate start from 0
    pos[idx] = newpos[idx] =data2[idx,2]+startp
  }
  if (class(data2[1,3])=="character") data2[,3]=as.numeric(data2[,3])
  p <- data2[, 3]

  #this part has problem, the orders is chr1,chr10,chr2,...; should change to chr1,chr2,...,chr10,... 
  #allchr <- as.vector(tablechr)
  
  allchr=NULL
  for (i in 1:length(uniqchr))
  {
    allchr=c(allchr,sum(chr==uniqchr[i]))
  }
  n.chr <- length(allchr)
  type <- control$type
  #usepos <- control$usepos
  logscale <- control$logscale
  base <- control$base
  cutoffs <- control$cutoffs
  colors <- control$colors
  labels <- control$labels
  srt <- control$srt
  gap <- control$gap
  pcex <- control$cex
  yline <- control$yline
  xline <- control$xline
  colorlist <- colors()
  if (is.null(colors)) 
    colors <- sample(colorlist, n.chr)
  
  if (is.null(labels)) 
    labels <- uniqchr
  if (is.null(gap)) 
    gap <- 0
  
  CMindex <- cumsum(allchr)
  #the overall cooridinate start from 0
  CM=newpos
  args <- list(...)
  if ("ylim" %in% names(args)) 
    dp <- seq(args$ylim[1], args$ylim[2], length = sum(allchr))
  else dp <- seq(min(p), max(p), length = sum(allchr))
  if (logscale) 
    y <- -log(dp, base)
  else y <- dp
  y1 <- min(y)
  y2 = max(y)
  par(xaxt = "n", yaxt = "n")
#   xy <- xy.coords(CM, y)
#   plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE, ...)
  plot(x=chrlen2[c(1,length(chrlen2))],y=c(floor(y1),ceiling(y2)),type="n",ann = FALSE, axes = FALSE)
  #axis(1)
  #axis(2)
  par(xaxt = "s", yaxt = "s")
  for (i in 1:n.chr) {
    u <- CMindex[i]
    l <- CMindex[i] - allchr[i] + 1
    idx <- l:u
    #cat("Plotting points ", l, "-", u, "\n")
    if (logscale) 
      y <- -log(p[idx], base)
    else y <- p[idx]
    col.chr <- colors[i]
    if (type == "l") 
      lines(CM[idx], y, col = col.chr, cex = pcex, ...)
    else 
      {
        points(CM[idx], y, col = col.chr, cex = pcex, ...)
      }
    #text(ifelse(i == 1, CM[1], CM[l]), y1, pos = 1, offset = 1, 
    #     labels[i], srt = srt, ...)
    # if (i<=15 | (i %% 2==1 & i>15))
    # {
    #   text(1/2*(chrlen2[i]+chrlen2[i+1]), y1, pos = 1, offset = 1, labels[i], srt = srt, cex=1.1)
    # }
    text(1/2*(chrlen2[i]+chrlen2[i+1]), y1, pos = 1, offset = 1, labels[i], srt = srt, cex=0.8)
  }

  if (!is.null(cutoffs)) 
    segments(0, cutoffs, n2 + gap * n.chr, cutoffs)
  if ("ylab" %in% names(args)) 
    mtext(args$ylab, 2, line = yline, las = 0, cex=1.3)
  else mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
                                    sep = ""), "Observed value"), 2, line = yline, las = 0, ...)
  if ("xlab" %in% names(args)) 
    xlabel <- args$xlab
  else xlabel <- ifelse(is.null(names(chr)), "Chromosome", names(chr))
  mtext(xlabel, 1, line = xline, las = 0, cex=1.3)
  #CM:coordinate of points, for idxnoNA rows in the original data which don't have NAs
  return(res=list(CM=CM,CMindex=CMindex,chrlen2=chrlen2,y1=y1,y2=y2,idxnoNA=idxnoNA)) # x coordinate and number of points in each chr
}

#main function
#if only provide pvalues (named numeric) then only the plot will be drawn
#if pvalues is a dataframe, and contains chr start columns, then pvalues should be 3-column dataframe: in the order of chr,start,pvalue
#if provide fdrs (pvalues on its first column, fdr on its second columns), fdr threshold will also be drawn.
draw_manhattan=function(pvalues=NULL,fdrs=NULL,fdrthreshold=0.05,maxy=6,chrs=NULL,keepres=F,logscale=T,main=NULL,ylab=NULL)
{
  #if fdrs is a dataframe, the first row is pvalue, second row is fdr
  if (is.null(pvalues))
  {
    pvalues=fdrs[,1]
    if(class(pvalues)=="character") pvalues=as.numeric(pvalues)
    names(pvalues)=rownames(fdrs)
  }else
  {
    if (class(pvalues)=="data.frame")
    {
      if (! "chr" %in% colnames(pvalues) | ! "start" %in% colnames(pvalues)) #if the dataframe have start and chr columns, it already have coordinates
      {
        tmp=rownames(pvalues)
        pvalues=pvalues[,ncol(pvalues)] #the last column is pvalues
        if(class(pvalues)=="character") pvalues=as.numeric(pvalues)
        names(pvalues)=tmp
      }
    }
  }
  if (class(pvalues)=="numeric") #without gene coordinate, needs to figure out
  {
    if (is.null(chrs))
    {
      chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
      #chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
    }
    geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
    #form a dataframe with gene's position and p value
    df_pvalues=data.frame(gene=names(pvalues),pvalue=pvalues)
    genetable=merge(df_pvalues,geneposition,all.x=T)
    idxkeep=which(genetable$chr %in% chrs)
    genetable=genetable[idxkeep,]
    genetable=sortgenetable(genetable)
    genenames=genetable$gene
    #form the genetable for mahanttan plot
    genetable=data.frame(chr=genetable$chr,start=genetable$start,pvalue=genetable$pvalue,stringsAsFactors=F)
    rownames(genetable)=genenames
  }else
  {
    if (class(pvalues)=="data.frame")
    {
      warning("the columns should be chr, start and pvalue")
      genetable=pvalues
      genetable=sortgenetable(genetable)
    }
  }
  
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  par(las=1, xpd=TRUE,cex.axis=1.2,cex=1,cex.lab=1.5)
  ops<-mht.control(colors=color2,yline=1,xline=2,usepos=T,srt=0,logscale=logscale)
  if (is.null(ylab))
  {
    if (logscale==T)
    {
      ylab=expression(Observed~~-log[10](italic(p)))
    }else
    {
      ylab=expression(Observed~~italic(p))
    }
  }
  
  res=mhtplot2(genetable,ops,pch=10,bg=color2,ylab=ylab)
  #res=mhtplot2(genetable,ops,pch=10,bg=color2,cex.axis=1.2,cex.lab=1.2)
  #allpos keeps the overall cooridinate
  genetable=cbind(genetable,allpos=rep(NA,nrow(genetable)))
  genetable$allpos[res$idxnoNA]=res$CM
  #lines(c(0,max(res$CM)),c(0,0),lwd=2)
  if (logscale==T)
  {
    axis(2,at=c(floor(res$y1):(ceiling(res$y2))),pos=0, cex.axis=1.2)
  }else
  {
    axis(2,pos=0, cex.axis=1.2)
  }
  xat=c(res$chrlen2)
  axis(1,pos=floor(res$y1),labels=FALSE,tick=T,at=xat)
  if (! is.null(main))
    title(main=main)
  
  #to add fdr line
  #draw the fdr cutoff line
  if (!is.null(fdrs))
  {
    dictfile="/fh/fast/dai_j/CancerGenomics/Tools/database/reference/compact/ucsc.hg19.compact.dict"
    chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
    chrlen=chrlen$V3
    chrlen=gsub("LN:","",chrlen,fixed=T)
    chrlen=as.numeric(chrlen)
    names(chrlen)=c(1:22,"X","Y")
    #available chrs length
    chrlen1=chrlen[names(chrlen) %in% genetable$chr]
   
    #change fdrs,pvalues to vector
    fdrs1=fdrs
    if (class(fdrs)=="data.frame")
    {
      tmp=rownames(fdrs)
      fdrs1=fdrs[,ncol(fdrs)] #assume fdr is in the last column
      names(fdrs1)=tmp
    }
    pvalues1=pvalues
    if (class(pvalues)=="data.frame")
    {
      tmp=rownames(pvalues)
      pvalues1=pvalues[,ncol(pvalues)] #assume pvalue is in the last column
      names(pvalues1)=tmp
    }

    smallfdrs=fdrs1[fdrs1<=fdrthreshold]
    if (length(smallfdrs)>0) #if have fdr<fdrthreshold
    {
      idxs=which(fdrs1 %in% smallfdrs)
      pvaluethreshold=max(pvalues1[idxs])
      segments(0,-log10(pvaluethreshold),par('usr')[2],-log10(pvaluethreshold),col="red")
      #text(sum(chrlen1)/2,-log10(pvaluethreshold)+0.5,paste0("FDR=",fdrthreshold))
      #text(sum(chrlen1)*0.9,-log10(pvaluethreshold)+0.25,paste0("FDR=",fdrthreshold),cex=1.1)
      text(sum(chrlen1)*0.9,-log10(pvaluethreshold)+0.25,paste0("FDR=",round(max(smallfdrs),digits = 2)),cex=1.1)
    }else
    {
      warning("no genes with small fdr were found")
      #find the gene with min fdr
   
      minfdr=round(min(fdrs1,na.rm=T),digits = 2)
      idxs=which(fdrs1==minfdr)
      
      #use the one with max p-value
      thepvalue=max(pvalues1[idxs])
      #abline(h=-log10(thepvalue),col="red")
      segments(0,-log10(thepvalue),par('usr')[2],-log10(thepvalue),col="red")
      if (class(fdrs)=="data.frame")
      {
        #text(sum(chrlen1)/2,-log10(thepvalue)+0.5,paste0("FDR=",minfdr))
        text(sum(chrlen1)*0.9,-log10(thepvalue)+0.25,paste0("FDR=",minfdr),cex=1.1)
      }else
      {
        #text(sum(chrlen1)/2,-log10(thepvalue)+0.5,paste0("FDR=",minfdr))
        text(sum(chrlen1)*0.9,-log10(thepvalue)+0.25,paste0("FDR=",minfdr),cex=1.1)
      }
    }
  }
  par(cex.axis=1,cex=1,cex.lab=1)
  if (keepres==T) return(genetable)
}

##wrap function for pvalue,fdr,manhattan
pvalue_fdr_manhattan=function(data=data_copynumber,platform="copynumber",usetrain=T)
{
  if (usetrain==T)
  {
    pvalues_permutation_file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,"_train_permutation100.txt")
    fdr_platform=paste0(platform,"_train")
  }else
  {
    pvalues_permutation_file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,"_permutation100.txt")
    fdr_platform=platform
  }
  Sys.time()
  if (file.exists(pvalues_permutation_file))
  {
    pvalues=compute_pvalues(data=data,platform=platform,runpermutation=F,usetrain=usetrain)
  }else
  {
    pvalues=compute_pvalues(data=data,platform=platform,runpermutation=T,usetrain=usetrain)
  }
  Sys.time()
  
  pvalues_permutation=read.table(file=pvalues_permutation_file,header=F)
  Sys.time()
  fdrs=sapply_compute_qvalue_permutation(pvalues=pvalues,pvalues_permutation=pvalues_permutation,
                                                    numit=100,platform=fdr_platform)
  Sys.time()
  draw_manhattan(fdrs=fdrs,maxy=5,fdrthreshold=0.05)
}


#function to generate tables for copynumber genes with small fdr----------------------------------------
#data is used to compute mean values of two groups (amp,del)
#set includecor=T to compute correlation between copynumber and mrna 
smallfdrtable=function(fdrs=fdrs1000p_copynumber,data=data_copynumber,threshold=0.05,includecor=T)
{
  #gistic qvalues generated by firehose
  gisticscore=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/scores.gistic",header=T,sep="\t")
  gisticscore$Chromosome=gsub(23,"X",gisticscore$Chromosome)
  gr_gisticscore=GRanges(seqnames=gisticscore$Chromosome,IRanges(start=gisticscore$Start,end=gisticscore$End),type=gisticscore$Type,
                         gisticqvalue=10^-gisticscore$X.log10.q.value.)
  
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
  geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
  if (class(fdrs)=="data.frame")
  {
    #first column of fdrs is pvalues
    df_fdrs=data.frame(gene=rownames(fdrs),pvalues=fdrs[,1],fdr=fdrs[,2])
  }else
  {
    #form a dataframe with gene's position and p value
    df_fdrs=data.frame(gene=names(fdrs),fdr=fdrs)
  }
  
  genetable=merge(df_fdrs,geneposition)
  idxkeep=which(genetable$chr %in% chrs)
  genetable=genetable[idxkeep,]
  genetable=sortgenetable(genetable)
  genenames=genetable$gene

  idx_keep=! is.na(genetable$fdr) & genetable$fdr<=threshold
  if (sum(idx_keep)==0)
  {
    warning("it has no genes with small fdr")
    #if pvalues provided, use pvalues
    if ("pvalues" %in% colnames(genetable))
    {
      tmp=order(genetable$pvalues)[100]
      idx_keep=which(genetable$pvalues<=genetable$pvalues[tmp])
    }else
    {
      #select about top 100 genes with smallest q values
      tmp=order(genetable$fdr)[100]
      idx_keep=which(genetable$fdr<=genetable$fdr[tmp])
    }
  }
  selecttable=genetable[idx_keep,]
  selecttable=sortgenetable(selecttable)
  #extract q values from gistic results
  delgisticqvalue=NULL
  ampgisticqvalue=NULL
  for (i in 1:nrow(selecttable))
  {
    gr_selecttable=GRanges(seqnames=selecttable$chr[i],ranges=IRanges(start=selecttable$start[i],end=selecttable$end[i]),qvalue=selecttable$qvalues[i])
    olap=subsetByOverlaps(gr_gisticscore,gr_selecttable)
    if (length(olap)>0)
    {
      idx=rep(F,length(olap))
      idx_amp=which(mcols(olap)$type=="Amp")
      idx[idx_amp]=T
      olap_amp=olap[idx,]
      olap_del=olap[!idx,]
      if (length(olap_amp)>0)
      {
        ampgisticqvalue=c(ampgisticqvalue,min(mcols(olap_amp)$gisticqvalue))
      }else
      {
        ampgisticqvalue=c(ampgisticqvalue,NA)
      }
      if (length(olap_del)>0)
      {
        delgisticqvalue=c(delgisticqvalue,min(mcols(olap_del)$gisticqvalue))
      }else
      {
        delgisticqvalue=c(delgisticqvalue,NA)
      }
    }else
    {
      ampgisticqvalue=c(ampgisticqvalue,NA)
      delgisticqvalue=c(delgisticqvalue,NA)
    }
  }
  selecttable=cbind(selecttable,ampgisticqvalue,delgisticqvalue)
#   #write the genes to a file
#   filename1=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_genes_smallqvalues.txt")
#   fcon=file(filename1,"w")
#   for (i in 1:nrow(selecttable))
#   {
#     writeLines(rownames(selecttable[i,]),fcon)
#   }
#   close(fcon)

  #update stage and residual and to form a dataframe
  alldata=formwholedataframe(data)
#   #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
#   alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
#   alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
#   alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
#   alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  idx_y=which(colnames(alldata)=="platinumclass")
  #compute the correlation between two consecutive genes
  bet_cor=NA
  for (i in 2:nrow(selecttable))
  {
    idx1=which(colnames(alldata)==selecttable$gene[i-1])
    idx2=which(colnames(alldata)==selecttable$gene[i])
    if (length(idx1)>0 & length(idx2)>0)
    {
      tmp=cor(alldata[,idx1],alldata[,idx2],use="complete")
    }else
    {
      tmp=NA
    }
    bet_cor=c(bet_cor,tmp)
  }
  selecttable=cbind(selecttable,betweencorrelation=bet_cor)
  #compute the mean in two classes
  resistantdata=alldata[tolower(alldata[,"platinumclass"])=="resistant",(idx_y+1):ncol(alldata)]
  sensitivedata=alldata[tolower(alldata[,"platinumclass"])=="sensitive",(idx_y+1):ncol(alldata)]
  mean_res=data.frame(mean_resi=rep(NA,nrow(selecttable)),mean_sens=rep(NA,nrow(selecttable)),pvalue=rep(NA,nrow(selecttable)))
  for (i in 1:nrow(selecttable))
  {
    gene=selecttable$gene[i]
    idx=which(colnames(resistantdata)==gene)
    #if (length(idx)>0)
    #{
      mean_res[i,1]=mean(resistantdata[,idx])
      mean_res[i,2]=mean(sensitivedata[,idx])
      mean_res[i,3]=t.test(resistantdata[,idx],sensitivedata[,idx])$p.value
      #info=paste0(selecttable$chr[i],":",selecttable$start[i],"-",selecttable$end[i])
      #hist(resistantdata[,idx],xlab="",main=paste0("resitant,",gene," ",info),breaks=20)
      #hist(sensitivedata[,idx],xlab="",main=paste0("sensitive,",gene," ",info),breaks=20)
    #}
  }
  selecttable=cbind(selecttable,mean_res)
  if (includecor==T)
  {
    load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")
    cor_mrna_copynumber=data.frame(gene=selecttable$gene,correlation=rep(NA,nrow(selecttable)))
    for (i in 1:nrow(selecttable))
    {
      gene=selecttable$gene[i]
      idx=which(rownames(copynumber_comcm)==gene)
      if (length(idx)>0)
      {
        cor_mrna_copynumber[i,2]=cor(as.numeric(mrna_comcm[idx,]),as.numeric(copynumber_comcm[idx,]),use="complete")
      }
    }
    selecttable=cbind(selecttable,cor_mrna_copynumber=cor_mrna_copynumber[,2])
  }
#   filename3=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneswithsmallqvaluestable_",platform,".txt")
#   write.table(selecttable,file=filename3,row.names=F,col.names=T,sep="\t",quote=F)
  return(selecttable)
}

##other functions-----------------------------------------------------------------
#map genes1 with genes2 considering alias genes name. symbol: current gene symbol name
mapaliasgenes=function(genes1,genes2,tolowerflag=F)
{
  res=data.frame(matrix(NA,ncol=4,nrow=length(genes1)))
  #idx2 is the index of found match in genes2
  colnames(res)=c("genes1","genes2","symbol","idx2")
  res[,1]=genes1
  matchidx=which(genes1 %in% genes2)
  if (length(matchidx)>0)
  {
    for (i in matchidx)
    {
      res[i,c(2,3)]=genes1[i]
      res[i,4]=which(genes2==genes1[i])[1]
    }
  }
  #search alias gene names through aliasgenetable
  genes2_aliasgenetable=data.frame(matrix(NA,ncol=2,nrow=length(genes2)))
  colnames(genes2_aliasgenetable)=c("gene","symbol")
  genes2_aliasgenetable[,1]=genes2_aliasgenetable[,2]=genes2
  idx_match=which(genes2 %in% aliasgenetable$gene)
  for (i in idx_match)
  {
    idx=which(aliasgenetable$gene==genes2[i])
    genes2_aliasgenetable[i,2]=aliasgenetable$symbol[idx[1]]
  }
  tofindidx=(1:length(genes1))[-matchidx]
  if (length(tofindidx)>0)
  {
    for (i in tofindidx)
    {
      idx=which(aliasgenetable$gene==genes1[i])
      if (length(idx)>0)
      {
        genes1_symbol=aliasgenetable$symbol[idx[1]]
        idx1=which(genes2_aliasgenetable$symbol==genes1_symbol)
        if (length(idx1)==0 & tolowerflag==T)
        {
          idx1=which(tolower(genes2_aliasgenetable$symbol)==tolower(genes1_symbol))
        }
        if (length(idx1)>0)
        {
          print(paste0("found genes through alias name for",genes1[i]))
          res[i,2]=genes2[idx1[1]]
          res[i,3]=genes1_symbol
          res[i,4]=idx1
        }
      }
    }
  }
  return(res)
}


#form copynumber data at gene level based on segment-------------------
#only for hg19
readcnasegs1=function(segments,opt=1)
{
  library(CNTools)
  colnames(segments)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  dim(segments)
  length(unique(segments$ID))
  segments=segments[!is.na(segments$seg.mean),]
  cnseg=CNSeg(segments)
  if (opt==1)
  {
    #     hg19_=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/copynumber_geneposition_biomart.txt",header=T,sep="\t")
    #     geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
    #     hg19=merge(hg19_,geneposition,by="gene")
    #     hg19=hg19[,c(1,5,6,7)]
    #use the one from firehose
    hg19=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
  }
  if (opt==2) #use the genesets defined by CNTools
  {
    load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")
    hg19=geneInfo_hg19[,c(4,1,2,3)]
  }
  colnames(hg19)=c("gene","chr","start","end")
  hg19=sortgenetable(hg19)
  
  hg19$chr=gsub("23","X",hg19$chr)
  hg19$chr=gsub("24","Y",hg19$chr)
  idxdup=duplicated(hg19$gene)
  hg19=hg19[! idxdup,]
  #change the column name of chr
  colnames(hg19)[2]="chrom"
  nstart=5
  genecn_max<- getRS(cnseg, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "max")
  genecn_max=rs(genecn_max)
  genecn_min<- getRS(cnseg, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "min")
  genecn_min=rs(genecn_min)
  genecn=genecn_max
  
  #pick the extreme one
  for (i in nstart:ncol(genecn))
  {
    pickmin=which(abs(genecn_min[,i])-abs(genecn_max[,i])>0)
    if (length(pickmin)>0) genecn[pickmin,i]=genecn_min[pickmin,i]
  }
  
  rownames(genecn)=genecn$gene
  genecn=genecn[,nstart:ncol(genecn)]
  
  #remove all 0 genes
  idxkeep=rep(T,nrow(genecn))
  idxkeep=sapply(1:nrow(genecn),function(i){
    res=F
    num0=sum(genecn[i,]==0)
    if (num0<ncol(genecn))
    {
      res=T
    }
    return(res)
  })
  genecn=genecn[idxkeep,]
  return(genecn)
}

#form TCGA plastimum definition
load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_firehose.RData")
tcga_platinumclass=data.frame(sample=rownames(data_copynumber),platinumclass=data_copynumber[,"platinumclass"])
tcga_platinumclass$sample=as.character(tcga_platinumclass$sample)

#compute 2df-pvalue (tcn and dh)------
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
two_df_pvalues_data=function(data1=data_copynumber,data2=data_copynumber_allele)
{
  #make sure the two datasets have the same order in row and columns
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
  pvalues=sapply(1:ncol(data1),function(i){
    x=data1[,i]
    y=data2[,i]
    res=two_df_pvalues(x=x,y=y,z=z)
  })
  names(pvalues)=colnames(data1)
  ##for some cases when x or y are mostly 0, which causes the problem of p=0
  idx_pvalue0=which(pvalues==0)
  if (length(idx_pvalue0)>0) pvalues[idx_pvalue0]=NA
  
  return(pvalues)
}

#calculate correlation between tcn and dh
cor_tcn_dh_data=function(data1=data_copynumber,data2=data_copynumber_allele)
{
  #make sure the two datasets have the same order in row and columns
  commongenes=intersect(colnames(data1)[2:ncol(data1)],colnames(data2)[2:ncol(data2)])
  commonsamples=intersect(rownames(data1)[1:nrow(data1)],rownames(data2)[1:nrow(data2)])
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data1)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data1)==mysample)
  })
  data1=data1[idx_row,idx_col]
  idx_col=sapply(commongenes, function(gene){
    res=which(colnames(data2)==gene)
  })
  idx_row=sapply(commonsamples, function(mysample){
    res=which(rownames(data2)==mysample)
  })
  data2=data2[idx_row,idx_col]
  corvar=sapply(1:ncol(data1),function(i){
    x=data1[,i]
    y=data2[,i]
    res=cor(x,y,use="complete")
  })
  names(corvar)=colnames(data1)
  return(corvar)
}
