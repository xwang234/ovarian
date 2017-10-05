#!/usr/bin/env Rscript
library(glmnet)
library(gap) #Manhattan plot
library(GenomicRanges) #find overlaps
library(ROCR)
library(pROC)
#include functions
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
class_copynumber=classify_platinum(data=data_rawmutation_hrgenes,platform="mutation_hrgenes",selgenes=NULL,includeclinical=T)
data=data_rawmutation_hrgenes
alldata=rbind(data$data1,data$data2)
#add indicator of train/test as the last column
train=rep(FALSE,nrow(alldata))
train[1:nrow(data$data1)]=TRUE
alldata=cbind(alldata,train=train)
alldata=formwholedataframe(alldata)
alldata$platinumclass=as.factor(alldata$platinumclass)
alldata$train=as.logical(alldata$train)
traindata=alldata[alldata[,"train"]==T,]
fm=glm(as.formula("platinumclass~hrgenes"),data=alldata[alldata[,"train"]==T,],family="binomial")
summary(fm)
pfit_train=predict(fm,type="response")
pfit_test=predict(fm,newdata=alldata[alldata[,"train"]==F,],type="response")
plotroc2(pfit_train,alldata[alldata[,"train"]==T,"platinumclass"],pfit_test,alldata[alldata[,"train"]==F,"platinumclass"],main=paste0("ROC"))

idx=traindata$hrgenes>0
table(pfit_train[idx])
table(pfit_train[!idx])
fm1=glm(as.formula("platinumclass~hrgenes"),data=alldata,family="binomial")
summary(fm1)

#get genes from glmnet on copynumber
class_copynumber=classify_platinum(data=data_copynumber,platform="copynumber",selgenes=NULL,includeclinical=T)

#add hrgenes to the first column
addhrgenes=function(data=data_copynumber,adddata=data_rawmutation_hrgenes)
{
  adddata=rbind(adddata$data1,adddata$data2)
  data1=data.frame(matrix(NA,nrow=nrow(data$data1),ncol=ncol(data$data1)+ncol(adddata)-13))
  #clinical info
  data1[,1:13]=data$data1[,1:13]
  #initialize as 0 matrix, missing value will be 0
  tmp=data.frame(matrix(0,nrow=nrow(data$data1),ncol=ncol(adddata)-13))
  colnames(tmp)=colnames(adddata)[14:ncol(adddata)]
  for (i in 1:nrow(data$data1))
  {
    idx=which(rownames(adddata)==rownames(data$data1)[i])
    if (length(idx)>0)
    {
      tmp[i,]=adddata[idx,14:ncol(adddata)]
    }
  }
  data1[,14:ncol(adddata)]=tmp
  data1[,(ncol(adddata)+1):ncol(data1)]=data$data1[,14:ncol(data$data1)]
  
  data2=data.frame(matrix(NA,nrow=nrow(data$data2),ncol=ncol(data$data2)+ncol(adddata)-13))
  #clinical info
  data2[,1:13]=data$data2[,1:13]
  #initialize as 0 matrix, missing value will be 0
  tmp=data.frame(matrix(0,nrow=nrow(data$data2),ncol=ncol(adddata)-13))
  colnames(tmp)=colnames(adddata)[14:ncol(adddata)]
  for (i in 1:nrow(data$data2))
  {
    idx=which(rownames(adddata)==rownames(data$data2)[i])
    if (length(idx)>0)
    {
      tmp[i,]=adddata[idx,14:ncol(adddata)]
    }
  }
  data2[,14:ncol(adddata)]=tmp
  data2[,(ncol(adddata)+1):ncol(data2)]=data$data2[,14:ncol(data$data2)]
  colnames(data1)=colnames(data2)=c(colnames(data$data1)[1:13],colnames(adddata)[14:ncol(adddata)],colnames(data$data1)[14:ncol(data$data1)])
  result=list(data1=data1,data2=data2)
}
data_copynumber_hrgenes=addhrgenes()

data_class_copynumber_hrgenes=addhrgenes(data=data_class_copynumber)
#when combine hrgnes, it was not selected finally
class_copynumber_hrgenes=classify_platinum(data=data_copynumber_hrgenes,platform="copynumber_hrgenes",selgenes=NULL,includeclinical=T)

#check the categorical copynumber data:
load("/fh/fast/dai_j/CancerGenomics/data/platinum_classification_copynumber_category.RData")
alldata=rbind(data_class_copynumber$data1,data_class_copynumber$data2)
trainflag=rep(T,nrow(alldata))
set.seed(222)
testidx=sample(1:nrow(alldata),ceiling(nrow(alldata)/3))
trainflag[testidx]=F
data_class_copynumber_rand=list(data1=alldata[trainflag,],data2=alldata[!trainflag,])
alldata1=rbind(data_class_copynumber_hrgenes$data1,data_class_copynumber_hrgenes$data2)
data_class_copynumber_hrgenes_rand=list(data1=alldata1[trainflag,],data2=alldata1[!trainflag,])
class_cat_copynumber=classify_platinum(data=data_class_copynumber_rand,platform="copynumber_cat",selgenes=NULL,includeclinical=T)
class_cat_copynumber1=classify_platinum(data=data_class_copynumber_rand,platform="copynumber_cat",selgenes=NULL,includeclinical=F)
class_cat_copynumber_hrgenes=classify_platinum(data=data_class_copynumber_hrgenes,platform="copynumber_cat_hrgenes",selgenes=NULL,includeclinical=T)

pvalues_train=sapply(14:ncol(data_class_copynumber_rand$data1),function(i){
  tmp=glm(as.formula(paste0(colnames(data_class_copynumber_rand$data1)[11],"~",colnames(data_class_copynumber_rand$data1)[i])),data=data_class_copynumber_rand$data1,family="binomial")
  tmp1=coef(summary(tmp))
  res=NA
  names(res)=colnames(data_class_copynumber_rand$data1)[i]
  if (nrow(tmp1)>1)
  {
    res=tmp1[2,4]
  }
  return(res)
})


selgenes=class_copynumber$coeff
selgenes=class_cat_copynumber1$coeff

selgenes=class_cat_copynumber_hrgenes$coeff
selgenes=selgenes[!names(selgenes) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease","age")]
selgenes=names(selgenes)


drawroc_onselectedgenes(selgenes = selgenes)
selgenes=c(selgenes,"hrgenes")
#return traindata,testdata:
data=drawroc_onselectedgenes(selgenes = selgenes,data=data_class_copynumber_hrgenes_rand)


#investigate selected genes using training data:
traindata=data$traindata
testdata=data$testdata

colnames(traindata)=colnames(testdata)=gsub("|","___",colnames(traindata),fixed=T)
selgenes=gsub("|","___",selgenes,fixed=T)
traindata=traindata[,c(11,14:ncol(traindata))]
# pvalues=sapply(2:ncol(traindata),function(i){
#   tmp=glm(as.formula(paste0(colnames(traindata)[1],"~",colnames(traindata)[i])),data=traindata,family="binomial")
#   tmp1=coef(summary(tmp))
#   res=NA
#   if (nrow(tmp1)>1)
#   {
#     res=tmp1[2,4]
#   }
#   return(res)
# })
# hist(traindata[traindata$platinumclass=="Sensitive",3])
# hist(traindata[traindata$platinumclass=="Resistant",3])

#convert to categorical data
traindata1=traindata

# #if the data is not categorical
# traindata1[,2:ncol(traindata1)]=0
# for (i in 2:ncol(traindata1))
# {
#   idx=traindata[,i]>=0.5
#   traindata1[idx,i]=1
#   idx=traindata[,i]<= -0.5
#   traindata1[idx,i]=-1
# }
# 
# hist(traindata1[traindata1$platinumclass=="Sensitive",6])
# hist(traindata1[traindata1$platinumclass=="Resistant",6])
#check deletions
pvalues_del=sapply(3:ncol(traindata1),function(i){
  M=table(traindata1$platinumclass,traindata1[,i]<0)
  if (ncol(M)>1)
  {
    res=fisher.test(M)$p.value
  }else
  {
    res=NA
  }
  names(res)=colnames(traindata1)[i]
  return(res)
})
pvalues_amp=sapply(3:ncol(traindata1),function(i){
  M=table(traindata1$platinumclass,traindata1[,i]>0)
  if (ncol(M)>1)
  {
    res=fisher.test(M)$p.value
  }else
  {
    res=NA
  }
  names(res)=colnames(traindata1)[i]
  return(res)
})

#compute probability of likelihood
getlikelihood=function(pvalues=pvalues_del,opt="del",data=traindata1)
{
  idx=which(pvalues<=0.1)
  res=data.frame(matrix(NA,nrow=length(idx),ncol=5))
  rownames(res)=names(pvalues)[idx]
  colnames(res)=c("Sensitive","Resistant","Type","pvalue","proportion")
  for (i in 1:length(idx))
  {
    idx1=which(colnames(data)==names(pvalues)[idx[i]])
    if (opt=="del")
    {
      res[i,"Sensitive"]=sum(data[data$platinumclass=="Sensitive",idx1]<0)/sum(data$platinumclass=="Sensitive")
      res[i,"Resistant"]=sum(data[data$platinumclass=="Resistant",idx1]<0)/sum(data$platinumclass=="Resistant")
      res[i,"Type"]="del"
      res[i,"pvalue"]=pvalues[idx[i]]
      res[i,"proportion"]=sum(data[,idx1]<0)/nrow(data)
    }else
    {
      res[i,"Sensitive"]=sum(data[data$platinumclass=="Sensitive",idx1]>0)/sum(data$platinumclass=="Sensitive")
      res[i,"Resistant"]=sum(data[data$platinumclass=="Resistant",idx1]>0)/sum(data$platinumclass=="Resistant")
      res[i,"Type"]="amp"
      res[i,"pvalue"]=pvalues[idx[i]]
      res[i,"proportion"]=sum(data[,idx1]>0)/nrow(data)
    }
  }
  return(res)
}
likelihood_del=getlikelihood(pvalues=pvalues_del,opt="del",data=traindata1)
likelihood_amp=getlikelihood(pvalues=pvalues_amp,opt="amp",data=traindata1)
#for hrgenes
likelihood_hrgenes=c(sum(traindata1[traindata1$platinumclass=="Sensitive",2]>0)/sum(traindata1$platinumclass=="Sensitive"),
                     sum(traindata1[traindata1$platinumclass=="Resistant",2]>0)/sum(traindata1$platinumclass=="Resistant"),
                     "hrgenes","0.07",sum(traindata1$hrgenes>0)/nrow(traindata1))
likelihoodmatrix=rbind(likelihood_hrgenes,likelihood_del,likelihood_amp)
likelihoodmatrix$Sensitive=as.numeric(likelihoodmatrix$Sensitive)
likelihoodmatrix$Resistant=as.numeric(likelihoodmatrix$Resistant)
likelihoodmatrix$proportion=as.numeric(likelihoodmatrix$proportion)
rownames(likelihoodmatrix)[1]="hrgenes"
prior_sensitive=sum(traindata1$platinumclass=="Sensitive")/nrow(traindata1)
prior_resistant=1-prior_sensitive

#for testdata

testdata=testdata[,c(11,14:ncol(testdata))]
testdata1=testdata

# testdata1[,2:ncol(testdata1)]=0
# for (i in 2:ncol(testdata1))
# {
#   idx=testdata[,i]>=0.5
#   testdata1[idx,i]=1
#   idx=testdata[,i]<= -0.5
#   testdata1[idx,i]=-1
# }

bayesclass=function(newdata=testdata1,usedlikelihoodmatrix=likelihoodmatrix)
{
  x=newdata[,3:ncol(newdata)]
  y=newdata$platinumclass
  
  res=data.frame(matrix(NA,ncol=3,nrow=nrow(newdata)))
  colnames(res)=c("Sensitive","Resistant","Class")
  for (i in 1:nrow(newdata))
  {
    prob_sensitive=prior_sensitive
    prob_resistant=prior_resistant
    #check hrgenes
    if ("hrgenes" %in% rownames(usedlikelihoodmatrix))
    {
      if (newdata[i,"hrgenes"]>0)
      {
        prob_sensitive=prob_sensitive*usedlikelihoodmatrix$Sensitive[1]
        prob_resistant=prob_resistant*usedlikelihoodmatrix$Resistant[1]
      }
    }
    
    #check dels:
    x1=as.numeric(x[i,])
    idx_del=which(x1<0)
    names(idx_del)=colnames(x)[idx_del]
    selectmatrix=data.frame(matrix(NA,nrow=0,ncol=ncol(usedlikelihoodmatrix)))
    colnames(selectmatrix)=colnames(usedlikelihoodmatrix)
    if (length(idx_del)>0)
    {
      idx2=which(rownames(usedlikelihoodmatrix) %in% names(idx_del) & usedlikelihoodmatrix$Type=="del")
      if (length(idx2)>0)
      {
        delmatrix=usedlikelihoodmatrix[idx2,]
        selectmatrix=rbind(selectmatrix,delmatrix)
      }
    }
    #check amps:
    idx_amp=which(x1>0)
    names(idx_amp)=colnames(x)[idx_amp]
    if (length(idx_amp)>0)
    {
      idx2=which(rownames(usedlikelihoodmatrix) %in% names(idx_amp) & usedlikelihoodmatrix$Type=="amp")
      if (length(idx2)>0)
      {
        ampmatrix=usedlikelihoodmatrix[idx2,]
        selectmatrix=rbind(selectmatrix,ampmatrix)
      }
    }
    if (nrow(selectmatrix)>0)
    {
      for (j in 1:nrow(selectmatrix))
      {
        prob_sensitive=prob_sensitive*selectmatrix$Sensitive[j]
        prob_resistant=prob_resistant*selectmatrix$Resistant[j]
      }
    }
    prob_evidence=(prob_sensitive+prob_resistant)
    res$Sensitive[i]=prob_sensitive/prob_evidence
    res$Resistant[i]=prob_resistant/prob_evidence
    res$Class[i]=ifelse((prob_sensitive > prob_resistant),"Sensitive","Resistant")
  }
  auc=calauc(res$Sensitive,y)
  result=list(mat=res,auc=auc)
  return(result)
}
likelihoodmatrix1=likelihoodmatrix[which(likelihoodmatrix$Resistant>likelihoodmatrix$Sensitive),]
train=bayesclass(newdata = traindata1,usedlikelihoodmatrix = likelihoodmatrix)$mat
train=cbind(train,traindata1$platinumclass)
table(train[,3],train[,4])

test=bayesclass(newdata = testdata1,usedlikelihoodmatrix = likelihoodmatrix)$mat
test=cbind(test,testdata1$platinumclass)
table(test[,3],test[,4])

pfit_train=bayesclass(newdata=traindata1,usedlikelihoodmatrix = likelihoodmatrix)$mat[,2]
pfit_test=bayesclass(newdata=testdata1,usedlikelihoodmatrix = likelihoodmatrix)$mat[,2]
plotroc2(pfit_train,traindata1$platinumclass,pfit_test,testdata1$platinumclass,plotflag=1,main="")

selectvariables=function(newdata=traindata1,usedlikelihoodmatrix=likelihoodmatrix,testdata=NULL)
{

  runflag=T
  count=1
  res=data.frame(matrix(NA,nrow=0,ncol=3))
  colnames(res)=c("maxauc","rmidx","testacu")
  while(runflag==T)
  {
    cat(count,"..")
    
    if (count==1)
    {
      auc=bayesclass(newdata,usedlikelihoodmatrix = usedlikelihoodmatrix)$auc
    }
    count=count+1
    if (nrow(usedlikelihoodmatrix)==1)
    {
      runflag=F
    }
    aucs=rep(0,nrow(usedlikelihoodmatrix))
    for (i in 1:nrow(usedlikelihoodmatrix))
    {
      aucs[i]=bayesclass(newdata,usedlikelihoodmatrix = usedlikelihoodmatrix[-i,])$auc
    }
    maxauc=max(aucs[1:length(aucs)],na.rm=T)
    cat(maxauc,"..")
    #remove it will result in max auc
    idx_minauc=which.max(aucs)
    if (!is.null(testdata))
    {
     testauc=bayesclass(newdata=testdata,usedlikelihoodmatrix=usedlikelihoodmatrix[-idx_minauc,])$auc
     cat(testauc,"\n")
    }
    #the reduced model is better
    if (maxauc>=auc)
    {
      rmidx=which(rownames(likelihoodmatrix)==rownames(usedlikelihoodmatrix)[idx_minauc])
      res=rbind(res,c(maxauc,rmidx,testauc))
      usedlikelihoodmatrix=usedlikelihoodmatrix[-idx_minauc,]

      # auc=bayesclass(newdata,usedlikelihoodmatrix=usedlikelihoodmatrix)$auc
    }else
    {
      runflag=F
    }
  }
  finallikelihoodmatrix=usedlikelihoodmatrix
  result=list(finallikelihoodmatrix=finallikelihoodmatrix,res=res)
}

pfit_train=bayesclass(newdata=traindata1,usedlikelihoodmatrix = usedlikelihoodmatrix)$mat[,2]
pfit_test=bayesclass(newdata=testdata1,usedlikelihoodmatrix = usedlikelihoodmatrix)$mat[,2]
plotroc2(pfit_train,traindata1$platinumclass,pfit_test,testdata1$platinumclass,plotflag=1,main="")

# #compute all pvalues:
# traindata=data_class_copynumber_hrgenes$data1
# traindata=traindata[,c(11,14:ncol(traindata))]
# colnames(traindata)=gsub("|","___",colnames(traindata),fixed=T)
# colnames(traindata)=gsub("-","____",colnames(traindata),fixed=T)
# colnames(traindata)[3:ncol(traindata)]=paste0("addxxxx",colnames(traindata)[3:ncol(traindata)])
# 
# pvalues_del=sapply(3:ncol(traindata),function(i){
#   M=table(traindata$platinumclass,traindata[,i]<0)
#   if (ncol(M)>1)
#   {
#     res=fisher.test(M)$p.value
#   }else
#   {
#     res=NA
#   }
#   names(res)=colnames(traindata)[i]
#   return(res)
# })
# pvalues_amp=sapply(3:ncol(traindata),function(i){
#   M=table(traindata$platinumclass,traindata[,i]>0)
#   if (ncol(M)>1)
#   {
#     res=fisher.test(M)$p.value
#   }else
#   {
#     res=NA
#   }
#   names(res)=colnames(traindata)[i]
#   return(res)
# })
# 
# names(pvalues_del)=gsub("addxxxx","",names(pvalues_del),fixed=T)
# names(pvalues_amp)=gsub("addxxxx","",names(pvalues_amp),fixed=T)
# colnames(traindata)=gsub("addxxxx","",colnames(traindata),fixed=T)
# sum(pvalues_amp<0.1,na.rm=T)
# sum(pvalues_del<0.1,na.rm=T)
# likelihood_del=getlikelihood(pvalues=pvalues_del,opt="del",data=traindata)
# likelihood_amp=getlikelihood(pvalues=pvalues_amp,opt="amp",data=traindata)
# #for hrgenes
# likelihood_hrgenes=c(sum(traindata1[traindata1$platinumclass=="Sensitive",2]>0)/sum(traindata1$platinumclass=="Sensitive"),
#                      sum(traindata1[traindata1$platinumclass=="Resistant",2]>0)/sum(traindata1$platinumclass=="Resistant"),
#                      "hrgenes","0.07",sum(traindata1$hrgenes>0)/nrow(traindata1))
# likelihoodmatrix=rbind(likelihood_hrgenes,likelihood_del,likelihood_amp)
# likelihoodmatrix$Sensitive=as.numeric(likelihoodmatrix$Sensitive)
# likelihoodmatrix$Resistant=as.numeric(likelihoodmatrix$Resistant)
# likelihoodmatrix$proportion=as.numeric(likelihoodmatrix$proportion)
# rownames(likelihoodmatrix)[1]="hrgenes"
# prior_sensitive=sum(traindata1$platinumclass=="Sensitive")/nrow(traindata1)
# prior_resistant=1-prior_sensitive

#boosting:

colnames(data_class_copynumber$data1)=colnames(data_class_copynumber$data2)=gsub("|","__",colnames(data_class_copynumber$data1),fixed=T)
alldada=rbind(data_class_copynumber$data1,data_class_copynumber$data2)
trainflag=rep(T,nrow(alldata))
set.seed(222)
trainflag[sample(1:nrow(alldata),ceiling(1/3*nrow(alldata)))]=F
traindata=alldada[trainflag,]
testdata=alldata[!trainflag,]

selgenes=NULL
for (i in 1:20)
{
  cat(i,"...")
  trainflag1=rep(T,nrow(traindata))
  set.seed(i)
  trainflag1[sample(1:nrow(traindata),ceiling(0.5*nrow(traindata)))]=F
  data_copynumber_rand1=list(data1=traindata[trainflag1,],data2=traindata[!trainflag1,])
  class_copynumber_rand1=classify_platinum(data=data_copynumber_rand1,platform="copynumber_rand",selgenes=NULL,includeclinical=F)
  selgenes1=class_copynumber_rand$coeff
  selgenes1=selgenes1[!names(selgenes1) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease","age")]
  selgenes1=names(selgenes1)
  selgenes=c(selgenes,selgenes1)
  selgenes=unique(selgenes)
}


library(adabag)
adaboost<-boosting(as.formula(paste0("platinumclass~",paste0(selgenes,collapse="+"))), data=traindata, boos=TRUE, mfinal=20,coeflearn='Zhu')
pfit_train=predict(adaboost,traindata)$prob[,1]
pfit_test=predict(adaboost,testdata)$prob[,1]
plotroc2(pfit_train,traindata[,"platinumclass"],pfit_test,testdata[,"platinumclass"],main=paste0("Adaba:ROC"))
t2<-adaboost$trees[[2]]
library(tree)
plot(t2)
text(t2,pretty=0)
adaboost$importance
