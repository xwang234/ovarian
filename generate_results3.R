#!/usr/bin/env Rscript
library(glmnet)
library(gap) #Manhattan plot
library(GenomicRanges) #find overlaps
library(ROCR)
library(pROC)
#include functions
source(file="./functions.R")

load("../data/platinum_classificationdata_stringent.RData")
load("../data/platinum_4classificationdata_stringent.RData")

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

data_class_copynumber_hrgenes=addhrgenes(data=data_class_copynumber)
#when combine hrgnes, it was not selected finally
class_copynumber_hrgenes=classify_platinum(data=data_copynumber_hrgenes,platform="copynumber_hrgenes",selgenes=NULL,includeclinical=T)

#check the categorical copynumber data:
load("./platinum_classification_copynumber_category.RData")
alldata=rbind(data_class_copynumber$data1,data_class_copynumber$data2)
trainflag=rep(T,nrow(alldata))
set.seed(222)
testidx=sample(1:nrow(alldata),ceiling(nrow(alldata)/3))
trainflag[testidx]=F
alltraindata=alldata[trainflag,]
testdata=alldata[!trainflag,]

compute2pvalues=function(data=alltraindata,it=1,opt="del")
{
  res=data.frame(matrix(NA,ncol=2,nrow=ncol(data)-13))
  colnames(res)=c("pvalue1","pvalue2")
  rownames(res)=colnames(data)[14:ncol(data)]
  
  trainflag1=rep(T,nrow(data))
  set.seed(it)
  trainflag1[sample(1:length(trainflag1),ceiling(1/2*length(trainflag1)))]=F              
  data1=data[trainflag1,]
  data2=data[!trainflag1,]
  pvalue1=sapply(14:ncol(data1),function(j){
      if (opt=="del")
      {
        M=table(data1[,"platinumclass"],data1[,j]<0)
      }else
      {
        M=table(data1[,"platinumclass"],data1[,j]>0)
      }
     
      if (ncol(M)>1)
      {
        res=fisher.test(M)$p.value
      }else
      {
        res=NA
      }
    })
  res$pvalue1=order(pvalue1)  
  pvalue2=sapply(14:ncol(data1),function(j){
    if (opt=="del")
    {
      M=table(data2[,"platinumclass"],data2[,j]<0)
    }else
    {
      M=table(data2[,"platinumclass"],data2[,j]>0)
    }
    if (ncol(M)>1)
    {
      res=fisher.test(M)$p.value
    }else
    {
      res=NA
    }
  })
  res$pvalue2=order(pvalue2)
  return(res)
}

obres=data.frame(matrix(NA,nrow=ncol(alltraindata)-13,ncol=0))
for (it in 1:10)
{
  res=compute2pvalues(data=alltraindata,it=it)
  obres=cbind(obres,res)
}

nullres=data.frame(matrix(NA,nrow=ncol(alltraindata)-13,ncol=0))
for (num in 1:10)
{
  cat(num,"...")
  alltraindata1=alltraindata
  set.seed(num)
  alltraindata1[,"platinumclass"]=alltraindata1[,"platinumclass"][sample(1:nrow(alltraindata1))]
  for (it in 1:100)
  {
    res=compute2pvalues(data=alltraindata1,it=it)
    nullres=cbind(nullres,res)
  }
}

data_class_copynumber_rand=list(data1=alldata[trainflag,],data2=alldata[!trainflag,])
alldata1=rbind(data_class_copynumber_hrgenes$data1,data_class_copynumber_hrgenes$data2)
data_class_copynumber_hrgenes_rand=list(data1=alldata1[trainflag,],data2=alldata1[!trainflag,])
class_cat_copynumber=classify_platinum(data=data_class_copynumber_rand,platform="copynumber_cat",selgenes=NULL,includeclinical=T)
class_cat_copynumber1=classify_platinum(data=data_class_copynumber_rand,platform="copynumber_cat",selgenes=NULL,includeclinical=F)
class_cat_copynumber_hrgenes=classify_platinum(data=data_class_copynumber_hrgenes,platform="copynumber_cat_hrgenes",selgenes=NULL,includeclinical=T)


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

data1=data4_copynumber$data1
data2=data4_copynumber$data2
data1=data4_mrna$data1
data2=data4_mrna$data2
data1=data_copynumber$data1
data2=data_copynumber$data2
drawmyroc=function(data=data4_copynumber,includeclinical=F,opt="rmborderline")
{
  data1=data$data1
  data2=data$data2
  if ("borderlinesensitive" %in% data1[,"platinumclass"])
  {
    if (opt=="rmborderline")
    {
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive",]
      data1[data1[,"platinumclass"]=="refractory","platinumclass"]="resistant"
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive",]
      data2[data2[,"platinumclass"]=="refractory","platinumclass"]="resistant"
      data1[,"platinumclass"]=as.factor(as.character(data1[,"platinumclass"]))
      data2[,"platinumclass"]=as.factor(as.character(data2[,"platinumclass"]))
    }else
    {
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive" & ! data1[,"platinumclass"] %in% "resistant",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive" & ! data2[,"platinumclass"] %in% "resistant",]
      data1[,"platinumclass"]=as.factor(as.character(data1[,"platinumclass"]))
      data2[,"platinumclass"]=as.factor(as.character(data2[,"platinumclass"]))
    }
    
  }
  data_copynumber_extreme=list(data1=data1,data2=data2)
  class_copynumber_extreme=classify_platinum(data=data_copynumber_extreme,platform=" ",selgenes=NULL,includeclinical=includeclinical)
  
  data1=data_copynumber_commc$data1
  data2=data_mrna_commc$data1
  selgenes=class_copynumber_extreme$coeff
  selgenes=selgenes[!names(selgenes) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease","age")]
  selgenes=names(selgenes)
  
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
  class_copynumber_extreme1=classify_platinum(data=data_copynumber_extreme,platform="copynumber_extreme",selgenes=selgenes1,includeclinical=includeclinical)
  test=drawroc_onselectedgenes(data=data_copynumber_extreme,selgenes = selgenes1)
  
  result=list(selgenes=selgenes,selgenes1=selgenes1)
}
res_copynumber_rmborderline_inc=drawmyroc(data=data4_copynumber,includeclinical = T)
res_copynumber_rmborderline=drawmyroc(data=data4_copynumber,includeclinical = F)
res_copynumber_extreme_inc=drawmyroc(data=data4_copynumber,includeclinical = T,opt="rm2")
res_copynumber_extreme=drawmyroc(data=data4_copynumber,includeclinical = F,opt="rm2")
res_copynumber_inc=drawmyroc(data=data_copynumber,includeclinical = T)
res_copynumber=drawmyroc(data=data_copynumber,includeclinical = F)
