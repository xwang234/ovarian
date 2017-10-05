#!/usr/bin/env Rscript

getlikelihoodmatrix=function(traindata=traindata)
{
  res=data.frame(matrix(NA,nrow=0,ncol=6))
  tmp=data.frame(matrix(NA,nrow=1,ncol=6))
  colnames(tmp)=colnames(res)=c("gene","Sensitive","Resistant","Type","pvalue","proportion")
  idx_y=which(colnames(traindata)=="platinumclass")
  idx_hrgenes=which(colnames(traindata)=="hrgenes")
  if ("hrgenes" %in% colnames(traindata))
  {
    tmp$gene="hrgenes"
    tmp$Sensitive=sum(traindata[traindata$platinumclass=="Sensitive",idx_hrgenes]>0)/sum(traindata$platinumclass=="Sensitive")
    tmp$Resistant=sum(traindata[traindata$platinumclass=="Resistant",idx_hrgenes]>0)/sum(traindata$platinumclass=="Resistant")
    tmp$Type="other"
    M=table(traindata[,"hrgenes"]>0,traindata[,"platinumclass"])
    tmp$pvalue=fisher.test(M)$p.value
    tmp$proportion=sum(traindata$hrgenes>0)/nrow(traindata)
    res=rbind(res,tmp)
    idx_y=idx_y+1 #the index of del/amps
  }
  
  for (i in (idx_y+1):ncol(traindata))
  {
    if (sum(traindata[,i]<0,na.rm=T)>0)
    {
      gene=colnames(traindata)[i]
      l_sensitive=sum(traindata[traindata$platinumclass=="Sensitive",i]<0)/sum(traindata$platinumclass=="Sensitive")
      l_resistant=sum(traindata[traindata$platinumclass=="Resistant",i]<0)/sum(traindata$platinumclass=="Resistant")
      M=table(traindata$platinumclass,traindata[,i]<0)
      pvalue=fisher.test(M)$p.value
      prop=sum(traindata[,i]<0)/nrow(traindata)
      if (pvalue<0.1)
      {
        tmp[1,]=c(gene,l_sensitive,l_resistant,"del",pvalue,prop)
        res=rbind(res,tmp)
      }
    }
    if (sum(traindata[,i]>0,na.rm=T)>0)
    {
      gene=colnames(traindata)[i]
      l_sensitive=sum(traindata[traindata$platinumclass=="Sensitive",i]>0)/sum(traindata$platinumclass=="Sensitive")
      l_resistant=sum(traindata[traindata$platinumclass=="Resistant",i]>0)/sum(traindata$platinumclass=="Resistant")
      M=table(traindata$platinumclass,traindata[,i]>0)
      pvalue=fisher.test(M)$p.value
      prop=sum(traindata[,i]>0)/nrow(traindata)
      if (pvalue<0.1)
      {
        tmp[1,]=c(gene,l_sensitive,l_resistant,"amp",pvalue,prop)
        res=rbind(res,tmp)
      }
    }
  }
  res$Sensitive=as.numeric(res$Sensitive)
  res$Resistant=as.numeric(res$Resistant)
  res$pvalue=as.numeric(res$pvalue)
  res$proportion=as.numeric(res$proportion)
  return(res)
}
bayesclass=function(newdata=testdata,usedlikelihoodmatrix=likelihoodmatrix,prior_sensitive,prior_resistant)
{
  idx_y=which(colnames(newdata)=="platinumclass")
  x=newdata[,(idx_y+1):ncol(newdata)]
  y=newdata$platinumclass
  
  res=data.frame(matrix(NA,ncol=3,nrow=nrow(newdata)))
  colnames(res)=c("Sensitive","Resistant","Class")
  for (i in 1:nrow(newdata))
  {
    prob_sensitive=prior_sensitive
    prob_resistant=prior_resistant
    #check hrgenes
    if ("hrgenes" %in% usedlikelihoodmatrix$gene)
    {
      idx_hrgenes=which(usedlikelihoodmatrix$gene=="hrgenes")
      if (newdata[i,"hrgenes"]>0)
      {
        prob_sensitive=prob_sensitive*usedlikelihoodmatrix$Sensitive[idx_hrgenes]
        prob_resistant=prob_resistant*usedlikelihoodmatrix$Resistant[idx_hrgenes]
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
      idx2=which(usedlikelihoodmatrix$gene %in% names(idx_del) & usedlikelihoodmatrix$Type=="del")
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
      idx2=which(usedlikelihoodmatrix$gene %in% names(idx_amp) & usedlikelihoodmatrix$Type=="amp")
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

bayesclassifier=function(data=data_class_copynumber,usedselgenes=selgenes,main)
{
  colnames(data$data1)=colnames(data$data2)=gsub("|","___",colnames(data$data1),fixed=T)
  usedselgenes=gsub("|","___",usedselgenes,fixed=T)
  usedselgenes=c("hrgenes",usedselgenes)
  subdata=drawroc_onselectedgenes(selgenes = usedselgenes,data=data)
  traindata=subdata$traindata
  testdata=subdata$testdata
  
  traindata=traindata[,c(11,14:ncol(traindata))]
  likelihoodmatrix=getlikelihoodmatrix(traindata=traindata)
  prior_sensitive=sum(traindata$platinumclass=="Sensitive")/nrow(traindata)
  prior_resistant=1-prior_sensitive
  testdata=testdata[,c(11,14:ncol(testdata))]
  trainres=bayesclass(newdata=traindata,usedlikelihoodmatrix=likelihoodmatrix,prior_sensitive=prior_sensitive,prior_resistant=prior_resistant)$mat
  trainres=cbind(trainres,traindata$platinumclass)
  table(trainres[,3],trainres[,4])
  testres=bayesclass(newdata=testdata,usedlikelihoodmatrix=likelihoodmatrix,prior_sensitive=prior_sensitive,prior_resistant=prior_resistant)$mat
  testres=cbind(testres,testdata$platinumclass)
  table(testres[,3],testres[,4])
  pfit_train=bayesclass(newdata=traindata,usedlikelihoodmatrix = likelihoodmatrix,prior_sensitive=prior_sensitive,prior_resistant=prior_resistant)$mat[,2]
  pfit_test=bayesclass(newdata=testdata,usedlikelihoodmatrix = likelihoodmatrix,prior_sensitive=prior_sensitive,prior_resistant=prior_resistant)$mat[,2]
  plotroc2(pfit_train,traindata$platinumclass,pfit_test,testdata$platinumclass,plotflag=1,main=main)
}
#check the categorical copynumber data:
load("/fh/fast/dai_j/CancerGenomics/data/platinum_classification_copynumber_category.RData")

class_cat_copynumber=classify_platinum(data=data_class_copynumber,platform="copynumber_cat",selgenes=NULL,includeclinical=F)
selgenes=class_cat_copynumber$coeff
selgenes=selgenes[!names(selgenes) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease","age")]
selgenes=names(selgenes)

alldata=rbind(data_class_copynumber$data1,data_class_copynumber$data2)
#only use train
alldata=data_copynumber$data1
res=data.frame(matrix(NA,nrow=500,ncol=100))
for (i in 1:100)
{
cat(i,"..")
trainflag=rep(T,nrow(alldata))
set.seed(i)
testidx=sample(1:nrow(alldata),ceiling(nrow(alldata)/3))
trainflag[testidx]=F
data_class_copynumber_rand=list(data1=alldata[trainflag,],data2=alldata[!trainflag,])
table(data_class_copynumber_rand$data1[,"platinumclass"])
table(data_class_copynumber_rand$data2[,"platinumclass"])

class_cat_copynumber_rand=classify_platinum(data=data_class_copynumber_rand,platform="copynumber_cat_rand",selgenes=NULL,includeclinical=F)
selgenes=class_cat_copynumber_rand$coeff
selgenes=selgenes[!names(selgenes) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease","age")]
selgenes=names(selgenes)
if (length(selgenes)>0)
{
  res[1:length(selgenes),i]=selgenes
}



# bayesclassifier(data=data_class_copynumber_rand,usedselgenes=selgenes,main=paste0("iteration:",i))
}
res1=unlist(res)
res1=res1[!is.na(res1)]
uniq_genes=unique(res1)
res2=data.frame(matrix(NA,nrow=length(uniq_genes),ncol=2))
colnames(res2)=c("gene","count")
res2$gene=uniq_genes
for (i in 1:length(uniq_genes))
{
  res2$count[i]=sum(res1==uniq_genes[i])
}
res2=res2[order(res2$count,decreasing=T),]
