#!/usr/bin/env Rscript

ngenes=20000
nsamples=380
sdata=data.frame(matrix(NA,nrow=nsamples,ncol=ngenes))
for (i in 1:ngenes)
{
  set.seed(i)
  sdata[,i]=rnorm(nsamples)
}


ngenesused=10
ngenesused=50
sy=rep(0,nsamples)
for (i in 1:nsamples)
{
  if (sum(sdata[i,1:ngenesused])>0)
  {
    sy[i]=1
  }
}

x=sdata[1:280,]
x1=sdata[281:nrow(sdata),]
y=sy[1:280]
y1=sy[281:nrow(sdata)]

numsel=10000
trainresult=selectlambda(x[,1:numsel],y,numlambda=100,ntime=20,opt="binomial",penalty=penalty,main=platform,nfolds=10)
#lambda.sel=trainresult$lambda.1se
lambda.sel=trainresult$lambda.min

fit=glmnet(as.matrix(x[,1:numsel]),y,family="binomial",nlambda = 100,penalty.factor=penalty)
tmp=predict(fit,type="coefficients",s=lambda.sel)
coeff=tmp@x
names(coeff)=c("intercept",colnames(x[,1:numsel])[tmp@i[tmp@i>0]])
selectedgenes=names(coeff)[! names(coeff) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease")]
sum(selectedgenes %in% colnames(x)[1:ngenesused])
pfit_train = predict(fit,as.matrix(x[,1:numsel]),s=lambda.sel,type="response")
pfit_test = predict(fit,as.matrix(x1[,1:numsel]),s=lambda.sel,type="response")

plotroc2(pfit_train,y,pfit_test,y1,main=paste0("ROC, ",length(selectedgenes), " genes selected"))
plotgenecor(selectedgenes,ngenesused)
plotgenecor=function(selectedgenes,ngenesused)
{
  x2=rbind(x,x1)
  idx=which(colnames(x2) %in% selectedgenes)
#   corr=matrix(NA,nrow=length(idx),ncol=length(idx))
#   colnames(corr)=rownames(corr)=colnames(x2)[idx]
#   for (i in 1:length(idx)) corr[i,i]=1
#   for (i in 1:(length(idx)-1))
#   {
#     for (j in (i+1):length(idx))
#     {
#       corr[i,j]=corr[j,i]=cor(x2[,idx[i]],x2[,idx[j]])
#     }
#   }
  idx1=1:ngenesused # the used genes
  corr=matrix(NA,nrow=length(idx),ncol=length(idx1))
  colnames(corr)=colnames(x2)[idx1]
  rownames(corr)=colnames(x2)[idx]
  for (i in 1:(length(idx)))
  {
    for (j in 1:length(idx1))
    {
      corr[i,j]=cor(x2[,idx[i]],x2[,idx1[j]])
    }
  }
  heatmap.2(corr,dendrogram = "none")
  
}
sapply_compute_pvalue=function(x=x,x1=x1,y=y,y1=y1)
{
  x=rbind(x,x1)
  y=c(y,y1)
  y=factor(y)
  res=sapply(1:ncol(x),function(i){
    tmp=glm(y~x[,i],family="binomial")
    summaryres=coef(summary(tmp))
    res1=NA
    if (nrow(summaryres)>1)
    {
      res1=summaryres[2,4]
    }
    return(res1)
  })
  return(res)
}

spvalues=sapply_compute_pvalue(x=x,x1=x1,y=y,y1=y1)
hist(spvalues,probability = T)
tmp=quantile(spvalues,c(0.01,0.05,0.1))
sum(spvalues[1:ngenesused]<=tmp[1])
sum(spvalues[1:ngenesused]<=tmp[2])
sum(spvalues[1:ngenesused]<=tmp[3])
tmp=colMeans(sdata)
sum(mean(sdata[,41])>tmp)
sfdrs=sapply(1:ngenesused,function(i){
  res=sum(spvalues[i]>=spvalues[(ngenesused+1):length(spvalues)])/(length(spvalues)-ngenesused)/(sum(spvalues[i]>=spvalues)/ngenesused)
})
