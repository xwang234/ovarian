
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_glmnetclassresult.RData")
xdat=x[,8:ncol(x)]
x1dat=x1[,8:ncol(x1)]
res=sapply(2:ncol(xdat),function(i){
  res1=cor(xdat[,i],xdat[,i-1],use="complete")
})
breakp=which(res<0.95)
breakp=c(0,breakp,ncol(xdat))

xdata=data.frame(matrix(NA,nrow=nrow(xdat),ncol=length(breakp)-1))
rownames(xdata)=rownames(x)
colnames(xdata)=colnames(xdat)[breakp[1:(length(breakp)-1)]+1]
x1data=data.frame(matrix(NA,nrow=nrow(x1dat),ncol=length(breakp)-1))
colnames(x1data)=colnames(x1dat)[breakp[1:(length(breakp)-1)]+1]
rownames(x1data)=rownames(x1)
#form segment data
rowMedians=function(xdat)
{
  res=apply(xdat,1,median)
  return(res)
}
for (i in 1:(length(breakp)-1))
{
  start=breakp[i]+1
  end=breakp[i+1]
  if (end>start)
  {
#     xdata[,i]=rowMeans(xdat[,start:end])
#     x1data[,i]=rowMeans(x1dat[,start:end])
    xdata[,i]=rowMedians(xdat[,start:end])
    x1data[,i]=rowMedians(x1dat[,start:end])
  }else
  {
    xdata[,i]=xdat[,start]
    x1data[,i]=x1dat[,start]
  }
}
xdata=cbind(residual=x[,7],xdata)
x1data=cbind(residual=x1[,7],x1data)
penalty=rep(1,ncol(xdata))
penalty[c(1)]=0
#trainresult=selectlambda(x,y,numlambda=200,ntime=20,opt="binomial",penalty=rep(1,ncol(x)),standardize=T,main=prefix,nfolds=10)
trainresult=selectlambda(xdata,y,numlambda=100,ntime=1000,opt="binomial",penalty=penalty,standardize=T,main=prefix,nfolds=10)
lambda.sel=trainresult$lambda.min
print(lambda.sel)
#[1] [1] 0.0206985

#fit=glmnet(as.matrix(x),y,family="binomial",standardize=T,nlambda = 200)
fit=glmnet(as.matrix(xdata),y,family="binomial",nlambda = 100,penalty.factor=penalty,standardize=T)
tmp=predict(fit,type="coefficients",s=lambda.sel)
coeff=tmp@x
names(coeff)=c("intercept",colnames(xdata)[tmp@i[tmp@i>0]])
pfit_train = predict(fit,as.matrix(xdata),s=lambda.sel,type="response")
pfit_test = predict(fit,as.matrix(x1data),s=lambda.sel,type="response")
plotroc2(pfit_train,y,pfit_test,y1,main="ROC")

