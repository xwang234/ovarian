#!/usr/bin/env Rscript

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")

alldata=rbind(data_copynumber$data1,data_copynumber$data2)
which(colnames(alldata)=="UNC5C")
#[1] 6016
hist(alldata[alldata[,"platinumclass"]=="Sensitive",6016])
hist(alldata[alldata[,"platinumclass"]=="Resistant",6016])

class_alldata=alldata
class_alldata[,14:ncol(alldata)]=0
for (i in 14:ncol(alldata))
{
  idx=alldata[,i]>=0.5
  class_alldata[idx,i]=1
  idx=alldata[,i]<= -0.5
  class_alldata[idx,i]=-1
}


data_class_copynumber=list(data1=class_alldata[1:nrow(data_copynumber$data1),],data2=class_alldata[(nrow(data_copynumber$data1)+1):nrow(class_alldata),])
save(data_class_copynumber,file="/fh/fast/dai_j/CancerGenomics/data/platinum_classification_copynumber_category.RData")
pvalues_class_alldata=compute_pvalues(data=data_class_copynumber,platform="class_copynumber",runpermutation=F)
hist(class_alldata[class_alldata[,"platinumclass"]=="Sensitive",6016])
hist(class_alldata[class_alldata[,"platinumclass"]=="Resistant",6016])
t.test(class_alldata[class_alldata[,"platinumclass"]=="Sensitive",6016],class_alldata[class_alldata[,"platinumclass"]=="Resistant",6016])


newdata=class_alldata[,14:ncol(class_alldata)]
selgenes=names(pvalues_class_alldata)[which(pvalues_class_alldata<=0.1)]
newdata=newdata[,which(colnames(newdata) %in% selgenes)]

idxsamples=rep(F,nrow(newdata))
set.seed(2)
idxsamples[sample(1:nrow(newdata),ceiling(nrow(newdata)/2))]=T
traindata=newdata[idxsamples,]
testdata=newdata[!idxsamples,]
trainy=class_alldata[idxsamples,"platinumclass"]
testy=class_alldata[!idxsamples,"platinumclass"]
test=glm(trainy~traindata[,413],family = "binomial")


library(caret)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
control <- rfeControl(functions=rfFuncs, method="boot", number=100)

control <- rfeControl(functions=treebagFuncs, method="cv", number=10)


rfemod <- rfe(traindata, trainy, sizes=c(1:50), rfeControl=control)
rfemod <- rfeIter(traindata, trainy,testX=traindata[1:100,],testY=testy[1:100],sizes=c(1:50), rfeControl=control)
predictors(rfemod)
plot(rfemod)

pfit_train=predict(rfemod,traindata)$Sensitive
pfit_test=predict(rfemod,testdata)$Sensitive
plotroc2(pfit_train,trainy,pfit_test,testy,main=paste0("ROC, rfe"))

library(randomForest)
library(mlbench)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
seed <- 7
metric <- "Accuracy"
mtry <- sqrt(ncol(traindata))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(traindata,trainy, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)

ctrl <- rfeControl(functions = nbFuncs,method = "repeatedcv",repeats = 5,verbose = FALSE)
rfemod <- rfe(traindata, trainy, sizes=c(1:50), rfeControl=ctrl)
