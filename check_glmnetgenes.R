#!/usr/bin/env Rscript

coeff=class_copynumber$coeff
coeff=coeff[! names(coeff) %in% c("intercept","residual_disease_largest_noduleNo_Macroscopic_disease")]
pvalues_train=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_train.txt",header=T)
pvalues_train_order=pvalues_train[order(pvalues_train[,1]),1]
names(pvalues_train_order)=rownames(pvalues_train)[order(pvalues_train[,1])]
pvalues=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
pvalues_order=pvalues_train[order(pvalues[,1]),1]
names(pvalues_order)=rownames(pvalues)[order(pvalues[,1])]

analysegenes=function(data=data_copynumber,coeff=coeff)
{
  data=dataselgenes(data,names(coeff))
  genes=smallfdr_copynumber$gene
  genes=names(pvalues_train_order)[1:50]
  data=dataselgenes(data,genes)
  #update stage and residual, combine training and testing into one dataframe
  alldata=formwholedataframe(data)
  
  #add indicator of train/test as the last column
  train=rep(FALSE,nrow(alldata))
  train[1:nrow(data$data1)]=TRUE
  alldata=cbind(alldata,train=train)
  
  #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
  alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
  alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  
  #clinicals to remove before analysis, they are not useful for the model. only age and residual_disease are kept
  remove_clinicals=c("tumor_grade","clinical_stage","drug_interval_computed","race","initial_pathologic_dx_year","vital_status","death_months_to","treatment_outcome_first_course","progression_free_survival","uselastcontact")
  alldata=alldata[,!colnames(alldata) %in% remove_clinicals]
  
  #number of clinical variables to keep (age,residual)
  num_clinical=2
  
  #some genes have name with "-" character,replace it to "_". Or it would cause problem in modeling
  colnames(alldata)=gsub("-","_",colnames(alldata))
  #the column indices of train indicator and outcome
  idx_train=which(colnames(alldata)=="train")
  #the column of output
  idx_y=which(colnames(alldata)=="platinumclass")
  
  #form training data
  ally=alldata[,idx_y]
  x=alldata[which(alldata[,idx_train]==1),]
  y=ally[which(alldata[,idx_train]==1)]
  x=x[,-c(idx_y,idx_train)]
  x=x[,(num_clinical+1):ncol(x)]
  x=scale(x)
  #form testing data
  x1=alldata[which(alldata[,idx_train]==0),]
  y1=ally[which(alldata[,idx_train]==0)]
  x1=x1[,-c(idx_y,idx_train)]
  x1=x1[,(num_clinical+1):ncol(x1)]
  x1=scale(x1)
  
  
  trainingdata=data.frame(y,x)
  glmfit_train=glm(as.formula(y~.),data=trainingdata,family="binomial")
  tmp=coef(summary(glmfit_train))

  idx=NULL
  for (i in 2:ncol(x))
  {
    idx=c(idx,which(colnames(x)==rownames(tmp)[i]))
  }
  b=sum(c(1,x[1,idx])*tmp[,1])
  c=exp(b)/(1+exp(b))
  pfit_train=predict(glmfit_train,newdata=trainingdata,type="response")
  testingdata=data.frame(y=y1,x1)
  glmfit_test=glm(as.formula(y~.),data=testingdata,family="binomial")
  pfit_test=predict(glmfit_test,newdata=testingdata,type="response")
  plotroc2(pfit_train,y,pfit_test,y1)
  
  x_=x
  x1_=x1
  for (i in 1:length(coeff))
  {
    x_[,i]=coeff[i]*x[,i]
    x1_[,i]=coeff[i]*x1[,i]
  }
  trainingdata_=data.frame(y,x_)
  testingdata_=data.frame(y=y1,x1_)
  
  #visualization
  pc=princomp(x)
  ycolor=rep("green",length(y))
  ycolor[which(as.character(y)=="Sensitive")]="red"
  plot(pc$scores[,1],pc$scores[,2],col=ycolor)
  pc1=princomp(x1)
  ycolor1=rep("green",length(y1))
  ycolor1[which(as.character(y1)=="Sensitive")]="red"
  plot(pc1$scores[,1],pc1$scores[,2],col=ycolor1)
  
  pc=princomp(x_)
  ycolor=rep("green",length(y))
  ycolor[which(as.character(y)=="Sensitive")]="red"
  plot(pc$scores[,3],pc$scores[,2],col=ycolor)
  pc1=princomp(x1_)
  ycolor1=rep("green",length(y1))
  ycolor1[which(as.character(y1)=="Sensitive")]="red"
  plot(pc1$scores[,3],pc1$scores[,2],col=ycolor1)
  
  #KNN
  library(class)
  set.seed(1)
  knn.pred=knn(x,x1,y,k=10)
  table(knn.pred,y1)
  mean(knn.pred==y1)
  set.seed(1)
  knn.pred1=knn(x_,x1_,y,k=10)
  table(knn.pred1,y1)
  mean(knn.pred1==y1)
  
  #SVM
  library(e1071)
  colnames(trainingdata)=gsub(".","_",colnames(trainingdata),fixed=T)
  colnames(testingdata)=gsub(".","_",colnames(testingdata),fixed=T)
  colnames(x)=gsub(".","_",colnames(x),fixed=T)
  colnames(x1)=gsub(".","_",colnames(x1),fixed=T)
  colnames(trainingdata)=gsub("|","_",colnames(trainingdata),fixed=T)
  colnames(testingdata)=gsub("|","_",colnames(testingdata),fixed=T)
  colnames(x)=gsub("|","_",colnames(x),fixed=T)
  colnames(x1)=gsub("|","_",colnames(x1),fixed=T)
  svm.model <- svm(y ~ ., data = trainingdata, cost = 100, gamma = 1)
  svm.model <- svm(y ~ ., data = trainingdata,degree=3,kernel="sigmoid")
  svm.pred <- predict(svm.model, x1)
  table(svm.pred,y1)
  mean(svm.pred==y1)
  svm.model1 <- svm(y ~ ., data = testingdata_)
  svm.pred1 <- predict(svm.model1, x1_)
  table(svm.pred1,y1)
  
  #decision tree
  library(party)
  tree=ctree(y~.,data=trainingdata)
  plot(tree)
  tree.pred=predict(tree,newdata=as.data.frame(x1))
  table(tree.pred,y1)
  
  library(rpart)
  tree=rpart(y~.,data=trainingdata)
  plot(tree)
  tree.pred=predict(tree,newdata=as.data.frame(x))
  tree.pred1=predict(tree,newdata=as.data.frame(x1))
  plotroc2(tree.pred[,1],y,tree.pred1[,1],y1)
  
  library(randomForest)
  rf=randomForest(y~.,data=trainingdata,ntree=100,proximity=T)
  rf.pred=predict(rf,newdata=x1)
  table(rf.pred,y1)
  varImpPlot(rf)
  plot(margin(rf,y))
}

