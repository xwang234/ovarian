#!/usr/bin/env Rscript
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")

hr30genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/HR30genes.txt",sep="\t",header=T,stringsAsFactors=F)
hr30genes[c(17,23),]
hr30genes=rbind(hr30genes,hr30genes[c(17,23),])
hr30genes[17,1]="MSH2"
hr30genes[31,1]="EPCAM"
hr30genes[23,1]="PTEN"
hr30genes[32,1]="KILLIN"
hr30genes=hr30genes$Gene

hrgenes=c("ATM","BARD1","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FAM175A","NBN","PALB2","RAD51C","RAD51D","MRE11A")

extractselgenesfromdatausealias=function(data=data_rawmutation,selgenes=hrgenes,num_allclinical=13)
{
  mapres=mapaliasgenes(selgenes,colnames(data$data1)[(num_allclinical+1):ncol(data$data1)])
  idxselgenes=which(!is.na(mapres$idx2))
  idxgenesindata=mapres$idx2[idxselgenes]

  res=list(data1=cbind(data$data1[,1:num_allclinical],data$data1[,idxgenesindata+num_allclinical]),
             data2=cbind(data$data2[,1:num_allclinical],data$data2[,idxgenesindata+num_allclinical]),
           idxselgenes=idxselgenes,idxgenesindata=idxgenesindata)
}

Sys.time()
data_hrgenes_mutation=extractselgenesfromdatausealias(data=data_rawmutation,selgenes=hrgenes,num_allclinical=13)
Sys.time()
# #remove constant columns
# tmp=sapply(14:ncol(data_hrgenes_mutation$data1),function(i){
#   res=T
#   if (var(data_hrgenes_mutation$data1[,i])==0)
#   {
#     res=F
#   }
#   return(res)
# } )
# tmp=c(rep(T,13),tmp)
# data_hrgenes_mutation$data1=data_hrgenes_mutation$data1[,tmp]
# data_hrgenes_mutation$data2=data_hrgenes_mutation$data2[,tmp]
# fm=as.formula(paste0("platinumclass~",paste0(colnames(data_hrgenes_mutation$data1)[14:ncol(data_hrgenes_mutation$data1)],collapse="+")))
# glm_hrgenes_mutation_train=glm(fm,data=data_hrgenes_mutation$data1,family="binomial")
# pfit_train=predict(glm_hrgenes_mutation_train,type="response")
# pfit_test=predict(glm_hrgenes_mutation_train,newdata=data_hrgenes_mutation$data2,type="response")
# plotroc2(pfit_train,data_hrgenes_mutation$data1$platinumclass,pfit_test,data_hrgenes_mutation$data2$platinumclass,main=paste0("ROC"))
# ldafit=lda(fm,data=data_hrgenes_mutation$data1)
# pfit_train=predict(ldafit,type="response")$posterior[,1]
# pfit_test=predict(ldafit,newdata=data_hrgenes_mutation$data2,type="response")$posterior[,1]
# plotroc2(pfit_train,data_hrgenes_mutation$data1$platinumclass,pfit_test,data_hrgenes_mutation$data2$platinumclass,main=paste0("ROC,lda"))

#first apply classfication, then form model based on seleged genes, and draw ROC
drawclasscalledrocs=function(data=data_hrgenes_mutation,includeclinical=T)
{
  class_hrgenesmutation=classify_platinum(data,platform="mutation",includeclinical=includeclinical)
  
  if (length(names(class_hrgenesmutation$coeff[-1]))>0)
  {
    fm=as.formula(paste0("platinumclass~",paste0(colnames(class_hrgenesmutation$x),collapse="+")))
    datatrain=cbind(platinumclass=class_hrgenesmutation$y,class_hrgenesmutation$x)
    glm_hrgenes_mutation_train=glm(fm,data=datatrain,family="binomial")
    print(summary(glm_hrgenes_mutation_train))
    pfit_train=predict(glm_hrgenes_mutation_train,type="response")
    datatest=cbind(platinumclass=class_hrgenesmutation$y1,class_hrgenesmutation$x1)
    pfit_test=predict(glm_hrgenes_mutation_train,newdata=datatest,type="response")
    plotroc2(pfit_train,class_hrgenesmutation$y,pfit_test,class_hrgenesmutation$y1,main=paste0("ROC,glm all variables"))
    glm_hrgenes_mutation_test=glm(fm,data=datatest,family="binomial")
    print(summary(glm_hrgenes_mutation_test))
    fm1=as.formula(paste0("platinumclass~",paste0(names(class_hrgenesmutation$coeff[-1]),collapse="+")))
    glm_hrgenes_mutation_train1=glm(fm1,data=datatrain,family="binomial")
    print(summary(glm_hrgenes_mutation_train1))
    pfit_train1=predict(glm_hrgenes_mutation_train1,type="response")
    pfit_test1=predict(glm_hrgenes_mutation_train1,newdata=datatest,type="response")
    plotroc2(pfit_train1,class_hrgenesmutation$y,pfit_test1,class_hrgenesmutation$y1,main=paste0("ROC,glm selected variables"))
    glm_hrgenes_mutation_test1=glm(fm1,data=datatest,family="binomial")
    print(summary(glm_hrgenes_mutation_test1))
  }
  if (includeclinical)
  {
    glm_hrgenes_mutation_clinical_train=glm(as.formula("platinumclass~age+residual_disease_largest_noduleNo_Macroscopic_disease"),data=datatrain,family="binomial")
    pfit_train=predict(glm_hrgenes_mutation_clinical_train,type="response")
    pfit_test=predict(glm_hrgenes_mutation_clinical_train,newdata=datatest,type="response")
    plotroc2(pfit_train,class_hrgenesmutation$y,pfit_test,class_hrgenesmutation$y1,main=paste0("ROC,glm clinical"))
    
    glm_hrgenes_mutation_residual_train=glm(as.formula("platinumclass~residual_disease_largest_noduleNo_Macroscopic_disease"),data=datatrain,family="binomial")
    pfit_train=predict(glm_hrgenes_mutation_residual_train,type="response")
    pfit_test=predict(glm_hrgenes_mutation_residual_train,newdata=datatest,type="response")
    plotroc2(pfit_train,class_hrgenesmutation$y,pfit_test,class_hrgenesmutation$y1,main=paste0("ROC,glm residual"))
    
    glm_hrgenes_mutation_age_train=glm(as.formula("platinumclass~age"),data=datatrain,family="binomial")
    pfit_train=predict(glm_hrgenes_mutation_age_train,type="response")
    pfit_test=predict(glm_hrgenes_mutation_age_train,newdata=datatest,type="response")
    plotroc2(pfit_train,class_hrgenesmutation$y,pfit_test,class_hrgenesmutation$y1,main=paste0("ROC,glm age"))
  }
}

drawrocs(data=data_hrgenes_mutation,includeclinical=T)
drawrocs(data=data_hrgenes_mutation,includeclinical=F)
#randomize samples into training set
alldata=rbind(data_hrgenes_mutation$data1,data_hrgenes_mutation$data2)
set.seed(2)
idx=sample(nrow(alldata))
alldata=alldata[idx,]
numtrain=as.integer(nrow(alldata)/2)+1
data_hrgenes_mutation1=list(data1=alldata[1:numtrain,],data2=alldata[(numtrain+1):nrow(alldata),])
drawrocs(data=data_hrgenes_mutation1,includeclinical=T)
drawrocs(data=data_hrgenes_mutation1,includeclinical=F)


data_hrgenes_copynumber=extractselgenesfromdatausealias(data=data_copynumber,selgenes=hrgenes,num_allclinical=13)
class_hrgenescopynumber=classify_platinum(data=data_hrgenes_copynumber,platform="copynumber",includeclinical=T)

data_hrgenes_mrna=extractselgenesfromdatausealias(data=data_mrna,selgenes=hrgenes,num_allclinical=13)
class_hrgenesmrna=classify_platinum(data=data_hrgenes_mrna,platform="mrna",includeclinical=T)

#t-test of data from 30 genes
alldata=data_rawmutation$data1
alldata=rbind(data_rawmutation$data1,data_rawmutation$data2)
hrgenes_sensitive=NULL
hrgenes_resistant=NULL

for (i in 1:nrow(alldata))
{
  if (alldata[i,"platinumclass"]=="Sensitive")
  {
    tmp=alldata[i,colnames(alldata) %in% hrgenes]
    tmp=sum(tmp>0,na.rm=T)
    names(tmp)=rownames(alldata)[i]
    hrgenes_sensitive=c(hrgenes_sensitive,tmp)
  }
  if (alldata[i,"platinumclass"]=="Resistant")
  {
    tmp=alldata[i,colnames(alldata) %in% hrgenes]
    tmp=sum(tmp>0,na.rm=T)
    names(tmp)=rownames(alldata)[i]
    hrgenes_resistant=c(hrgenes_resistant,tmp)
  }
}
length(hrgenes_sensitive)
length(hrgenes_resistant)
hist(hrgenes_sensitive)
hist(hrgenes_resistant)
table(hrgenes_sensitive)
table(hrgenes_resistant)
t.test(hrgenes_sensitive,hrgenes_resistant)

test=glm(c(rep(1,length(hrgenes_sensitive)),rep(0,length(hrgenes_resistant)))~c(hrgenes_sensitive,hrgenes_resistant),family="binomial")
summary(test)
#dataframe to keep the frequency
idxgenes=which(colnames(alldata) %in% hrgenes)
df_hrgenes=data.frame(matrix(NA,ncol=1+length(idxgenes),nrow=length(hrgenes_sensitive)+length(hrgenes_resistant)))
colnames(df_hrgenes)=c("type",colnames(alldata)[idxgenes])
rownames(df_hrgenes)=c(names(hrgenes_sensitive),names(hrgenes_resistant))
df_hrgenes$type=c(rep("sensitive",length(hrgenes_sensitive)),rep("resistant",length(hrgenes_resistant)))
for (i in 1:nrow(df_hrgenes))
{
  for (j in 2:ncol(df_hrgenes))
  {
    tmp1=which(rownames(alldata)==rownames(df_hrgenes)[i])
    tmp2=which(colnames(alldata)==colnames(df_hrgenes)[j])
    df_hrgenes[i,j]=alldata[tmp1,tmp2]
  }
}
table(df_hrgenes[df_hrgenes$type=="sensitive",2])
table(df_hrgenes[df_hrgenes$type=="resistant",2])
table(df_hrgenes[df_hrgenes$type=="sensitive",3])
table(df_hrgenes[df_hrgenes$type=="resistant",3])
t.test(df_hrgenes[df_hrgenes$type=="sensitive",2],df_hrgenes[df_hrgenes$type=="resistant",2])
colMeans(df_hrgenes[,2:ncol(df_hrgenes)]>0)
#      BRCA2       BRCA1       PALB2         ATM         NBN       CHEK2       BRIP1      MRE11A 
#0.034591195 0.037735849 0.006289308 0.022012579 0.003144654 0.003144654 0.000000000 0.003144654 
t.test(hrgenes_sensitive,hrgenes_resistant)
