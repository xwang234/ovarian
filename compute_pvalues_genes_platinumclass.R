#!/usr/bin/env Rscript

#calculate p-values using Rmpi
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
#salloc -t 2-1 -n 100 mpirun -n 1 ./compute_pvalues_genes_platinumclass.R
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
njob=100
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_4classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_4classification_otherdata_stringent.RData")


library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)

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
mpi.bcast.Robj2slave(calpvalues)

#parallel computing
mpi_calpvalues=function(data,y,filename=NULL,njob=100)
{
  
  res1 <- NULL
  nrun <- ceiling(ncol(data)/1000)
  for (j in 1:nrun){
    #cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):ncol(data)
    z=data[,cseq]
    res=mpi.parCapply(X=z,FUN=calpvalues,y=y,job.num=njob)
    res1=c(res1,res)
  }
  res1=res1[1:ncol(data)]
  res2=data.frame(pvalues=res1)
  rownames(res2)=names(res1)
  if (!is.null(filename))
  {
    write.table(res2,file=filename,quote=F,sep="\t",row.names=T,col.names=T)
  }
  return(res2)
}

calpvalues4=function(z,y)
{
  fit=glm(z~y)
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
mpi.bcast.Robj2slave(calpvalues4)

mpi4_calpvalues=function(data,y,filename=NULL,njob=100)
{
  
  res1 <- NULL
  nrun <- ceiling(ncol(data)/1000)
  for (j in 1:nrun){
    #cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):ncol(data)
    z=data[,cseq]
    res=mpi.parCapply(X=z,FUN=calpvalues4,y=y,job.num=njob)
    res1=c(res1,res)
  }
  res1=res1[1:ncol(data)]
  res2=data.frame(pvalues=res1)
  rownames(res2)=names(res1)
  if (!is.null(filename))
  {
    write.table(res2,file=filename,quote=F,sep="\t",row.names=T,col.names=T)
  }
  return(res2)
}

#update cilinical items of "clinical stage" and "residual disease"
updateclinicalitem1=function(data=data_mrna)
{
  if (is.data.frame(data)==F)
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
task="2class"
if (task=="2class")
{
  # compute p values on different data sets
  #platforms=c("copynumber","mrna","copynumber_CGH1M","copynumber_1MDUO","mrna_exon","mrna_u133")
  #platforms=c("copynumber","mrna","copynumber_1MDUO","copynumber_CGH1M","copynumber_snp6cna","copynumber_snp6cna_germline","mrna_exon","mrna_u133")
  #platforms=c("copynumber_extreme","mrna_extreme")
  platforms=c("copynumber_filtered","copynumber_tangent_filtered")
  platforms="copynumber_aocs_filtered"
  platforms=c("copynumber_tangent_gistic_filtered","copynumber_tangent_filtered")
  platforms=c("copynumber_firehosesegment","copynumber_gistic")
  platforms=c("copynumber_filtered","copynumber_tangent_filtered","copynumber_tangent_gistic_filtered","copynumber_aocs_filtered")
  for (platform in platforms)
  {
    print(platform)
    if (platform=="copynumber")
    {
      data=data_copynumber
    }
    if (platform=="copynumber_CGH1M")
    {
      data=data_copynumber_CGH1M
    }
    if (platform=="copynumber_1MDUO")
    {
      data=data_copynumber_1MDUO
    }
    if (platform=="mrna")
    {
      data=data_mrna
    }
    if (platform=="mrna_exon")
    {
      data=data_mrna_exon
    }
    if (platform=="mrna_u133")
    {
      data=data_mrna_u133
    }
    if (platform=="copynumber_snp6cna")
    {
      data=data_copynumber_snp6cna
    }
    if (platform=="copynumber_snp6cna_germline")
    {
      data=data_copynumber_snp6cna_germline
    }
    if (platform=="mrna_extreme")
    {
      data1=data4_mrna$data1
      data2=data4_mrna$data2
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive" & ! data1[,"platinumclass"] %in% "resistant",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive" & ! data2[,"platinumclass"] %in% "resistant",]
      data1[,"platinumclass"]=factor(data1[,"platinumclass"])
      data2[,"platinumclass"]=factor(data2[,"platinumclass"])
      data=list(data1=data1,data2=data2)
    }
    if (platform=="copynumber_extreme")
    {
      data1=data4_copynumber$data1
      data2=data4_copynumber$data2
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive" & ! data1[,"platinumclass"] %in% "resistant",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive" & ! data2[,"platinumclass"] %in% "resistant",]
      data1[,"platinumclass"]=factor(data1[,"platinumclass"])
      data2[,"platinumclass"]=factor(data2[,"platinumclass"])
      data=list(data1=data1,data2=data2)
    }
    if (platform=="copynumber_filtered")
    {
      #load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
      #use qthreshold=0.01 instead of 0.1
      load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.q0.01.RData")
      data=data_copynumber_filtered
    }
    if (platform=="copynumber_tangent_filtered")
    {
      #load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
      load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.q0.01.RData")
      data=data_copynumber_tangent_filtered
    }
    if (platform=="copynumber_aocs_filtered")
    {
      #load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
      load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.q0.01.RData")
      data=data_aocs_copynumber_filtered
    }
    if (platform=="copynumber_tangent_gistic_filtered")
    {
      #load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
      load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.q0.01.RData")
      data=data_copynumber_tangent_gistic_filtered
    }
    if (platform=="copynumber_firehosesegment")
    {
      load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_firehose.RData")
      data=data_copynumber_firehosesegment
    }
    if (platform=="copynumber_gistic")
    {
      load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_firehose.RData")
      data=data_copynumber_gistic
    }
    #outputprefix=paste0(platform,"_1000permutation")
    outputprefix=paste0(platform,"_1000permutation.q0.01")
    
    #update stage and residual
#     alldata=updateclinicalitem(data)
#     train=rep(T,nrow(alldata))
#     train[(nrow(data$data1)+1):nrow(alldata)]=F
#     alldata=cbind(alldata,train=train)
#     #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
#     alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
#     alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
#     alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
#     alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
#     #number of training
#     numtrain=sum(alldata[,"train"]==T)
#     #remove train column
#     alldata=alldata[,colnames(alldata)!="train"]
#     idx_y=which(colnames(alldata)=="platinumclass")
#     #outcome
#     y=alldata[,idx_y]
#     #gene matrix
#     data=alldata[,(idx_y+3):ncol(alldata)]
#     colnames(data)=gsub("-","_",colnames(data))
    alldata=formwholedataframe(data)
    idx_y=which(colnames(alldata)=="platinumclass")
    #     #outcome
    y=alldata[,idx_y]
    data=alldata[,(idx_y+1):ncol(alldata)]
    mpi.bcast.Robj2slave(data)
    #mpi.bcast.Robj2slave(y)
    
    
    #to compute permutations
    #permutation from iteration 2:
    permutationp=NULL
    numit=1000
    #numit=200
    for (i in 1:(numit+1))
    {
      y1=y
      if (i>1)
      {
        #permutation on outcome
        set.seed(i+1000)
        y1=y[sample(length(y))]
        names(y1)=names(y)
      }
      mpi.bcast.Robj2slave(y1)
      if (i==1)
      {
        #when i==1, calculate p values without permutation
        #pvalues_copynumber_platinum=mpi_calpvalues(data,y1,filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",njob=njob)
        filename1=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",platform,"_platinum.txt")
        pvalues_mrna_platinum=mpi_calpvalues(data,y1,filename=filename1,njob=njob)
      }else
      {
        #when i>1, calculate p values with permutation
        tmp=mpi_calpvalues(data,y1,njob=njob)
        #order the p values
        tmp1=tmp[order(tmp[,1]),1]
        tmp2=data.frame(pvalues=tmp1)
        rownames(tmp2)=paste0(rownames(tmp)[order(tmp[,1])],".",i-1)
        permutationp=rbind(permutationp,tmp2)
      }
      if (i %% 50==0)
      {
        cat(i,"..")
        #write.table(permutationp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_permutation.txt")
        filename2=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",outputprefix,".txt")
        write.table(permutationp,file=filename2)
      }
    }
    #write.table(permutationp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_permutation_top100.txt")
    #write.table(permutationp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_permutation.txt")
    write.table(permutationp,file=filename2)
    
    #to compute pvalues based on only training data
    #pvalues_copynumber_platinum_train=mpi_calpvalues(data[1:numtrain,],y[1:numtrain],filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_train.txt")
    #pvalues_mrna_platinum_train=mpi_calpvalues(data[1:numtrain,],y[1:numtrain],filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum_train.txt")
  }
}

if (task=="4class")
{
  # compute p values on different data sets
  #platforms=c("copynumber","mrna","copynumber_CGH1M","copynumber_1MDUO","mrna_exon","mrna_u133")
  #platforms=c("copynumber","mrna","copynumber_1MDUO","copynumber_CGH1M","copynumber_snp6cna","copynumber_snp6cna_germline","mrna_exon","mrna_u133")
  platforms=c("copynumber_rmborderline","mrna_rmborderline","copynumber_extreme","mrna_extreme") #remove borderline sensitive
  for (platform in platforms)
  {
    print(platform)
    if (platform=="copynumber")
    {
      data=data4_copynumber
    }
    if (platform=="copynumber_CGH1M")
    {
      data=data4_copynumber_CGH1M
    }
    if (platform=="copynumber_1MDUO")
    {
      data=data4_copynumber_1MDUO
    }
    if (platform=="mrna")
    {
      data=data4_mrna
    }
    if (platform=="mrna_exon")
    {
      data=data4_mrna_exon
    }
    if (platform=="mrna_u133")
    {
      data=data4_mrna_u133
    }
    if (platform=="copynumber_snp6cna")
    {
      data=data4_copynumber_snp6cna
    }
    if (platform=="copynumber_snp6cna_germline")
    {
      data=data4_copynumber_snp6cna_germline
    }
    if (platform=="copynumber_rmborderline")
    {
      data1=data4_copynumber$data1
      data2=data4_copynumber$data2
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive",]
      data=list(data1=data1,data2=data2)
    }
    if (platform=="mrna_rmborderline")
    {
      data1=data4_mrna$data1
      data2=data4_mrna$data2
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive",]
      data=list(data1=data1,data2=data2)
    }
    if (platform=="copynumber_extreme")
    {
      data1=data4_copynumber$data1
      data2=data4_copynumber$data2
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive" & ! data1[,"platinumclass"] %in% "resistant",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive" & ! data2[,"platinumclass"] %in% "resistant",]
      data=list(data1=data1,data2=data2)
    }
    if (platform=="mrna_extreme")
    {
      data1=data4_mrna$data1
      data2=data4_mrna$data2
      data1=data1[! data1[,"platinumclass"] %in% "borderlinesensitive" & ! data1[,"platinumclass"] %in% "resistant",]
      data2=data2[! data2[,"platinumclass"] %in% "borderlinesensitive" & ! data2[,"platinumclass"] %in% "resistant",]
      data=list(data1=data1,data2=data2)
    }
    outputprefix=paste0(platform,"_4class_1000permutation")
    
    #update stage and residual
    alldata=updateclinicalitem(data)
    train=rep(T,nrow(alldata))
    train[(nrow(data$data1)+1):nrow(alldata)]=F
    alldata=cbind(alldata,train=train)
    #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
    alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
    alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
    alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
    alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
    #number of training
    numtrain=sum(alldata[,"train"]==T)
    #remove train column
    alldata=alldata[,colnames(alldata)!="train"]
    idx_y=which(colnames(alldata)=="platinumclass")
    #outcome
    y=alldata[,idx_y]
    y[y=="refractory"]=1
    y[y=="resistant"]=2
    y[y=="borderlinesensitive"]=3
    y[y=="sensitive"]=4
    y=as.integer(y)
    #gene matrix
    data=alldata[,(idx_y+3):ncol(alldata)]
    colnames(data)=gsub("-","_",colnames(data))
    mpi.bcast.Robj2slave(data)
    #mpi.bcast.Robj2slave(y)
    
    #to compute permutations
    #permutation from iteration 2:
    permutationp=NULL
    numit=200
    for (i in 1:(numit+1))
    {
      y1=y
      if (i>1)
      {
        #permutation on outcome
        set.seed(i+1000)
        y1=y[sample(length(y))]
      }
      mpi.bcast.Robj2slave(y1)
      if (i==1)
      {
        #when i==1, calculate p values without permutation
        #pvalues_copynumber_platinum=mpi_calpvalues(data,y1,filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",njob=njob)
        filename1=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_4class_",platform,"_platinum.txt")
        pvalues_mrna_platinum=mpi4_calpvalues(data,y1,filename=filename1,njob=njob)
      }else
      {
        #when i>1, calculate p values with permutation
        tmp=mpi4_calpvalues(data,y1,njob=njob)
        #order the p values
        tmp1=tmp[order(tmp[,1]),1]
        tmp2=data.frame(pvalues=tmp1)
        rownames(tmp2)=paste0(rownames(tmp)[order(tmp[,1])],".",i-1)
        permutationp=rbind(permutationp,tmp2)
      }
      if (i %% 50==0)
      {
        cat(i,"..")
        #write.table(permutationp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_permutation.txt")
        filename2=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_",outputprefix,".txt")
        write.table(permutationp,file=filename2)
      }
    }
    #write.table(permutationp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_permutation_top100.txt")
    #write.table(permutationp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_permutation.txt")
    write.table(permutationp,file=filename2)
    
    #to compute pvalues based on only training data
    #pvalues_copynumber_platinum_train=mpi_calpvalues(data[1:numtrain,],y[1:numtrain],filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum_train.txt")
    #pvalues_mrna_platinum_train=mpi_calpvalues(data[1:numtrain,],y[1:numtrain],filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum_train.txt")
  }
}


print("done")
mpi.close.Rslaves()
#program ends here
mpi.quit()

#pvalues_copynumber_platinum=mpi_calpvalues(data,y,filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt")
#pvalues_mrna_platinum=mpi_calpvalues(data,y,filename="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum.txt")

#pvalues_copynumber_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
#pvalues_mrna_platinum=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_platinum.txt",header=T)
