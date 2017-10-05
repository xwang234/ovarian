#!/usr/bin/env Rscript

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")
updateclinicalitem=function(data=data_mrna)
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

sortgenedata=function(genedata)
{
  genedata$chr=as.character(genedata$chr)
  genedata$chr=gsub("23","X",genedata$chr)
  genedata$chr=gsub("24","Y",genedata$chr)
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genedata)))
  for (chr in chrs)
  {
    tmptable=genedata[which(genedata$chr==chr),]
    tmptable=tmptable[order(tmptable$start),]
    res=rbind(res,tmptable)
  }
  return(res)
}

extractqvaluesforgenesbylocation=function(qvaluetable,genetable=HR30genes,filename=NULL)
{
  genetable$Gene=as.character(genetable$Gene)
  genetable$Regions.targeted.for.capture...UCSC.hg19.=as.character(genetable$Regions.targeted.for.capture...UCSC.hg19.)
  qvaluetable=qvaluetable[!is.na(qvaluetable$chr),]
  gr_qvaluetable=GRanges(seqnames=qvaluetable$chr,ranges=IRanges(start=qvaluetable$start,end=qvaluetable$end),qvalue=qvaluetable$qvalues)
  
  res=data.frame(matrix(NA,nrow=nrow(genetable),ncol=5))
  colnames(res)=c("gene","chr","start","end","qvalue")
  res[,1]=genetable[,1]
  for (i in 1:nrow(genetable))
  {
    tmp=unlist(strsplit(genetable[i,2],":"))
    res[i,2]=gsub("chr","",tmp[1])
    tmp1=unlist(strsplit(tmp[2],"-",fixed=T))
    res[i,3]=as.integer(tmp1[1])
    res[i,4]=as.integer(tmp1[2])
    gr_res=GRanges(seqnames=res[i,2],ranges=IRanges(start=res[i,3],end=res[i,4]))
    olap=subsetByOverlaps(gr_qvaluetable,gr_res)
    if (length(olap)>0)
    {
      res[i,5]=min(mcols(olap)$qvalue)
    }
  }
  if (! is.null(filename))
  {
    write.table(res,file=filename,row.names=F,col.names=T,sep="\t",quote=F)
  }
  return(res)
}
gisticscore=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/scores.gistic",header=T,sep="\t")
gisticscore$Chromosome=gsub(23,"X",gisticscore$Chromosome)
library(GenomicRanges)
gr_gisticscore=GRanges(seqnames=gisticscore$Chromosome,IRanges(start=gisticscore$Start,end=gisticscore$End),type=gisticscore$Type,
                       gisticqvalue=10^-gisticscore$X.log10.q.value.)

HR30genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/HR30genes.txt",header=T,sep="\t")
#form gene position table considering genes from all the platforms
platforms=c("copynumber","copynumber_CGH1M","copynumber_1MDUO","mrna","mrna_u133","mrna_exon")
addtable=function(res,tmptable)
{
  if (nrow(tmptable)>0)
  {
    idxnew=which(! tmptable$gene %in% res$gene)
    if (length(idxnew)>0)
    {
      res=rbind(res,tmptable[idxnew,])
    }
  }
  return(res)
}
generate_geneposition=function(platforms)
{
  res=NULL
  for (platform in platforms)
  {
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_biomart.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    if (platform == platforms[1])
    {
      res=tmptable[idxnoNA,]
    }else
    {
      tmptable=tmptable[idxnoNA,]
      res=addtable(res,tmptable)
    }
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_gisticknowntable.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    tmptable=tmptable[idxnoNA,]
    res=addtable(res,tmptable)
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_ucscknowntable.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    tmptable=tmptable[idxnoNA,]
    res=addtable(res,tmptable)
  }
  return(res)
}

geneposition=generate_geneposition(platforms)
geneposition=sortgenedata(geneposition)
write.table(geneposition,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",row.names=F,col.names=T,sep="\t",quote=F)

checkgenes=function(geneqvaluefile,genepositionfile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",platform="copynumber")
{
  #qvaluetable=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_platinum.txt",header=T)
  qvaluetable=read.table(geneqvaluefile,header=T)
  #gene position, hg19
  #positiontable=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_geneposition_gisticknowntable.txt",header=T)
  positiontable=read.table(genepositionfile,header=T)
  #align genes in positiontable with those in qvaluetable
  idx=sapply(1:nrow(qvaluetable),function(i){
    res=which(positiontable$gene==rownames(qvaluetable)[i])
    if (length(res)==0)
    {
      res=NA
    }
    return(res)
  })
  #combine gene position and q value
  genenames=rownames(qvaluetable)
  qvaluetable=cbind(positiontable[idx,],qvaluetable)
  rownames(qvaluetable)=genenames
  
  idx_keep=qvaluetable$qvalues<=0.1
  if (sum(idx_keep)==0)
  {
    warning(paste0(platform," has no genes with small qvalues"))
    #select top 100 genes with smallest q values
    tmp=order(qvaluetable$qvalues)[100]
    idx_keep=qvaluetable$qvalues<=qvaluetable$qvalues[tmp]
  }
  selecttable=qvaluetable[idx_keep,]
  selecttable=sortgenedata(selecttable)
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
  filename1=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/",platform,"_","genes_smallqvalues.txt")
  fcon=file(filename1,"w")
  for (i in 1:nrow(selecttable))
  {
    writeLines(rownames(selecttable[i,]),fcon)
  }
  close(fcon)
  
  #prepare the data
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
    data=data_mrna
  }
  if (platform=="mrna_u133")
  {
    data=data_mrna_u133
  }
  #update stage and residual
  alldata=updateclinicalitem(data)
  #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
  alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
  alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
  alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  idx_y=which(colnames(alldata)=="platinumclass")
  #compute the correlation between two consecutive genes
  bet_cor=NA
  for (i in 2:nrow(selecttable))
  {
    idx1=which(colnames(alldata)==rownames(selecttable)[i-1])
    idx2=which(colnames(alldata)==rownames(selecttable)[i])
    tmp=cor(alldata[,idx1],alldata[,idx2],use="complete")
    bet_cor=c(bet_cor,tmp)
  }
  selecttable=cbind(selecttable,betweencorrelation=bet_cor)
  #compute the mean in two classes
  resistantdata=alldata[alldata[,"platinumclass"]=="Resistant",(idx_y+3):ncol(alldata)]
  sensitivedata=alldata[alldata[,"platinumclass"]=="Sensitive",(idx_y+3):ncol(alldata)]
  mean_res=data.frame(mean_resi=rep(NA,nrow(selecttable)),mean_sens=rep(NA,nrow(selecttable)),pvalue=rep(NA,nrow(selecttable)))
  for (i in 1:nrow(selecttable))
  {
    gene=rownames(selecttable)[i]
    idx=which(colnames(alldata)==gene)
    mean_res[i,1]=mean(resistantdata[,idx])
    mean_res[i,2]=mean(sensitivedata[,idx])
    mean_res[i,3]=t.test(resistantdata[,idx],sensitivedata[,idx])$p.value
    #info=paste0(selecttable$chr[i],":",selecttable$start[i],"-",selecttable$end[i])
    #hist(resistantdata[,idx],xlab="",main=paste0("resitant,",gene," ",info),breaks=20)
    #hist(sensitivedata[,idx],xlab="",main=paste0("sensitive,",gene," ",info),breaks=20)
  }
  selecttable=cbind(selecttable,mean_res[,c(1,2)])
  #search HR30genes
  filename2=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/HR30genesqvalues_",platform,".txt")
  qvalues_HR30genes=extractqvaluesforgenesbylocation(qvaluetable,genetable=HR30genes,filename=filename2)
  if (platform=="copynumber")
  {
    load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")
    cor_mrna_copynumber=data.frame(gene=selecttable$gene,correlation=rep(NA,nrow(selecttable)))
    for (i in 1:nrow(selecttable))
    {
      gene=rownames(selecttable)[i]
      idx=which(rownames(copynumber_comcm)==gene)
      if (length(idx)>0)
      {
        cor_mrna_copynumber[i,2]=cor(as.numeric(mrna_comcm[idx,]),as.numeric(copynumber_comcm[idx,]),use="complete")
      }
    }
    selecttable=cbind(selecttable,cor_mrna_copynumber=cor_mrna_copynumber[,2])
  }
  filename3=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneswithsmallqvaluestable_",platform,".txt")
  write.table(selecttable,file=filename3,row.names=F,col.names=T,sep="\t",quote=F)
  return(selecttable)
}

copynumber_res=checkgenes(geneqvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_platinum.txt",
                          platform="copynumber")

mrna_res=checkgenes(geneqvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_mrna_platinum.txt",
                          platform="mrna")

copynumber_CGH1M_res=checkgenes(geneqvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_CGH1M_platinum.txt",
                          platform="copynumber_CGH1M")

#check genes selected by glmnet
geneqvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_platinum.txt"
qvaluetable=read.table(geneqvaluefile,header=T)
#gene position, hg19
positiontable=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,stringsAsFactors=F)
positiontable=read.table(genepositionfile,header=T)
glmnettgenes=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_glmnetselectedgenes.txt",header=F)
glmnettgenes=as.character(glmnettgenes[2:nrow(glmnettgenes),1])
glmnetres=data.frame(matrix(NA,ncol=6,nrow=length(glmnettgenes)))
colnames(glmnetres)=c("gene","chr","start","end","pvalue","qvalue")
glmnetres[,"gene"]=glmnettgenes
for (i in 1:nrow(glmnetres))
{
  idx=which(positiontable$gene==glmnettgenes[i])
  if (length(idx)>0)
  {
    glmnetres[i,2]=as.character(positiontable[idx,2])
    glmnetres[i,3:4]=positiontable[idx,3:4]
  }
  idx=which(rownames(qvaluetable)==glmnettgenes[i])
  if (length(idx)>0)
  {
    glmnetres[i,5:6]=qvaluetable[idx,]
  }
}

glmnetres=sortgenedata(glmnetres)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_glmnetclassresult.RData")
glmnetcor=data.frame(matrix(NA,ncol=nrow(glmnetres),nrow=nrow(glmnetres)))
for (i in 1:(nrow(glmnetres)-1))
{
  gene1=glmnetres$gene[i]
  idx1=which(colnames(alldata2)==gene1)
  glmnetcor[i,i]=1
  for (j in (i+1):nrow(glmnetres))
  {
    gene2=glmnetres$gene[j]
    idx2=which(colnames(alldata2)==gene2)
    glmnetcor[j,i]=glmnetcor[i,j]=cor(alldata2[,idx1],alldata2[,idx2],use="complete")
  }
}
colnames(glmnetcor)=rownames(glmnetcor)=glmnetres$gene

#plot means of all genes on chromosomes
library(gap)
copynumbergenepos=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_geneposition_gisticknowntable.txt",header=T)
chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
idx=which(copynumbergenepos[,1] %in% chrs)
meanresis=colMeans(resistantdata)
meansens=colMeans(sensitivedata)
pvalues=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",header=T)
pvalues=pvalues[,1]
sortgenedata=function(genedata)
{
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genedata)))
  for (chr in chrs)
  {
    tmptable=genedata[which(genedata$chr==chr),]
    tmptable=tmptable[order(tmptable$start),]
    res=rbind(res,tmptable)
  }
  return(res)
}
plotchrdata=function(genepos=copynumbergenepos,data=meanresis,ylab="")
{
  data1=cbind(genepos[,c("chr","start")],data=data)
  data1$chr=as.character(data1$chr)
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  idx=which(genepos[,"chr"] %in% chrs)
  data1=data1[idx,]
  data1=sortgenedata(data1)
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
 plot(data1[,3])
}
plotchrdata(data=meansens)
plotchrdata(data=pvalues)




#compute the pvalues of cor_mrna_copynumber based on permutation
cor_mrna_copynumber_permutation=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/cor_mrna_copynumber_permutation.txt",header=T,sep="\t")
cor_mrna_copynumber_permutation_abs=abs(cor_mrna_copynumber_permutation)
cor_mrna_copynumber_pvalues=rep(NA,nrow(cor_mrna_copynumber))
for (i in 1:nrow(cor_mrna_copynumber))
{
  if (! is.na(cor_mrna_copynumber[i,"correlation"]))
  {
    cor_mrna_copynumber_pvalues[i]=sum(abs(cor_mrna_copynumber[i,"correlation"])<cor_mrna_copynumber_permutation_abs)/nrow(cor_mrna_copynumber_permutation)/ncol(cor_mrna_copynumber_permutation)
  }
}  
  
cor_mrna_copynumber=cbind(cor_mrna_copynumber,pvalues=cor_mrna_copynumber_pvalues)

#check if they are cancer genes
checkcancergenes=function(cancerlistfile,genes=as.character(smallfdr_copynumber$gene))
{
  cancergenes=read.table(file=cancerlistfile,header=F,sep="\t",fill=T,stringsAsFactors = F)
  cancergenelist=NULL
  for (i in 1:nrow(cancergenes))
  {
    cancergenelist=c(cancergenelist,cancergenes[i,])
  }
  cancergenelist=cancergenelist[!cancergenelist ==""]
  cancergenelist=gsub(" ","",cancergenelist,fixed = T)
  idx=which(genes %in% cancergenelist)
  res=genes[idx]
  return(res)
}

cancergenes=checkcancergenes(cancerlistfile="/fh/fast/dai_j/CancerGenomics/Tools/database/other/Cancergenes.txt")
candidatecancergenes=checkcancergenes(cancerlistfile="/fh/fast/dai_j/CancerGenomics/Tools/database/other/Cancercandidategenes.txt")
