#!/usr/bin/env Rscript

segfile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/copy_number_somatic_mutation.tsv"
allsegs=read.table(file=segfile,header = T,sep="\t")
uniq_segsamples=as.character(unique(allsegs$submitted_sample_id))
uniq_segpatient=sapply(uniq_segsamples,function(mysample){
  res=unlist(strsplit(mysample,"-",fixed=T))
  res=paste0(res[1:2],collapse="-")
})
uniq_segpatient=unique(uniq_segpatient)

#get primary samples
samplefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/sample.OV-AU.1474317952273.tsv"
allsamples=read.table(file=samplefile,header=T,sep="\t",stringsAsFactors=F)
primarysamples=allsamples[allsamples$specimen_type=="Primary tumour - solid tissue",]
primarysamples=primarysamples[primarysamples$sequencing_strategy=="non-NGS",]
primarysamples=primarysamples[primarysamples$repository=="GEO",]
uniq_primarysamples=unique(primarysamples$submitted_sample_id)
uniq_primarypatient=sapply(uniq_primarysamples,function(mysample){
  res=unlist(strsplit(mysample,"-",fixed=T))
  res=paste0(res[1:2],collapse="-")
})

#find the primary samples in segfile
uniq_primarysegsample=uniq_segsamples[uniq_segsamples %in% uniq_primarysamples]
  

#work on outcome def
outcomedeffile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/nature14410-Supplementary Table 5.xlsx"
library(gdata)
outcomedef=read.xls(outcomedeffile,1,skip=1)
#the last two lines are comments
outcomedef=outcomedef[1:(nrow(outcomedef)-2),]
outcomedef$Chemotherapy.response=as.character(outcomedef$Chemotherapy.response)
outcomedef$Case.ID=as.character(outcomedef$Case.ID)
outcomedef$Case.ID=gsub("_","-",outcomedef$Case.ID,fixed=T)
outcomedef=outcomedef[outcomedef$Sample.type=="7PrimaryTumour",]
#form map of sample ->resistance def
sample_platinumstatus=data.frame(matrix(NA,nrow=nrow(outcomedef),ncol=3))
colnames(sample_platinumstatus)=c("patient","sample","status")
sample_platinumstatus$patient=outcomedef$Case.ID
sample_platinumstatus$status=outcomedef$Chemotherapy.response
#include samples in segfile
for (i in 1:nrow(sample_platinumstatus))
{
  idx=which(grepl(sample_platinumstatus$patient[i],uniq_primarysegsample)==T)
  if (length(idx)>0)
  {
    sample_platinumstatus$sample[i]=uniq_primarysegsample[idx]
  }
}

#extract segs from segfile only from primary tumors
usedsegs=allsegs[allsegs$submitted_sample_id %in% sample_platinumstatus$sample,]

#modified function readcnasegs() in readdata_firehose.R
readcnasegs1=function(hg19=T,segments,opt=1)
{
  library(CNTools)
  colnames(segments)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  dim(segments)
  length(unique(segments$ID))
  segments=segments[!is.na(segments$seg.mean),]
  cnseg_snp6=CNSeg(segments)
  #data("geneInfo"): it includes duplicate genes, need to be removed
  #load 
  # load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")

    #hg19=read.table("/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/refgenefiles/hg19.txt",header=F,sep="\t")
  if (opt==1)
  {
#     hg19_snp6=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/copynumber_geneposition_biomart.txt",header=T,sep="\t")
#     geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
#     hg19=merge(hg19_snp6,geneposition,by="gene")
#     hg19=hg19[,c(1,5,6,7)]
    #use the one from firehose
    hg19=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
  }
    if (opt==2)
    {
      load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")
      hg19=geneInfo_hg19[,c(4,1,2,3)]
    }
    colnames(hg19)=c("gene","chr","start","end")
    hg19=sortgenetable(hg19)
    
    hg19$chr=gsub("23","X",hg19$chr)
    hg19$chr=gsub("24","Y",hg19$chr)
    idxdup=duplicated(hg19$gene)
    hg19=hg19[! idxdup,]
    #change the column name of chr
    colnames(hg19)[2]="chrom"
    nstart=5
    snp6_max<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "max")
    snp6_max=rs(snp6_max)
    snp6_min<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "min")
    snp6_min=rs(snp6_min)
    snp6=snp6_max
  
  #pick the extreme one
  for (i in nstart:ncol(snp6))
  {
    pickmin=which(abs(snp6_min[,i])-abs(snp6_max[,i])>0)
    if (length(pickmin)>0) snp6[pickmin,i]=snp6_min[pickmin,i]
  }
  
  rownames(snp6)=snp6$gene
  snp6=snp6[,nstart:ncol(snp6)]
  #remove all 0 genes
  idxkeep=rep(T,nrow(snp6))
  idxkeep=sapply(1:nrow(snp6),function(i){
    res=F
    num0=sum(snp6[i,]==0)
    if (num0<ncol(snp6))
    {
      res=T
    }
    return(res)
  })
  snp6=snp6[idxkeep,]
  return(snp6)
}
#transform the segfile format to make it consistent with function readcnasegs() in readdata_firehose.R
transformseg=function(usedsegs,colnum=11)
{
  res=data.frame(matrix(NA,nrow=nrow(usedsegs),ncol=6))
  colnames(res)=c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")
  res$Sample=usedsegs$submitted_sample_id
  res$Chromosome=usedsegs$chromosome
  res$Start=usedsegs$chromosome_start
  res$End=usedsegs$chromosome_end
  res$Num_Probes=usedsegs$end_probe_id-usedsegs$start_probe_id+1
  res$Segment_Mean=usedsegs[,colnum]
  return(res)
}

transformedsegs=transformseg(usedsegs,colnum = 11)
aocs_copynumber=readcnasegs1(hg19=T,segments=transformedsegs)
aocs_platinumclass=rep(NA,ncol(aocs_copynumber))
for (i in 1:ncol(aocs_copynumber))
{
  idx=which(sample_platinumstatus$sample==colnames(aocs_copynumber)[i])
  if (length(idx)>0)
  {
    aocs_platinumclass[i]=sample_platinumstatus$status[idx]
  }
}
aocs_platinumclass=gsub("refractory","resistant",aocs_platinumclass)
#aocs_platinumclass=as.factor(aocs_platinumclass)
smallfdrgenes=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)
smallfdrgenes=smallfdrgenes$gene
for (i in 1:length(smallfdrgenes))
{
  idx=which(rownames(aocs_copynumber)==smallfdrgenes[i])
  if (length(idx)>0)
  {
    print(i)
    fm=glm(as.factor(aocs_platinumclass)~as.numeric(aocs_copynumber[idx,]),family = "binomial")
    print(coef(summary(fm))[2,4])
  }
}


aocs_copynumber_pvalues=sapply(1:nrow(aocs_copynumber),function(i){
  fm=glm(as.factor(aocs_platinumclass)~as.numeric(aocs_copynumber[i,]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=rownames(aocs_copynumber)[i]
  return(res)
})

draw_manhattan(pvalues=aocs_copynumber_pvalues)
platform="copynumber"
fdrs=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",platform,"_1000permutation.txt"),header=T)
draw_manhattan(fdrs=fdrs,fdrthreshold=0.05)

data_aocs_copynumber=cbind(platinumclass=aocs_platinumclass,as.data.frame(t(aocs_copynumber)))

#save(data_aocs_copynumber,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_copynumber.RData")

#check distribution of characteristics
sample_age=data.frame(matrix(NA,nrow=80,ncol=2))
colnames(sample_age)=c("patient","age")
sample_age[,1]=paste0("AOCS-",c("064","094","106","108","139","091","090","171","149","145",
                                "122","116","112","058","125","095","093","147","152","088",
                                "170","158","114","128","107","113","105","065","146","153",
                                "056","131","034","148","092","086","126","104","063","079",
                                "164","057","096","004","078","130","001","144","160","143",
                                "083","165","080","168","081","076","115","161","163","061",
                                "111","084","159","169","005","133","123","124","075","162",
                                "109","157","166","085","132","059","055","097","077","060"))
sample_age[,2]=c(67.5,59.0,64.6,61.2,62.4,39.2,53.8,52.1,56.6,62.2,
                 59.0,63.6,70.0,58.4,54.6,50.6,57.7,55.7,73.0,56.5,
                 72.8,52.5,55.3,45.1,47.9,77.4,73.9,46.3,52.1,57.9,
                 72.5,67.1,51.8,45.7,68.7,64.7,58.9,48.5,62.3,72.8,
                 55.1,56.1,61.6,53.8,59.9,65.6,53.8,54.4,66.1,49.7,
                 58.4,74.0,67.6,59.8,59.8,58.1,56.4,60.6,65.6,59.3,
                 65.8,60.9,78.4,74.1,51.6,66.8,55.8,58.1,54.0,74.7,
                 63.9,45.1,58.0,44.5,48.7,76.7,59.6,55.5,66.5,69.7)
sample_platinumstatus_age=merge(sample_platinumstatus,sample_age)
TCGA_clinical=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/clinical_vairables.txt",header=T,sep="\t")
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
TCGA_clinical=updateclinicalitem(TCGA_clinical)
TCGA_clinical=TCGA_clinical[is.na(TCGA_clinical[,"tumor_grade"]) | (! is.na(TCGA_clinical[,"tumor_grade"]) & TCGA_clinical[,"tumor_grade"]!="G1"),]
TCGA_clinical[,"tumor_grade"]=as.character(TCGA_clinical[,"tumor_grade"])
TCGA_clinical=TCGA_clinical[is.na(TCGA_clinical[,"clinical_stage"]) | (! is.na(TCGA_clinical[,"clinical_stage"]) & TCGA_clinical[,"clinical_stage"]!="Stage1"),]
TCGA_clinical[,"clinical_stage"]=as.character(TCGA_clinical[,"clinical_stage"])
t.test(TCGA_clinical[,"age"],sample_platinumstatus_age[,"age"])
t.test(TCGA_clinical[TCGA_clinical$platinumclass=="Sensitive","age"],sample_platinumstatus_age[sample_platinumstatus_age$status=="sensitive","age"])
t.test(TCGA_clinical[TCGA_clinical$platinumclass=="Resistant","age"],sample_platinumstatus_age[sample_platinumstatus_age$status %in% c("resistant","refractory"),"age"])
M=matrix(c(sum(TCGA_clinical$tumor_grade=="G2",na.rm=T),sum(TCGA_clinical$tumor_grade=="G3",na.rm=T),
           11,69),nrow=2,byrow=T)
fisher.test(M)
M=matrix(c(sum(TCGA_clinical$clinical_stage=="Stage3",na.rm=T),sum(TCGA_clinical$clinical_stage=="Stage4",na.rm=T),
           12+31+25,6+6),nrow=2,byrow=T)
fisher.test(M)
M=matrix(c(sum(TCGA_clinical$residual_disease_largest_nodule=="No_Macroscopic_disease",na.rm=T),sum(TCGA_clinical$residual_disease_largest_nodule %in% c("1_10mm","11_20mm","20mm_"),na.rm=T),
           1+2+1,11+35+30),nrow=2,byrow=T)
fisher.test(M)

#RNA-seq
tmp=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/sequencing_expression/exp_seq.tsv",header=T,sep="\t")
uniq_geneid=as.character(unique(tmp$gene_id))
df_uniq_geneid=data.frame(gene_id=uniq_geneid)
rnaseqsamples=allsamples$submitted_sample_id[allsamples$specimen_type=="Primary tumour - solid tissue" 
                                             & allsamples$sequencing_strategy=="RNA-Seq"]
rnaseqpatients=sapply(rnaseqsamples,function(mysample){
  res=unlist(strsplit(mysample,"-",fixed=T))
  res=paste0(res[1:2],collapse="-")
})
df_rnaseqpatients=data.frame(patient=rnaseqpatients,rnasample=rnaseqsamples)
sample_platinumstatus_rnaseq=merge(sample_platinumstatus,df_rnaseqpatients)
#1 sample without resistance def, removed
rnaseqsamples=as.character(sample_platinumstatus_rnaseq$rnasample)

aocs_rnaseq=data.frame(matrix(NA,nrow=length(uniq_geneid),ncol=length(rnaseqsamples)))
colnames(aocs_rnaseq)=rnaseqsamples
rownames(aocs_rnaseq)=uniq_geneid
for (i in 1:length(rnaseqsamples))
{
  mysample=rnaseqsamples[i]
  tmptable=tmp[tmp$submitted_sample_id==mysample,]
  if (nrow(tmptable)>0)
  {
    if (nrow(tmptable)==length(uniq_geneid))
    {
      tmptable1=merge(df_uniq_geneid,tmptable)
      aocs_rnaseq[,i]=tmptable1$normalized_read_count
      tmptable1=tmptable1[,c("gene_id","normalized_read_count")]
    }else
    {
      print(i)
    }
  }else
  {
    print(i)
  }
}

data_aocs_rnaseq=cbind(platinumclass=sample_platinumstatus_rnaseq$status,t(aocs_rnaseq))
data_aocs_rnaseq[,"platinumclass"]=gsub("refractory","resistant",data_aocs_rnaseq[,"platinumclass"])
data_aocs_rnaseq$platinumclass=as.factor(data_aocs_rnaseq$platinumclass)
save(data_aocs_copynumber,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_rnaseq.RData")
rnaseq_pvalues=sapply(2:ncol(data_aocs_rnaseq),function(i){
  res=NA
  if (sum(is.na(data_aocs_rnaseq[,i]))<40)
  {
    res=calpvalues(as.numeric(data_aocs_rnaseq[,i]),as.factor(data_aocs_rnaseq[,1]))  
  }
  return(res)
})

hist(rnaseq_pvalues,probability=T)

#process the raw data:
#ID convert:
GSM_AOCS_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSM_AOCS_ID.txt",sep="\t",header=F,stringsAsFactors = F)
AOCS_GEO_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSE65821_GEO_to_DCC_IDs.txt",header=T,sep="\t",stringsAsFactors = F)

#platforms:
omni_1.0=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni2.5-8v1_C_Gene_Annotation.txt",header=T,sep="\t",stringsAsFactors=F)
omni_1.1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni25M-8v1-1_B.annotated.txt",header=T,sep="\t",fill=T,stringsAsFactors=F)
commonprobes=intersect(omni_1.0$Name,omni_1.1$Name)
length(commonprobes)
#[1] 2338047
#find out probeset with GC score >0.7 across all the samples
njob=100
require(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
mpi.bcast.Robj2slave(commonprobes)
removedprobes=function(myfile)
{
  mytable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                     stringsAsFactors=F)
  print(nrow(mytable))
  mytable=mytable[mytable$SNP.Name %in%commonprobes,]
  mytable=mytable[! is.na(mytable$GC.Score),]
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".removedpobes.txt")
  res=mytable$SNP.Name[which(mytable$GC.Score<=0.7)]
  res1=data.frame(matrix(NA,nrow=length(res)))
  colnames(res1)="Name"
  res1[,1]=res
  write.table(res1,file=output,row.names=F,quote=F)
  return(0)
}
mpi.bcast.Robj2slave(removedprobes)
allfiles=list.files("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/")
myfiles=c()
for (mysample in uniq_primarysamples)
{
  idx=which(AOCS_GEO_ID$submitted_sample_id==mysample)
  if (length(idx)==1)
  {
    mygeo=AOCS_GEO_ID$GEO_sample_id[idx]
    idx1=which(GSM_AOCS_ID[,2]==mygeo)
    mygsm=GSM_AOCS_ID[idx1,1]
    mygsm=gsub(" ","",mygsm,fixed=T)
    idx2=which(grepl(paste0(mygsm,"\\w*.txt"),allfiles,perl=T)==T)
    myfile=allfiles[idx2]
    myfiles=c(myfiles,myfile)
  }
}
res=mpi.parSapply(X=myfiles,FUN=removedprobes,job.num=njob)

probes2remove=c()
for (myfile in myfiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".removedpobes.txt"),header=T)
  probes2remove=c(probes2remove,tmp$Name)
  probes2remove=unique(probes2remove)
}
length(probes2remove)
#[1] 538087
probes2keep=commonprobes[!commonprobes %in% probes2remove]
length(probes2keep)
save(probes2keep,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/probes2keep.RData")
#[1] 2338047
library(DNAcopy)
mpi.bcast.cmd(library(DNAcopy))
chrs = paste0('chr',c(1:22,'X','Y'))
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/probes2keep.RData")
mpi.bcast.Robj2slave(probes2keep)
mpi.bcast.Robj2slave(sortgenetable)
CBS=function(myfile)
{
  #to find sample id:
  GSM_AOCS_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSM_AOCS_ID.txt",sep="\t",header=F,stringsAsFactors = F)
  AOCS_GEO_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSE65821_GEO_to_DCC_IDs.txt",header=T,sep="\t",stringsAsFactors = F)
  gsm=unlist(strsplit(myfile,"_",fixed=T))[1]
  GSM_AOCS_ID$V1=gsub(" ","",GSM_AOCS_ID$V1,fixed=T)
  idx1=which(GSM_AOCS_ID$V1==gsm)
  aocs=GSM_AOCS_ID$V2[idx1]
  idx2=which(AOCS_GEO_ID$GEO_sample_id==aocs)
  mysample=AOCS_GEO_ID$submitted_sample_id[idx2]
  
  chrs = paste0('chr',c(1:22,'X','Y'))
  mytable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                     stringsAsFactors=F)
  mytable=mytable[mytable$SNP.Name %in% probes2keep,c("SNP.Name","Chr","Position" ,"B.Allele.Freq","Log.R.Ratio")]
  colnames(mytable)[2]="chr"
  colnames(mytable)[3]="start"
  mytable=sortgenetable(mytable)
  CNA.obj = CNA(mytable$Log.R.Ratio, 
                mytable$chr, 
                mytable$start, 
                data.type = "logratio",
                sampleid=mysample)
  smoothed.CNA.obj = smooth.CNA(CNA.obj)
  segsall=c()
  
  for (j in 1:length(chrs))
  {
    chr= chrs[j]
    if (grepl('chr',chr)) chr=substr(chr,4,nchar(chr))
    CNAchr.smoothed=smoothed.CNA.obj[smoothed.CNA.obj[,1]==chr,]
    if (nrow(CNAchr.smoothed)>0)
    {
      segs <- segment(CNAchr.smoothed)
      segsall=rbind(segsall,segs$output)
    }
  }
  output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.segment.txt")
  write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  return(nrow(mytable))
}
mpi.bcast.Robj2slave(CBS)
res=mpi.parSapply(X=myfiles,FUN=CBS,job.num=njob)

mpi.close.Rslaves()
mpi.quit()

segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.cbs.txt")
for (myfile in myfiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

segsall=read.table(file=output,header=T,sep="\t")
#form marker file
res=data.frame(matrix(NA,nrow=length(probes2keep),ncol=3))
colnames(res)=c('unitName','chromosome','position')
res[,1]=probes2keep
df_probes2keep=data.frame(Name=probes2keep)
tmp=merge(df_probes2keep,omni_1.0)
colnames(tmp)[2]="chr"
colnames(tmp)[3]="start"
res=sortgenetable(tmp)
res=res[,1:3]
colnames(res)=c('unitName','chromosome','position')
output="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.markers.txt"
write.table(res,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#form copynumber data based on segment result from DNAcopy:
aocs_mean_copynumber=readcnasegs1(hg19=T,segments=segsall)
colnames(aocs_mean_copynumber)=gsub(".","-",colnames(aocs_mean_copynumber),fixed=T)
aocs_mean_platinumclass=rep(NA,ncol(aocs_mean_copynumber))
for (i in 1:ncol(aocs_mean_copynumber))
{
  idx=which(sample_platinumstatus$sample==colnames(aocs_mean_copynumber)[i])
  if (length(idx)>0)
  {
    aocs_mean_platinumclass[i]=sample_platinumstatus$status[idx]
  }
}

aocs_mean_platinumclass=gsub("refractory","resistant",aocs_mean_platinumclass)
data_aocs_mean_copynumber=cbind(platinumclass=aocs_mean_platinumclass,as.data.frame(t(aocs_mean_copynumber)))

#compare copynumber data generated for 1 sample:
plot(as.numeric(data_aocs_copynumber[1,2:ncol(data_aocs_copynumber)]))
points(1:(ncol(data_aocs_copynumber)-1),as.numeric(data_aocs_mean_copynumber[1,2:ncol(data_aocs_copynumber)]),col="red")
legend("bottomleft",legend=c("median","mean"),col=c("black","red"),pch=1)

save(data_aocs_copynumber,data_aocs_mean_copynumber,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_copynumber.RData")
aocs_mean_platinumclass=as.factor(aocs_mean_platinumclass)
#aocs_platinumclass=as.factor(aocs_platinumclass)
smallfdrgenes=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)
smallfdrgenes=smallfdrgenes$gene
for (i in 1:length(smallfdrgenes))
{
  idx=which(rownames(aocs_mean_copynumber)==smallfdrgenes[i])
  if (length(idx)>0)
  {
    print(i)
    fm=glm(as.factor(aocs_mean_platinumclass)~as.numeric(aocs_mean_copynumber[idx,]),family = "binomial")
    print(coef(summary(fm))[2,4])
  }
}

aocs_mean_copynumber_pvalues=sapply(1:nrow(aocs_mean_copynumber),function(i){
  fm=glm(aocs_mean_platinumclass~as.numeric(aocs_mean_copynumber[i,]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=rownames(aocs_mean_copynumber)[i]
  return(res)
})

draw_manhattan(pvalues=aocs_mean_copynumber_pvalues)

hist(aocs_copynumber_pvalues)
hist(aocs_mean_copynumber_pvalues)

