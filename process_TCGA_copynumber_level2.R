#! /usr/bin/env Rscript
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
# geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
# 
# sum(! rownames(copynumber) %in% geneposition$gene)
# load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_copynumber.RData")
# sum(! colnames(data_aocs_copynumber) %in% geneposition$gene)
# hg19=read.table("/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/refgenefiles/hg19.txt",header=F,sep="\t")
# colnames(hg19)=c("chr1","start1","end1","gene")
# test=merge(hg19,geneposition)
# sum(! colnames(data_aocs_copynumber) %in% hg19[,4])
getmap_ID_file_TCGA=function(manifestfile,metadatafile)
{
  manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
  idxraw=which(grepl(".raw.copynumber.data.txt",manifest$filename))
  res1=data.frame(matrix(NA,nrow=length(idxraw),ncol=4))
  colnames(res1)=c("fileid","filename","sampleid","type")
  res1$fileid=manifest$id[idxraw]
  res1$filename=manifest$filename[idxraw]
  res1$type="raw"
  metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
  for (i in 1:nrow(res1))
  {
    fileid=res1$fileid[i]
    res1$sampleid[i]=get_sampleid(fileid,metadata)
  }
  
  idxtan=which(grepl(".tangent.copynumber.data.txt",manifest$filename))
  res2=data.frame(matrix(NA,nrow=length(idxtan),ncol=4))
  colnames(res2)=c("fileid","filename","sampleid","type")
  res2$fileid=manifest$id[idxtan]
  res2$filename=manifest$filename[idxtan]
  res2$type="tangent"
  for (i in 1:nrow(res2))
  {
    fileid=res2$fileid[i]
    res2$sampleid[i]=get_sampleid(fileid,metadata)
  }
  
  idxallele=which(grepl(".byallele.copynumber.data.txt",manifest$filename))
  res3=data.frame(matrix(NA,nrow=length(idxallele),ncol=4))
  colnames(res3)=c("fileid","filename","sampleid","type")
  res3$fileid=manifest$id[idxallele]
  res3$filename=manifest$filename[idxallele]
  res3$type="allele"
  for (i in 1:nrow(res3))
  {
    fileid=res3$fileid[i]
    res3$sampleid[i]=get_sampleid(fileid,metadata)
  }
  res=rbind(res1,res2,res3)
  return(res)
}

get_sampleid=function(fileid,metadata)
{
  res=NA
  idxrow=which(grepl(fileid,metadata[,1]))
  flag=F
  i=idxrow
  while(i<idxrow+100 & flag==F)
  {
    flag=grepl("entity_submitter_id",metadata[i,1])
    if (flag==T)
    {
      tmp=unlist(strsplit(metadata[i,1],": ",fixed=T))
      res=gsub(", ","",tmp[2])
    }
    i=i+1
  }
  return(res)
}
manifestfile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/gdc_manifest_20161007_182518.txt"
metadatafile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/metadata.cart.2016-10-07T18-07-28.399762.json"
filelist=getmap_ID_file_TCGA(manifestfile,metadatafile)
# write.table(filelist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/filelist.txt",
#             row.names=F,col.names=T,sep="\t",quote=F)

filelist=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/filelist.txt",header=T,sep="\t")

cn_anno=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.cn.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
cn_anno=cn_anno[,1:5]
cn_anno1=cbind.data.frame(cn_anno$Probe.Set.ID,cn_anno$Chromosome,floor(1/2*(as.numeric(cn_anno$Chromosome.Start)+as.numeric(cn_anno$Chromosome.Stop))))


snp_anno=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno=snp_anno[,1:5]
colnames(cn_anno1)=colnames(snp_anno)[c(1,3,4)]
snp_anno1=snp_anno[,c(1,3,4)]

anno=rbind(cn_anno1,snp_anno1)
keepcolnames=colnames(anno)
colnames(anno)[2:3]=c("chr","start")
anno1=sortgenetable(anno)
colnames(anno1)=keepcolnames

#remove probes in cnvs:
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
njob=100
require(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
mpi.bcast.Robj2slave(sortgenetable)
cnvfile="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/examplefiles/SNP6.merged.151117.hg19.CNV.txt"
cnvs=read.table(file=cnvfile,header=T,sep="\t")
cnvs$Chromosome=gsub("23","X",cnvs$Chromosome)
cnvs$Chromosome=gsub("24","Y",cnvs$Chromosome)
library(GenomicRanges)
mpi.bcast.cmd(library(GenomicRanges))
gr_cnvs=GRanges(seqnames=cnvs$Chromosome,ranges=IRanges(start=cnvs$Start,end=cnvs$End))
gr_anno1=GRanges(seqnames=anno1$Chromosome,ranges=IRanges(start=anno1$Physical.Position,width=1))
mpi.bcast.Robj2slave(gr_cnvs)
mpi.bcast.Robj2slave(gr_anno1)
probes_in_cnv=function(i)
{
  olap=subsetByOverlaps(gr_anno1[i,],gr_cnvs)
  res=F
  if (length(olap)>0) res=T
  return(res)
}
mpi.bcast.Robj2slave(probes_in_cnv)
idx_incnv=mpi.parSapply(X=1:length(gr_anno1),FUN=probes_in_cnv,job.num=njob)

anno1_rmcnv=anno1[!idx_incnv,]
save(anno1,anno1_rmcnv,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/snp6_anno.RData")

#check data
i=1115
i=512
myfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[i],"/",filelist$filename[i])
cnadata0=read.table(file=myfile,header=T,skip=1,sep="\t",stringsAsFactors=F)
cnadata=read.table(file=myfile,header=T,skip=1,sep="\t",stringsAsFactors=F)
colnames(cnadata)[1]=keepcolnames[1]
test=merge(cnadata,anno1)
test1=merge(anno1,cnadata)
#marker file:
marker=sortgenetable(test1)
marker=marker[,-4]
colnames(marker)=c('unitName','chromosome','position')
output="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA.markers.txt"
# write.table(marker,file=output,row.names=F,col.names=T,sep="\t",quote=F)


load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")
allsamples=colnames(copynumber)

mpi.bcast.Robj2slave(anno1)
mpi.bcast.Robj2slave(anno1_rmcnv)

myfiles=c()
count=1
for (mysample in allsamples)
{
  #idx=which(grepl(paste0(mysample,"-01"),filelist$sampleid) & filelist$type=="raw")
  idx=which(grepl(paste0(mysample,"-01"),filelist$sampleid) & filelist$type=="tangent")
  #a sample may have dulplicates with multiple plates, pick the one with higher number (https://confluence.broadinstitute.org/display/GDAC/FAQ)
  if (length(idx)>1) {
    tmp=as.character(filelist$sampleid[idx])
    tmp1=sapply(tmp,function(sn){
      res=as.numeric(unlist(strsplit(sn,"-",fixed=T))[6])
    })
    idx1=which(tmp1==max(tmp1))
    tmp=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idx[idx1]],"/",filelist$filename[idx[idx1]])
    myfiles=c(myfiles,tmp)
  }else
  {
    if (length(idx)==1)
    {
      tmp=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idx],"/",filelist$filename[idx])
      myfiles=c(myfiles,tmp)
    }
  }
  if (length(idx)==0) print(paste0("no: ",count))
  count=count+1
}
mpi.bcast.Robj2slave(myfiles)
mpi.bcast.Robj2slave(allsamples)

library(DNAcopy)
mpi.bcast.cmd(library(DNAcopy))


CBS=function(i,type="tangent",removecnv=F)
{
  chrs = paste0('chr',c(1:22,'X','Y'))
  myfile=myfiles[i]
  mysample=allsamples[i]
  cnadata=read.table(file=myfile,header=T,skip=1,sep="\t",stringsAsFactors=F)
  colnames(cnadata)=c("Probe.Set.ID","Log.R.Ratio")
  cnadata$Log.R.Ratio[which(cnadata$Log.R.Ratio==0)]=0.001
  cnadata$Log.R.Ratio=log2(cnadata$Log.R.Ratio/2) #CN to logR
  if (removecnv==T)
  {
    cnadata=merge(cnadata,anno1_rmcnv)
  }else
  {
    cnadata=merge(cnadata,anno1)
  }
  colnames(cnadata)=c("Probe.Set.ID","Log.R.Ratio","chr","start")
  cnadata=sortgenetable(cnadata)
  CNA.obj = CNA(cnadata$Log.R.Ratio, 
                cnadata$chr, 
                cnadata$start, 
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
      segs <- segment(CNAchr.smoothed,alpha=0.0001)
      segsall=rbind(segsall,segs$output)
    }
  }
  if (removecnv==T)
  {
    output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,"_",type,".tcga.rmcnv.segment.txt")
  }else
  {
    output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,"_",type,".tcga.segment.txt")
  }
  write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  return(nrow(cnadata))
}

mpi.bcast.Robj2slave(CBS)
res=mpi.parSapply(X=1:length(allsamples),FUN=CBS,job.num=njob)
#removecnv
res=mpi.parSapply(X=1:length(allsamples),removecnv=T,FUN=CBS,job.num=njob)

mpi.close.Rslaves()
mpi.quit()

type="raw"
#type="tangent"
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA.",type,".cbs.txt")
for (mysample in allsamples)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,"_",type,".tcga.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
segsall=read.table(file=output,header=T,sep="\t")

#form copynumber data based on segment result from DNAcopy:
copynumber_raw=readcnasegs1(hg19=T,segments=segsall)
colnames(copynumber_raw)=gsub(".","-",colnames(copynumber_raw),fixed=T)
alldata=rbind(data_copynumber$data1,data_copynumber$data2)
alldata1=formwholedataframe(alldata)
tcga_platinumclass=data.frame(sample=rownames(alldata1),platinumclass=alldata1[,"platinumclass"])
tcga_platinumclass$sample=as.character(tcga_platinumclass$sample)
data_copynumber_raw=cbind.data.frame(sample=colnames(copynumber_raw),t(copynumber_raw))
data_copynumber_raw[,"sample"]=as.character(data_copynumber_raw[,"sample"])
data_copynumber_raw=merge(tcga_platinumclass,data_copynumber_raw,by="sample")
rownames(data_copynumber_raw)=data_copynumber_raw[,"sample"]
data_copynumber_raw=data_copynumber_raw[,-which(colnames(data_copynumber_raw)=="sample")]

smallfdrgenes=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)
smallfdrgenes=smallfdrgenes$gene
for (i in 1:length(smallfdrgenes))
{
  idx=which(colnames(data_copynumber_raw)==smallfdrgenes[i])
  if (length(idx)>0)
  {
    print(i)
    fm=glm(as.factor(data_copynumber_raw[,"platinumclass"])~as.numeric(as.character(data_copynumber_raw[,idx])),family = "binomial")
    print(coef(summary(fm))[2,4])
  }
}

copynumber_raw_pvalues=sapply(2:ncol(data_copynumber_raw),function(i){
  fm=glm(as.factor(data_copynumber_raw[,"platinumclass"])~as.numeric(as.character(data_copynumber_raw[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_raw)[i-1]
  return(res)
})

draw_manhattan(pvalues=copynumber_raw_pvalues)

test=copynumber_raw_pvalues[!is.na(copynumber_raw_pvalues)]
draw_manhattan(test)

#work on tangent
type="tangent"
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA.",type,".cbs.txt")
for (mysample in allsamples)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,"_",type,".tcga.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
segsall=read.table(file=output,header=T,sep="\t")

#form copynumber data based on segment result from DNAcopy:
copynumber_tangent=readcnasegs1(hg19=T,segments=segsall)
colnames(copynumber_tangent)=gsub(".","-",colnames(copynumber_tangent),fixed=T)
data_copynumber_tangent=cbind.data.frame(sample=colnames(copynumber_tangent),t(copynumber_tangent))
data_copynumber_tangent[,"sample"]=as.character(data_copynumber_tangent[,"sample"])
data_copynumber_tangent=merge(tcga_platinumclass,data_copynumber_tangent,by="sample")
rownames(data_copynumber_tangent)=data_copynumber_tangent[,"sample"]
data_copynumber_tangent=data_copynumber_tangent[,-which(colnames(data_copynumber_tangent)=="sample")]

smallfdrgenes=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)
smallfdrgenes=smallfdrgenes$gene
for (i in 1:length(smallfdrgenes))
{
  idx=which(colnames(data_copynumber_tangent)==smallfdrgenes[i])
  if (length(idx)>0)
  {
    print(i)
    fm=glm(as.factor(data_copynumber_tangent[,"platinumclass"])~as.numeric(as.character(data_copynumber_tangent[,idx])),family = "binomial")
    print(coef(summary(fm))[2,4])
  }
}

copynumber_tangent_pvalues=sapply(2:ncol(data_copynumber_tangent),function(i){
  fm=glm(as.factor(data_copynumber_tangent[,"platinumclass"])~as.numeric(as.character(data_copynumber_tangent[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_tangent)[i-1]
  return(res)
})

test=copynumber_tangent_pvalues[!is.na(copynumber_tangent_pvalues)]
draw_manhattan(test)
for (i in 2:ncol(data_copynumber_raw))
{
  data_copynumber_raw[,i]=as.numeric(as.character(data_copynumber_raw[,i]))
}
for (i in 2:ncol(data_copynumber_tangent))
{
  data_copynumber_tangent[,i]=as.numeric(as.character(data_copynumber_tangent[,i]))
}

#removecnv
type="tangent"
segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA.",type,".rmcnv.cbs.txt")
for (mysample in allsamples)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,"_",type,".tcga.rmcnv.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
segsall=read.table(file=output,header=T,sep="\t")

#form copynumber data based on segment result from DNAcopy:
copynumber_tangent_rmcnv=readcnasegs1(hg19=T,segments=segsall)
colnames(copynumber_tangent_rmcnv)=gsub(".","-",colnames(copynumber_tangent_rmcnv),fixed=T)
data_copynumber_tangent_rmcnv=cbind.data.frame(sample=colnames(copynumber_tangent_rmcnv),t(copynumber_tangent_rmcnv))
data_copynumber_tangent_rmcnv[,"sample"]=as.character(data_copynumber_tangent_rmcnv[,"sample"])
data_copynumber_tangent_rmcnv=merge(tcga_platinumclass,data_copynumber_tangent_rmcnv,by="sample")
rownames(data_copynumber_tangent_rmcnv)=data_copynumber_tangent_rmcnv[,"sample"]
data_copynumber_tangent_rmcnv=data_copynumber_tangent_rmcnv[,-which(colnames(data_copynumber_tangent_rmcnv)=="sample")]


#save(data_copynumber_tangent,file="../data/platinum_classificationdata_stringent_tangent.RData")

#work from firehose segments:
segment_scna_minus_germline=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
                                         header=T,sep="\t",stringsAsFactors=F)
allsamples=colnames(copynumber)
segmentsamples=sapply(segment_scna_minus_germline$Sample,function(mysample){
  tmp=unlist(strsplit(mysample,"-",fixed=T))
  tmp1=substr(tmp[4],1,2)
  res=paste0(paste0(tmp[1:3],collapse="-"),"-",tmp1)
})

#check if there are duplicated tumor samples
segmentsamples1=unique(segment_scna_minus_germline$Sample)
count=1
for (mysample in allsamples)
{
  idx=which(grepl(paste0(mysample,"-01"),segmentsamples1))
  if (length(idx)>1)
  {
    print(count)
  }
  if (length(idx)==0)
  {
    print(paste0("no: ",count))
  }
  count=count+1
}
#no output, so there are no multiple tumor samples included in the segment
#pick tumor samples
idx=which(segmentsamples %in% paste0(allsamples,"-01"))
segment_scna_minus_germline=segment_scna_minus_germline[idx,]
segment_scna_minus_germline$Sample=segmentsamples[idx]
segment_scna_minus_germline$Chromosome=gsub("23","X",segment_scna_minus_germline$Chromosome)
#save segment only in tumor samples
#write.table(segment_scna_minus_germline,file="../result/firehose_segment.txt",sep="\t",quote=F,row.names = F,col.names = T)
copynumber_scna_minus_germline=readcnasegs1(hg19=T,segments=segment_scna_minus_germline)
data_copynumber_scna_minus_germline=cbind.data.frame(sample=colnames(copynumber_scna_minus_germline),t(copynumber_scna_minus_germline))
data_copynumber_scna_minus_germline[,"sample"]=as.character(data_copynumber_scna_minus_germline[,"sample"])
data_copynumber_scna_minus_germline[,"sample"]=gsub("-01$","",data_copynumber_scna_minus_germline[,"sample"],perl=T)
data_copynumber_scna_minus_germline=merge(tcga_platinumclass,data_copynumber_scna_minus_germline,by="sample")
rownames(data_copynumber_scna_minus_germline)=data_copynumber_scna_minus_germline[,"sample"]
data_copynumber_scna_minus_germline=data_copynumber_scna_minus_germline[,-which(colnames(data_copynumber_scna_minus_germline)=="sample")]
copynumber_scna_minus_germline_pvalues=sapply(2:ncol(data_copynumber_scna_minus_germline),function(i){
  fm=glm(as.factor(data_copynumber_scna_minus_germline[,"platinumclass"])~as.numeric(as.character(data_copynumber_scna_minus_germline[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_scna_minus_germline)[i-1]
  return(res)
})

test=copynumber_scna_minus_germline_pvalues[!is.na(copynumber_scna_minus_germline_pvalues)]
draw_manhattan(test)


#use the firehose segment:
copynumber_firehosesegment=readcnasegs(hg19=T,segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
data_copynumber_firehosesegment=cbind.data.frame(sample=colnames(copynumber_firehosesegment),t(copynumber_firehosesegment))
data_copynumber_firehosesegment=merge(tcga_platinumclass,data_copynumber_firehosesegment,by="sample")
data_copynumber_firehosesegment[,"sample"]=as.character(data_copynumber_firehosesegment[,"sample"])
data_copynumber_firehosesegment[,"sample"]=gsub("-01$","",data_copynumber_firehosesegment[,"sample"],perl=T)
rownames(data_copynumber_firehosesegment)=data_copynumber_firehosesegment[,"sample"]
data_copynumber_firehosesegment=data_copynumber_firehosesegment[,-which(colnames(data_copynumber_firehosesegment)=="sample")]
copynumber_firehosesegment_pvalues=sapply(2:ncol(data_copynumber_firehosesegment),function(i){
  fm=glm(as.factor(data_copynumber_firehosesegment[,"platinumclass"])~as.numeric(data_copynumber_firehosesegment[,i]),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_firehosesegment)[i-1]
  return(res)
})
test=copynumber_firehosesegment_pvalues[!is.na(copynumber_firehosesegment_pvalues)]
draw_manhattan(test)


#work on gistic data


#segment from firehose, redo gistic
copynumber_firehosesegment_gistic=read.table(file="../../Tools/GISTIC/TCGA_rx0_conf0.99_armpeel1_brlen0.7_broad1/all_data_by_genes.txt",header=T,
                                                 sep="\t",quote="")
rownames(copynumber_firehosesegment_gistic)=copynumber_firehosesegment_gistic$Gene.Symbol
copynumber_firehosesegment_gistic=copynumber_firehosesegment_gistic[,4:ncol(copynumber_firehosesegment_gistic)]
colnames(copynumber_firehosesegment_gistic)=gsub(".","-",colnames(copynumber_firehosesegment_gistic),fixed=T)
tmp=sapply(1:ncol(copynumber_firehosesegment_gistic),function(x){
  tmp1=unlist(strsplit(colnames(copynumber_firehosesegment_gistic)[x],"-",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
})
length(unique(tmp))
colnames(copynumber_firehosesegment_gistic)=tmp
data_copynumber_firehosesegment_gistic=cbind.data.frame(sample=colnames(copynumber_firehosesegment_gistic),t(copynumber_firehosesegment_gistic))
data_copynumber_firehosesegment_gistic[,"sample"]=as.character(data_copynumber_firehosesegment_gistic[,"sample"])
data_copynumber_firehosesegment_gistic=merge(tcga_platinumclass,data_copynumber_firehosesegment_gistic,by="sample")
rownames(data_copynumber_firehosesegment_gistic)=data_copynumber_firehosesegment_gistic[,"sample"]
data_copynumber_firehosesegment_gistic=data_copynumber_firehosesegment_gistic[,-which(colnames(data_copynumber_firehosesegment_gistic)=="sample")]

copynumber_firehosesegment_gistic_pvalues=sapply(2:ncol(data_copynumber_firehosesegment_gistic),function(i){
  fm=glm(as.factor(data_copynumber_firehosesegment_gistic[,"platinumclass"])~as.numeric(as.character(data_copynumber_firehosesegment_gistic[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_firehosesegment_gistic)[i-1]
  return(res)
})
test=copynumber_firehosesegment_gistic_pvalues[!is.na(copynumber_firehosesegment_gistic_pvalues)]
draw_manhattan(test)
#contains data from firehose
#save(data_copynumber,data_copynumber_firehosesegment_gistic,data_copynumber_firehosesegment,file="../data/platinum_classificationdata_stringent_firehose.RData")
copynumber_fdrs=read.table(file="../result/qvalues_copynumber_1000permutation.txt")
draw_manhattan(fdrs=copynumber_fdrs)


#the generated gistic data:
copynumber_tangent_gistic=read.table(file="../../Tools/GISTIC/TCGA_new_rx0_conf0.99_armpeel1_brlen0.7_broad1/all_data_by_genes.txt",header=T,
                                                 sep="\t",quote="")
rownames(copynumber_tangent_gistic)=copynumber_tangent_gistic$Gene.Symbol
copynumber_tangent_gistic=copynumber_tangent_gistic[,4:ncol(copynumber_tangent_gistic)]
colnames(copynumber_tangent_gistic)=gsub(".","-",colnames(copynumber_tangent_gistic),fixed=T)
tmp=sapply(1:ncol(copynumber_tangent_gistic),function(x){
  tmp1=unlist(strsplit(colnames(copynumber_tangent_gistic)[x],"-",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
})
length(unique(tmp))
data_copynumber_tangent_gistic=cbind.data.frame(sample=colnames(copynumber_tangent_gistic),t(copynumber_tangent_gistic))
data_copynumber_tangent_gistic[,"sample"]=as.character(data_copynumber_tangent_gistic[,"sample"])
data_copynumber_tangent_gistic=merge(tcga_platinumclass,data_copynumber_tangent_gistic,by="sample")
rownames(data_copynumber_tangent_gistic)=data_copynumber_tangent_gistic[,"sample"]
data_copynumber_tangent_gistic=data_copynumber_tangent_gistic[,-which(colnames(data_copynumber_tangent_gistic)=="sample")]

copynumber_tangent_gistic_pvalues=sapply(2:ncol(data_copynumber_tangent_gistic),function(i){
  fm=glm(as.factor(data_copynumber_tangent_gistic[,"platinumclass"])~as.numeric(as.character(data_copynumber_tangent_gistic[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_tangent_gistic)[i-1]
  return(res)
})
test=copynumber_tangent_gistic_pvalues[!is.na(copynumber_tangent_gistic_pvalues)]
draw_manhattan(test)
#save(data_copynumber_tangent,data_copynumber_tangent_rmcnv,data_copynumber_tangent_gistic,file="../data/platinum_classificationdata_stringent_tangent.RData")


#check with previous result:
copynumber_CGH1M=readcnasegs(hg19=F,segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.Level_3.2016012800.0.0/OV.cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.seg.txt")
data_copynumber_CGH1M=cbind.data.frame(sample=colnames(copynumber_CGH1M),t(copynumber_CGH1M))
data_copynumber_CGH1M=as.data.frame(data_copynumber_CGH1M)
data_copynumber_CGH1M[,"sample"]=as.character(data_copynumber_CGH1M[,"sample"])
data_copynumber_CGH1M[,"sample"]=gsub("-01$","",data_copynumber_CGH1M[,"sample"],perl=T)
data_copynumber_CGH1M=merge(tcga_platinumclass,data_copynumber_CGH1M,by="sample")
rownames(data_copynumber_CGH1M)=data_copynumber_CGH1M[,"sample"]
data_copynumber_CGH1M=data_copynumber_CGH1M[,-which(colnames(data_copynumber_CGH1M)=="sample")]
copynumber_CGH1M_pvalues=sapply(2:ncol(data_copynumber_CGH1M),function(i){
  fm=glm(as.factor(data_copynumber_CGH1M[,"platinumclass"])~as.numeric(as.character(data_copynumber_CGH1M[,i])),family = "binomial")
  res=NA
  if (nrow(coef(summary(fm)))>1)
  {
    res=coef(summary(fm))[2,4]
  }
  names(res)=colnames(data_copynumber_CGH1M)[i-1]
  return(res)
})

test=copynumber_CGH1M_pvalues[!is.na(copynumber_CGH1M_pvalues)]
draw_manhattan(test)



