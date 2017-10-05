#!/usr/bin/env Rscript
source("functions.R")

geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)

#first work on aocs data
load("../data/aocs_data.RData")
oacs_allsamples=allsamples
df_probes2keep=data.frame(Name=probes2keep)
aocs_anno=merge(df_probes2keep,omni_1.1)
aocs_anno=aocs_anno[,1:3]
colnames(aocs_anno)=c("probe","chr","start")
myfile=primarysamplefiles[1]
aocs_probedata=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                          stringsAsFactors=F)
aocs_probedata=aocs_probedata[,c("SNP.Name","Log.R.Ratio")]
colnames(aocs_probedata)=c("probe","value")
aocs_probedata=aocs_probedata[aocs_probedata$probe %in% aocs_anno$probe,]
aocs_probedata=merge(aocs_anno,aocs_probedata)
aocs_probedata=sortgenetable(aocs_probedata)
#salloc -t  1-1 -n 100 mpirun -n 1 R --interactive
njob=100
library(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
mpi.bcast.Robj2slave(aocs_probedata)
library(GenomicRanges)
mpi.bcast.cmd(library(GenomicRanges))
gr_genepostion=GRanges(seqnames=geneposition$chr,ranges=IRanges(start=geneposition$start,end=geneposition$end),gene=geneposition$gene)
mpi.bcast.Robj2slave(gr_genepostion)

map_probe_gene=function(i,probedata=aocs_probedata)
{
  res=NA
  gr_probe=GRanges(seqnames=probedata$chr[i],ranges=IRanges(start=probedata$start[i],width=1))
  olap=subsetByOverlaps(gr_genepostion,gr_probe)
  if (length(olap)>0)
  {
    res=mcols(olap)$gene[1]
  }
  return(res)
}
mpi.bcast.Robj2slave(map_probe_gene)
aocs_probe_genes=mpi.parSapply(X=1:nrow(aocs_probedata),probedata=aocs_probedata,FUN=map_probe_gene,job.num=njob)
#save(aocs_probe_genes,file="../data/aocs_level2data.Rdata")
#check the probes are the same for all the files
for (i in 20:25)
{
  myfile=primarysamplefiles[1]
  aocs_probedata=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                            stringsAsFactors=F)
  print(aocs_probedata$SNP.Name[1000000])
  print(aocs_probedata$SNP.Name[2000000])
  print(nrow(aocs_probedata))
}

mpi.bcast.Robj2slave(GSM_AOCS_ID)
mpi.bcast.Robj2slave(AOCS_GEO_ID)
mpi.bcast.Robj2slave(aocs_anno)
mpi.bcast.cmd(source("functions.R"))
form_level2data_fromaocsfile=function(myfile)
{
  gsm=unlist(strsplit(myfile,"_",fixed=T))[1]
  GSM_AOCS_ID$V1=gsub(" ","",GSM_AOCS_ID$V1,fixed=T)
  idx1=which(GSM_AOCS_ID$V1==gsm)
  aocs=GSM_AOCS_ID$V2[idx1]
  idx2=which(AOCS_GEO_ID$GEO_sample_id==aocs)
  #sample id:
  mysample=AOCS_GEO_ID$submitted_sample_id[idx2]
  aocs_probedata=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                            stringsAsFactors=F)
  aocs_probedata=aocs_probedata[,c("SNP.Name","Log.R.Ratio")]
  colnames(aocs_probedata)=c("probe","value")
  aocs_probedata=aocs_probedata[aocs_probedata$probe %in% aocs_anno$probe,]
  aocs_probedata=merge(aocs_anno,aocs_probedata)
  aocs_probedata=sortgenetable(aocs_probedata)
  res=data.frame(value=aocs_probedata$value)
  colnames(res)=mysample
  return(res)
}
mpi.bcast.Robj2slave(form_level2data_fromaocsfile)
num_aocsprobes=nrow(aocs_probedata)
res=mpi.parSapply(X=primarysamplefiles,FUN=form_level2data_fromaocsfile,job.num=njob)
aocs_level2data=data.frame(matrix(unlist(res),nrow=num_aocsprobes,byrow=F))
samplenames=rep(NA,length(res))
for (i in 1:length(res))
{
  tmp=unlist(strsplit(names(res[i]),".",fixed=T))
  samplenames[i]=tmp[length(tmp)]
}
colnames(aocs_level2data)=samplenames
rownames(aocs_level2data)=aocs_probedata$probe
#save(aocs_probe_genes,aocs_level2data,file="../data/aocs_level2data.Rdata")
#check the order:
sum(colnames(aocs_level2data)==aocs_platinumclass$sample)
data_aocs_level2data=cbind.data.frame(platinumclass=aocs_platinumclass$platinumclass,t(aocs_level2data))
#save(aocs_probe_genes,data_aocs_level2data,file="../data/aocs_level2data.Rdata")
y=data_aocs_level2data[,1]
mpi.bcast.Robj2slave(y)
data=data_aocs_level2data[,2:ncol(data_aocs_level2data)]
#the data is too big to use mpi.bcast, so split the jobs
pvalues_aocs_level2data=NULL
numprobes2run=10000
numrun <- ceiling(ncol(data)/numprobes2run)
print(paste0("total runs: ",numrun))
for (j in 1:numrun){
  print(j)
  if (j < numrun) cseq <- ((j-1)*numprobes2run+1):(j*numprobes2run)  else  cseq <- ((j-1)*numprobes2run+1):ncol(data)
  data1=data[,cseq]
  mpi.bcast.Robj2slave(data1)
  tmp=mpi_calpvalues(data=data1,y=y)
  pvalues_aocs_level2data=rbind.data.frame(pvalues_aocs_level2data,tmp)
}

#pvalues_aocs_level2data=mpi_calpvalues(data=data_aocs_level2data[,2:ncol(data_aocs_level2data)],y=data_aocs_level2data[,1])
df_pvalues_aocs_level2data=data.frame(gene=aocs_probe_genes,pvalue=pvalues_aocs_level2data)
#save(aocs_probe_genes,data_aocs_level2data,df_pvalues_aocs_level2data,file="../data/aocs_level2data.Rdata")

pvalues_gene_aocs_level2data=get_pvalues_gene(df_pvalues=df_pvalues_aocs_level2data)
save(aocs_probe_genes,data_aocs_level2data,df_pvalues_aocs_level2data,pvalues_gene_aocs_level2data,aocs_probe2gene_mean,file="../data/aocs_level2data.Rdata")

#work on tcga data
load("../data/tcga_data.RData")
tcga_allsamples=allsamples
tcga_anno=anno1
colnames(tcga_anno)=c("probe","chr","start")
#tangent probedata
myfiles=c()
count=1
for (mysample in tcga_allsamples)
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

#check probes for different samples
for (i in 1:1)
{
myfile=myfiles[i]
tcga_probedata=read.table(file=myfile,header=T,skip=1,sep="\t",stringsAsFactors=F)
colnames(tcga_probedata)=c("probe","value")
tcga_probedata=tcga_probedata[tcga_probedata$probe %in% tcga_anno$probe,]
tcga_probedata=merge(tcga_anno,tcga_probedata)
tcga_probedata=sortgenetable(tcga_probedata)
print(nrow(tcga_probedata))
print(tcga_probedata$probe[800000])
print(tcga_probedata$probe[1400000])
}
mpi.bcast.Robj2slave(tcga_probedata)
mpi.bcast.cmd(source("functions.R"))

tcga_probe_genes=mpi.parSapply(X=1:nrow(tcga_probedata),probedata=tcga_probedata,FUN=map_probe_gene,job.num=njob)
#save(tcga_probe_genes,file="../data/tcga_level2data.Rdata")

mpi.bcast.Robj2slave(tcga_anno)
mpi.bcast.Robj2slave(myfiles)
mpi.bcast.Robj2slave(tcga_allsamples)

form_level2data_fromtcgafile=function(myfile)
{
  tcga_probedata=read.table(file=myfile,header=T,skip=1,sep="\t",stringsAsFactors=F)
  colnames(tcga_probedata)=c("probe","value")
  tcga_probedata=tcga_probedata[tcga_probedata$probe %in% tcga_anno$probe,]
  tcga_probedata=merge(tcga_anno,tcga_probedata)
  tcga_probedata=sortgenetable(tcga_probedata)
  #tcga_probedata$probe=as.character(tcga_probedata$probe)
  idx=which(myfiles ==myfile)
  mysample=tcga_allsamples[idx]
  res=data.frame(value=tcga_probedata$value)
  colnames(res)=mysample
  return(res)
}
mpi.bcast.Robj2slave(form_level2data_fromtcgafile)
num_tcgaprobes=nrow(tcga_probedata)
res=mpi.parSapply(X=myfiles,FUN=form_level2data_fromtcgafile,job.num=njob)
tcga_level2data=data.frame(matrix(NA,nrow=num_tcgaprobes),ncol=length(myfiles))
for (i in 1:length(myfiles))
{
  tcga_level2data[,i]=res[[i]]
}
#tcga_level2data=data.frame(matrix(unlist(res),nrow=num_tcgaprobes,byrow=F))
tcga_samplenames=rep(NA,length(res))
for (i in 1:length(res))
{
  tmp=unlist(strsplit(names(res[i]),".",fixed=T))
  tcga_samplenames[i]=tmp[length(tmp)]
}
colnames(tcga_level2data)=tcga_samplenames
rownames(tcga_level2data)=tcga_probedata$probe
#save(tcga_probe_genes,tcga_level2data,file="../data/tcga_level2data.Rdata")
#check the order:
sum(colnames(tcga_level2data)==tcga_platinumclass$sample)
data_tcga_level2data=cbind.data.frame(platinumclass=tcga_platinumclass$platinumclass,t(tcga_level2data))
#save(tcga_probe_genes,data_tcga_level2data,file="../data/tcga_level2data.Rdata")
data=data_tcga_level2data[,2:ncol(data_tcga_level2data)]
y=data_tcga_level2data[,1]
mpi.bcast.Robj2slave(y)
#the data is too big to use mpi.bcast, so split the jobs
pvalues_tcga_level2data=NULL
numprobes2run=10000
numrun <- ceiling(ncol(data)/numprobes2run)
print(paste0("total runs: ",numrun))
for (j in 1:numrun){
  print(j)
  if (j < numrun) cseq <- ((j-1)*numprobes2run+1):(j*numprobes2run)  else  cseq <- ((j-1)*numprobes2run+1):ncol(data)
  data1=data[,cseq]
  mpi.bcast.Robj2slave(data1)
  tmp=mpi_calpvalues(data=data1,y=y)
  pvalues_tcga_level2data=rbind.data.frame(pvalues_tcga_level2data,tmp)
}

df_pvalues_tcga_level2data=data.frame(gene=tcga_probe_genes,pvalue=pvalues_tcga_level2data)
#save(tcga_probe_genes,data_tcga_level2data,df_pvalues_tcga_level2data,file="../data/tcga_level2data.Rdata")

get_pvalues_gene=function(df_pvalues=df_pvalues_tcga_level2data)
{
  df_pvalues=df_pvalues[! is.na(df_pvalues$gene),]
  df_pvalues$gene=as.character(df_pvalues$gene)
  genes=unique(df_pvalues$gene)
  res=data.frame(matrix(NA,nrow=length(genes),ncol=3))
  colnames(res)=c("gene","numprobes","pvalue")
  res$gene=genes
  count=1
  for (gene in genes)
  {
    idx=which(df_pvalues$gene==gene)
    res$numprobes[count]=length(idx)
    res$pvalue[count]=min(df_pvalues$pvalues[idx],na.rm=T)
    count=count+1
  }
  return(res)
}

pvalues_gene_tcga_level2data=get_pvalues_gene(df_pvalues=df_pvalues_tcga_level2data)
save(tcga_probe_genes,data_tcga_level2data,df_pvalues_tcga_level2data,pvalues_gene_tcga_level2data,tcga_probe2gene_mean,file="../data/tcga_level2data.Rdata")

load("../result/tmp_level2data.Rdata")
rownames(pvalues_gene_tcga_level2data)=pvalues_gene_tcga_level2data$gene
pvalues_gene_tcga_level2data$pvalue=as.numeric(pvalues_gene_tcga_level2data$pvalue)
draw_manhattan(pvalues=pvalues_gene_tcga_level2data,main="tcga_probe2gene")

rownames(pvalues_gene_aocs_level2data)=pvalues_gene_aocs_level2data$gene
pvalues_gene_aocs_level2data$pvalue=as.numeric(pvalues_gene_aocs_level2data$pvalue)
draw_manhattan(pvalues=pvalues_gene_aocs_level2data,main="aocs_probe2gene")

commongenes=intersect(pvalues_gene_tcga_level2data$gene,pvalues_gene_aocs_level2data$gene)
df_commongenes=data.frame(gene=commongenes)
pvalues_gene_level2data=merge(df_commongenes,pvalues_gene_tcga_level2data)
pvalues_gene_level2data=merge(pvalues_gene_level2data,pvalues_gene_aocs_level2data,by="gene")
colnames(pvalues_gene_level2data)=c("gene","numprobe_tcga","pvalue_tcga","numprobe_aocs","pvalue_aocs")
rownames(pvalues_gene_level2data)=pvalues_gene_level2data$gene
pvalues_gene_level2data=merge(pvalues_gene_level2data,geneposition,by="gene")
pvalues_gene_level2data=sortgenetable(pvalues_gene_level2data)
plot(pvalues_gene_level2data$numprobe_tcga,pvalues_gene_level2data$numprobe_aocs,xlab="tcga_numofprobes",ylab="aocs_numofprobes")
cor(pvalues_gene_level2data$numprobe_tcga,pvalues_gene_level2data$numprobe_aocs)
plot(pvalues_gene_level2data$pvalue_tcga,pvalues_gene_level2data$pvalue_aocs,xlab="tcga_pvalue",ylab="aocs_pvalue")
cor(pvalues_gene_level2data$pvalue_tcga,pvalues_gene_level2data$pvalue_aocs)


#find genes with small pvalue in both data
cutoff=0.001
idx=which(pvalues_gene_level2data$pvalue_tcga<=cutoff & pvalues_gene_level2data$pvalue_aocs<=cutoff)
View(pvalues_gene_level2data[idx,])
write.table(pvalues_gene_level2data[idx,],file="../result/p_valueslessthan0.001_level2data.txt",col.names=T,row.names=F,sep="\t",quote=F)

pvalues_small_tcga=data.frame(pvalue=pvalues_gene_level2data[idx,3])
rownames(pvalues_small_tcga)=rownames(pvalues_gene_level2data)[idx]
draw_manhattan(pvalues_small_tcga,main=paste0("tcga_p<=",cutoff))

pvalues_small_aocs=data.frame(pvalue=pvalues_gene_level2data[idx,5])
rownames(pvalues_small_aocs)=rownames(pvalues_gene_level2data)[idx]
draw_manhattan(pvalues_small_aocs,main=paste0("aocs_p<=",cutoff))



load("../data/tcga_level2data.Rdata")
which(rownames(data_tcga_level2data)=="TCGA-04-1331")
tmp_tcga_04_1331_level2data=t(data.frame(value=data_tcga_level2data[173,2:ncol(data_tcga_level2data)]))
tmp_tcga_04_1331_level2data=as.data.frame(tmp_tcga_04_1331_level2data)
rownames(tmp_tcga_04_1331_level2data)=colnames(data_tcga_level2data)[2:ncol(data_tcga_level2data)]
df_tcga_04_1331_level2data=data.frame(probe=rownames(tmp_tcga_04_1331_level2data),value=tmp_tcga_04_1331_level2data)
colnames(df_tcga_04_1331_level2data)[1]="Probe.Set.ID"
dim(df_tcga_04_1331_level2data)
df_tcga_04_1331_level2data=merge(df_tcga_04_1331_level2data,anno1)
dim(df_tcga_04_1331_level2data)
colnames(df_tcga_04_1331_level2data)=c("probe","value","chr","start")
df_tcga_04_1331_level2data=sortgenetable(df_tcga_04_1331_level2data)
df_tcga_04_1331_level2data$value=log2(df_tcga_04_1331_level2data$value/2)
tmp=draw_manhattan(pvalues=pvalues_copynumber_tangent,keepres = T)
df_tmp=cbind.data.frame(gene=rownames(tmp),tmp)
df_copynumber_tangent_tcga_041331=data.frame(gene=rownames(copynumber_tangent),value=copynumber_tangent[,1])
df_copynumber_tangent_tcga_041331=merge(df_tmp,df_copynumber_tangent_tcga_041331,by="gene")
idx1=which(df_tcga_04_1331_level2data$chr==1)
plot(df_tcga_04_1331_level2data$start[idx1],df_tcga_04_1331_level2data$value[idx1],col="red",xlab="coordinate",ylab="cn")
idx=which(df_copynumber_tangent_tcga_041331$chromosome==1)
points(df_copynumber_tangent_tcga_041331$position[idx],df_copynumber_tangent_tcga_041331$value[idx])
legend("bottomright",legend=c("raw","segmented"),col=c("red","black"),pch=1)
hist(pvalues_gene_tcga_level2data$numprobes[pvalues_gene_tcga_level2data$numprobes<=500],breaks=50)
quantile(pvalues_gene_tcga_level2data$numprobes)
df_pvalues_tcga_level2data=cbind.data.frame(Probe.Set.ID=rownames(df_pvalues_tcga_level2data),df_pvalues_tcga_level2data)
df_pvalues_tcga_level2data=merge(df_pvalues_tcga_level2data,anno1,by="Probe.Set.ID")
colnames(df_pvalues_tcga_level2data)=c("probe","gene","pvalue","chr","start")
#draw manhattan at probe level
draw_manhattan(pvalues=df_pvalues_tcga_level2data[,c("chr","start","pvalue")],main="tcga_probelevel")

df_pvalues_tcga_level2data=sortgenetable(df_pvalues_tcga_level2data)
idx=which(df_pvalues_tcga_level2data$gene=="C15orf32")
plot(df_pvalues_tcga_level2data$start[idx],-log10(as.numeric(df_pvalues_tcga_level2data$pvalue[idx])),xlab="coordinate",ylab="-log10(pvalue)",main="tcga_C15orf32")
df_pvalues_aocs_level2data=cbind.data.frame(probe=rownames(df_pvalues_aocs_level2data),df_pvalues_aocs_level2data)
df_pvalues_aocs_level2data=merge(df_pvalues_aocs_level2data,aocs_anno,by="probe")
colnames(df_pvalues_aocs_level2data)=c("probe","gene","pvalue","chr","start")
df_pvalues_aocs_level2data=sortgenetable(df_pvalues_aocs_level2data)
idx=which(df_pvalues_aocs_level2data$gene=="C15orf32")
plot(df_pvalues_aocs_level2data$start[idx],-log10(as.numeric(df_pvalues_aocs_level2data$pvalue[idx])),xlab="coordinate",ylab="-log10(pvalue)",main="aocs_C15orf32")

#draw manhattan at probe level
draw_manhattan(pvalues=df_pvalues_aocs_level2data[,c("chr","start","pvalue")],main="aocs_probelevel")

#save(df_pvalues_aocs_level2data,pvalues_gene_aocs_level2data,df_pvalues_tcga_level2data,pvalues_gene_tcga_level2data,pvalues_gene_level2data,file="../result/pvalues_level2data.Rdata")



idx=which(! smallfdr$gene %in% pvalues_gene_tcga_level2data$gene)
View(smallfdr[idx,])

#check the distribution of cn distribution among probes
gene="UNC5C"
gene="C15orf32"
idx=which(df_pvalues_tcga_level2data$gene==gene)
probe=as.character(df_pvalues_tcga_level2data$probe[idx])
start=df_pvalues_tcga_level2data$start[idx]
Mean_sens=Mean_resi=sd_sens=sd_resi=rep(NA,length(probe))
for (i in 1:length(probe))
{
  idx=which(colnames(data_tcga_level2data)==probe[i])
  Mean_sens[i]=mean(data_tcga_level2data[data_tcga_level2data[,1]=="Sensitive",idx],na.rm=T)
  Mean_resi[i]=mean(data_tcga_level2data[data_tcga_level2data[,1]=="Resistant",idx],na.rm=T)
  sd_sens[i]=sd(data_tcga_level2data[data_tcga_level2data[,1]=="Sensitive",idx],na.rm=T)
  sd_resi[i]=sd(data_tcga_level2data[data_tcga_level2data[,1]=="Resistant",idx],na.rm=T)
}

plot(start,Mean_sens,xlab="coordinate",ylab="cn",main=paste0("tcga_",gene),las=1,ylim=c(min(Mean_sens,Mean_resi)-0.1,max(Mean_sens,Mean_resi)+0.2))
points(start,Mean_resi,col="red")
legend("topleft",legend=c("Sensitive","Resistant"),col=c("black","red"),pch=1)
plot(start,sd_sens,xlab="coordinate",ylab="sd(cn)",main=paste0("tcga_",gene),las=1,ylim=c(min(sd_sens,sd_resi)-0.1,max(sd_sens,sd_resi)+0.2))
points(start,sd_resi,col="red")
legend("topleft",legend=c("Sensitive","Resistant"),col=c("black","red"),pch=1)

#load("../data/aocs_level2data.Rdata")

idx=which(df_pvalues_aocs_level2data$gene==gene)
which.min(df_pvalues_aocs_level2data$pvalue[idx])
#[1] 196
probe=as.character(df_pvalues_aocs_level2data$probe[idx])
start=df_pvalues_aocs_level2data$start[idx]
Mean_sens=Mean_resi=sd_sens=sd_resi=rep(NA,length(probe))
for (i in 1:length(probe))
{
  idx=which(colnames(data_aocs_level2data)==probe[i])
  cn_sens=2*2^data_aocs_level2data[data_aocs_level2data[,1]=="sensitive",idx]
  Mean_sens[i]=mean(cn_sens,na.rm=T)
  cn_resi=2*2^data_aocs_level2data[data_aocs_level2data[,1]=="resistant",idx]
  Mean_resi[i]=mean(cn_resi,na.rm=T)
  sd_sens[i]=sd(cn_sens,na.rm=T)
  sd_resi[i]=sd(cn_resi,na.rm=T)
}

plot(start,Mean_sens,xlab="coordinate",ylab="cn",main=paste0("aocs_",gene),las=1,ylim=c(min(Mean_sens,Mean_resi)-0.1,max(Mean_sens,Mean_resi)+0.2))
points(start,Mean_resi,col="red")
legend("topleft",legend=c("Sensitive","Resistant"),col=c("black","red"),pch=1)
plot(start,sd_sens,xlab="coordinate",ylab="sd(cn)",main=paste0("aocs_",gene),las=1,ylim=c(min(sd_sens,sd_resi)-0.1,max(sd_sens,sd_resi)+0.2))
points(start,sd_resi,col="red")
legend("topleft",legend=c("Sensitive","Resistant"),col=c("black","red"),pch=1)

#get gene level data using different methods
formgenedatafromprobe=function(probedata=data_tcga_level2data,pvalues=df_pvalues_tcga_level2data,opt="mean")
{
  pvalues$gene=as.character(pvalues$gene)
  genes=unique(pvalues$gene)
  genes=genes[! is.na(genes)]
  genes=genes[genes != "NA"]
  res=data.frame(matrix(NA,nrow=nrow(probedata),ncol=1+length(genes)))
  rownames(res)=rownames(probedata)
  colnames(res)=c("platinumclass",genes)
  res[,1]=probedata[,1]
  for (i in 2:ncol(res))
  {
    idx=which(pvalues$gene==colnames(res)[i])
    probes=pvalues$probe[idx]
    idx1=which(colnames(probedata) %in% probes)
    if (length(idx)==1)
    {
      res[,i]=probedata[,idx1]
    }else
    {
      if (opt=="mean")
      {
        res[,i]=rowMeans(probedata[,idx1],na.rm=T)
      }
    }
  }
  pvalues1=compute_pvalues1(res)
  draw_manhattan(pvalues=pvalues1)
  return(result=list(data=res,pvalues=pvalues1))
}

tcga_probe2gene_mean=formgenedatafromprobe(probedata=data_tcga_level2data,pvalues=df_pvalues_tcga_level2data,opt="mean")
aocs_probe2gene_mean=formgenedatafromprobe(probedata=data_aocs_level2data,pvalues=df_pvalues_aocs_level2data,opt="mean")

#work on those 158 genes:
load("../result/tmp_pvalues.RData")
smallfdr=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)
smallfdr1=cbind.data.frame(smallfdr,num_tcgaprobes=rep(0,nrow(smallfdr)),num_aocsprobes=rep(0,nrow(smallfdr)),mean_cp_aocs_sense=rep(NA,nrow(smallfdr)),mean_cp_aocs_resi=rep(NA,nrow(smallfdr)),
                           pvalues_tcga_cn=rep(NA,nrow(smallfdr)),pvalues_tcga_dh=rep(NA,nrow(smallfdr)),pvalues_tcga_2df=rep(NA,nrow(smallfdr)),
                           pvalues_aocs_cn=rep(NA,nrow(smallfdr)),pvalues_aocs_dh=rep(NA,nrow(smallfdr)),pvalues_aocs_2df=rep(NA,nrow(smallfdr)),
                          pvalue_tcga_probe=rep(NA,nrow(smallfdr)),pvalue_aocs_probe=rep(NA,nrow(smallfdr)))
for (i in 1:nrow(smallfdr))
{
  gene=smallfdr$gene[i]
  idx=which(pvalues_gene_tcga_level2data$gene==gene)
  if (length(idx)>0) smallfdr1$num_tcgaprobes[i]=pvalues_gene_tcga_level2data$numprobes[idx]
  idx=which(pvalues_gene_aocs_level2data$gene==gene)
  if (length(idx)>0) smallfdr1$num_aocsprobes[i]=pvalues_gene_aocs_level2data$numprobes[idx]
  idx=which(colnames(data_aocs_mean_copynumber)==gene)
  if (length(idx)>0)
  {
    cn_sens=2*2^data_aocs_mean_copynumber[data_aocs_mean_copynumber[,1]=="sensitive",idx]
    smallfdr1$mean_cp_aocs_sense[i]=mean(cn_sens,na.rm=T)
    cn_resi=2*2^data_aocs_mean_copynumber[data_aocs_mean_copynumber[,1]=="resistant",idx]
    smallfdr1$mean_cp_aocs_resi[i]=mean(cn_resi,na.rm=T)
  }
  idx=which(names(pvalues_copynumber_tangent)==gene)
  if (length(idx)>0) smallfdr1$pvalues_tcga_cn[i]=pvalues_copynumber_tangent[idx]
  idx=which(names(pvalues_copynumber_allele)==gene) #generated based on tangent data
  if (length(idx)>0) smallfdr1$pvalues_tcga_dh[i]=pvalues_copynumber_allele[idx]
  idx=which(names(pvalues_2df_tangent)==gene)
  if (length(idx)>0) smallfdr1$pvalues_tcga_2df[i]=pvalues_2df_tangent[idx]
  
  idx=which(names(pvalues_aocs_mean_copynumber)==gene)
  if (length(idx)>0) smallfdr1$pvalues_aocs_cn[i]=pvalues_aocs_mean_copynumber[idx]
  idx=which(names(pvalues_aocs_mean_copynumber_allele)==gene)
  if (length(idx)>0) smallfdr1$pvalues_aocs_dh[i]=pvalues_aocs_mean_copynumber_allele[idx]
  idx=which(names(pvalues_aocs_mean_2df)==gene)
  if (length(idx)>0) smallfdr1$pvalues_aocs_2df[i]=pvalues_aocs_mean_2df[idx]
  
  idx=which(names(tcga_probe2gene_mean$pvalues)==gene)
  if (length(idx)>0) smallfdr1$pvalue_tcga_probe[i]=tcga_probe2gene_mean$pvalues[idx]
  idx=which(names(aocs_probe2gene_mean$pvalues)==gene)
  if (length(idx)>0) smallfdr1$pvalue_aocs_probe[i]=aocs_probe2gene_mean$pvalues[idx]
  
}
smallfdr1$cp_mean_resi=2+smallfdr1$cp_mean_resi
smallfdr1$cp_mean_sens=2+smallfdr1$cp_mean_sens
write.table(smallfdr1,file="../result/copynumber_result_fdrlessthan0.05_update1.txt",col.names = T,row.names = F,sep="\t",quote=F)
smallfdr1=merge(smallfdr,pvalues_gene_level2data,by="gene")
colnames(smallfdr1)[3:5]=c("chr","start","end")
idx=which(colnames(smallfdr1) %in% c("chr.y","start.y","end.y"))
smallfdr1=smallfdr1[,-idx]
smallfdr1=sortgenetable(smallfdr1)
write.table(smallfdr1,file="../result/copynumber_result1_fdrlessthan0.05.txt",col.names=T,row.names=F,sep="\t",quote=F)

findcommonsmallpvalues=function(pvalues1=tcga_probe2gene_mean$pvalues,pvalues2=aocs_probe2gene_mean$pvalues)
{
  commongenes=intersect(names(pvalues1),names(pvalues2))
  df_pvalues1=data.frame(gene=names(pvalues1),pvalue1=pvalues1)
  df_pvalues2=data.frame(gene=names(pvalues2),pvalue2=pvalues2)
  df_pvalues=merge(df_pvalues1,df_pvalues2)
  maxofrow=apply(df_pvalues[,2:3],1,max)
  idx=order(maxofrow)
  cutoff=0.01
  print(paste0(sum(maxofrow<cutoff,na.rm=T)," genes with pvalue less than ", cutoff))
  res=df_pvalues[idx,]
}

tmp=findcommonsmallpvalues()
tmp1=tmp$pvalue1
names(tmp1)=tmp$gene
cutoff=0.01
idx=which(tmp$pvalue1<=cutoff & tmp$pvalue2<=cutoff)
draw_manhattan(pvalues = tmp1,main="tcga_cn")
draw_manhattan(pvalues = tmp1[idx],main=paste0("tcga_cn_with pvalue <=",cutoff))

tmp1=tmp$pvalue2
names(tmp1)=tmp$gene
cutoff=0.01
idx=which(tmp$pvalue1<=cutoff & tmp$pvalue2<=cutoff)
draw_manhattan(pvalues = tmp1,main="aocs_cn")
draw_manhattan(pvalues = tmp1[idx],main=paste0("aocs_cn_with pvalue <=",cutoff))

tmp=findcommonsmallpvalues(pvalues1=pvalues_copynumber_tangent,pvalues2=pvalues_aocs_mean_copynumber)
tmp=findcommonsmallpvalues(pvalues1=pvalues_copynumber,pvalues2=pvalues_aocs_copynumber)
tmp=findcommonsmallpvalues(pvalues1=pvalues_2df,pvalues2=pvalues_aocs_2df)
tmp=findcommonsmallpvalues(pvalues1=pvalues_2df_tangent,pvalues2=pvalues_aocs_mean_2df)
tmp=findcommonsmallpvalues(pvalues1=pvalues_copynumber_allele,pvalues2=pvalues_aocs_copynumber_allele)

tmp1=tmp$pvalue1
names(tmp1)=tmp$gene
cutoff=0.01
idx=which(tmp$pvalue1<=cutoff & tmp$pvalue2<=cutoff)
draw_manhattan(pvalues = tmp1,main="tcga_allele")
draw_manhattan(pvalues = tmp1[idx],main=paste0("tcga_allele_with pvalue <=",cutoff))

tmp1=tmp$pvalue2
names(tmp1)=tmp$gene
cutoff=0.01
idx=which(tmp$pvalue1<=cutoff & tmp$pvalue2<=cutoff)
draw_manhattan(pvalues = tmp1,main="aocs_allele")
draw_manhattan(pvalues = tmp1[idx],main=paste0("aocs_allele_with pvalue <=",cutoff))
