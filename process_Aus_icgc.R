#! /usr/bin/env Rscript
library(data.table)
source("functions.R")
load("../data/aocs_data.RData")

#find germline vcf from sequencing file
allsamples=rownames(data_aocs_copynumber)
germlinevcffiles=rep(NA,length(allsamples))
donortable=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/clinical/donor.tsv",header=T,sep="\t",stringsAsFactors = F)
manifesttable=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/manifest.collaboratory.1480527008951.tsv",header=T,sep="\t",stringsAsFactors = F)
manifestable_london=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/manifest.pcawg-london.1480448530358.xml",header=F,sep="|",stringsAsFactors = F)
idx_analysis=which(grepl("<analysis_id>",manifestable_london[,1])==T)
for (i in 1:length(allsamples))
{
  cat(i,"..")
  mysample=allsamples[i]
  mysample1=unlist(strsplit(mysample,"-",fixed=T))
  mysample1=paste0(mysample1[1:2],collapse = "-")
  idx=which(donortable$submitted_donor_id==mysample1)
  mydonor=donortable$icgc_donor_id[idx]
  idx1=which(manifesttable$donor_id.donor_count==mydonor)
  for (j in idx1)
  {
    md5=manifesttable$md5_sum[j]
    idx2=which(grepl(md5,manifestable_london[,1])==T)
    if (length(idx2)>0)
    {
      myfile=manifestable_london[idx2-2,1]
      myfile=gsub(" ","",myfile)
      myfile=gsub("<filename>","",myfile)
      myfile=gsub("</filename>","",myfile)
      if (grepl("germline.snv_mnv.vcf",myfile)==T)
      {
        idx3=which(idx_analysis<idx2)
        idx3=idx_analysis[length(idx3)]
        myfolder=manifestable_london[idx3,1]
        myfolder=gsub(" ","",myfolder)
        myfolder=gsub("<analysis_id>","",myfolder)
        myfolder=gsub("</analysis_id>","",myfolder)
        cnvfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/",myfolder,"/",myfile)
        cnvfile=gsub(".gz$","",cnvfile)
        if (!file.exists(cnvfile))
          system(paste0("gunzip ",paste0(cnvfile),".gz"))
        
        germlinevcffiles[i]=cnvfile
        break
      }
    }
  }
}
names(germlinevcffiles)=allsamples

samplefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/sample.OV-AU.1474317952273.tsv"
allsamples=read.table(file=samplefile,header=T,sep="\t",stringsAsFactors=F)
primarysamples=allsamples[allsamples$specimen_type=="Primary tumour - solid tissue",]
primarysamples=primarysamples[primarysamples$sequencing_strategy=="non-NGS",]
primarysamples=primarysamples[primarysamples$repository=="GEO",]
uniq_primarysamples=unique(primarysamples$submitted_sample_id)

allfiles=list.files("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/")
primarysamplefiles=c()
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
    primarysamplefiles=c(primarysamplefiles,myfile)
  }
}

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/probes2keep.RData")
omni_1.1=fread(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni25M-8v1-1_B.annotated.txt",header=T,sep="\t",fill=T,stringsAsFactors=F)
omni_1.1=as.data.frame(omni_1.1)
df_probes2keep=data.frame(Name=probes2keep)
aocs_anno=merge(df_probes2keep,omni_1.1)
aocs_anno=aocs_anno[,1:3]
colnames(aocs_anno)=c("probe","chr","start")
aocs_anno$probe=as.character(aocs_anno$probe)

findhtprobesinnormal=function(data1=data1,data2=data2)
{
  tmp=merge(data1,data2,by=c("chr","start"))
  gtflag=sapply(1:nrow(tmp),function(i){
    res1=F
    res=unlist(strsplit(tmp$control[i],":",fixed=T))[1]
    if (res=="0/1" | res=="1/0") res1=T
    return(res1)
  })
  leftprobes=tmp$probe.y[gtflag]
}

#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
njob=100
require(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
mpi.bcast.Robj2slave(probes2keep)
mpi.bcast.Robj2slave(probes2keep_rmcnv)
mpi.bcast.Robj2slave(aocs_anno)
mpi.bcast.Robj2slave(findhtprobesinnormal)
mpi.bcast.Robj2slave(germlinevcffiles)
library(DNAcopy)
mpi.bcast.cmd(library(DNAcopy))
#parallel computing part:
#function to apply CBS on each level 2 copy ratio file:
#To see how it works, you can run CBS(myfile=primarysamplefiles[1])
CBS=function(myfile,type="tcn",removecnv=F)
{
  #to find sample id:
  GSM_AOCS_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSM_AOCS_ID.txt",sep="\t",header=F,stringsAsFactors = F)
  AOCS_GEO_ID=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/GSE65821_GEO_to_DCC_IDs.txt",header=T,sep="\t",stringsAsFactors = F)
  gsm=unlist(strsplit(myfile,"_",fixed=T))[1]
  GSM_AOCS_ID$V1=gsub(" ","",GSM_AOCS_ID$V1,fixed=T)
  idx1=which(GSM_AOCS_ID$V1==gsm)
  aocs=GSM_AOCS_ID$V2[idx1]
  idx2=which(AOCS_GEO_ID$GEO_sample_id==aocs)
  #sample id:
  mysample=AOCS_GEO_ID$submitted_sample_id[idx2]
  
  chrs = paste0('chr',c(1:22,'X','Y'))
  #read the log ratio data:
  mytable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                     stringsAsFactors=F)
  #Apply CBS on the data: only keep a few columns, and probes after filtering
  if (removecnv==T)
  {
    mytable=mytable[mytable$SNP.Name %in% probes2keep_rmcnv,c("SNP.Name","Chr","Position" ,"B.Allele.Freq","Log.R.Ratio")]
  }else
  {
    mytable=mytable[mytable$SNP.Name %in% probes2keep,c("SNP.Name","Chr","Position" ,"B.Allele.Freq","Log.R.Ratio")]
  }
  
  if (type=="tcn") #total copy number
  {
    CNA.obj = CNA(mytable$Log.R.Ratio, 
                  mytable$Chr, 
                  mytable$Position, 
                  data.type = "logratio",
                  sampleid=mysample)
  }else #DH
  {
    idx3=which(names(germlinevcffiles)==mysample)
    if (length(idx3)==0)
    {
      cutoff=0.2
      mytable=mytable[mytable$B.Allele.Freq>=cutoff & mytable$B.Allele.Freq<=1-cutoff,]
    }else
    {
      germtable=read.table(file=germlinevcffiles[idx3],skip=42,header=T,sep="\t",comment.char = "",stringsAsFactors = F)
      #overlap by position
      data1=germtable[,c(1,2,3,7,10)]
      data1=data1[data1$FILTER=="PASS",]
      data1=data1[,-4]
      colnames(data1)=c("chr","start","probe","control")
      leftprobes_normal=findhtprobesinnormal(data1=data1,data2=aocs_anno)
      mytable=mytable[mytable$SNP.Name %in% leftprobes_normal,]
      
    }
    
    CNA.obj = CNA(abs(0.5-mytable$B.Allele.Freq), 
                  mytable$Chr, 
                  mytable$Position, 
                  data.type = "logratio",
                  sampleid=mysample)
    
  }
  smoothed.CNA.obj = smooth.CNA(CNA.obj)
  
  segsall=c()
  #segment on each chromsome
  for (j in 1:length(chrs))
  {
    chr= chrs[j]
    if (grepl('chr',chr)) chr=substr(chr,4,nchar(chr))
    CNAchr.smoothed=smoothed.CNA.obj[smoothed.CNA.obj[,1]==chr,]
    if (nrow(CNAchr.smoothed)>0)
    {
      #Do the segment:
      if (type=="tcn")
      {
        segs <- segment(CNAchr.smoothed,alpha=0.0001)
      }else
      {
        segs <- segment(CNAchr.smoothed,alpha=0.05)
      }
      segsall=rbind(segsall,segs$output)
    }
  }
  #write segment results to tmp folder
  if (removecnv==T)
  {
    if (type=="tcn")
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.rmcnv.segment.txt")
    }else
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.rmcnv.segment.txt")
    }
  }else
  {
    if (type=="tcn")
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".cbs.segment.txt")
    }else
    {
      output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.segment.txt")
      dh_probe_file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.txt")
      write.table(mytable,file=dh_probe_file,row.names=F,col.names=T,sep="\t",quote=F)
    }
  }
  write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  
  return(nrow(mytable))
}

mpi.bcast.Robj2slave(CBS)
res=mpi.parSapply(X=primarysamplefiles,type="allele",FUN=CBS,job.num=njob)

segsall=c()
output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.allele.cbs.txt")
for (myfile in primarysamplefiles)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.cbs.segment.txt"),header=T)
  segsall=rbind(segsall,tmp)
}
#write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)

#once the segment file was written, we can start from here by simply reading the above results.
segsall=read.table(file=output,header=T,sep="\t")
aocs_mean_copynumber_allele=readcnasegs1(segments=segsall)
colnames(aocs_mean_copynumber_allele)=gsub(".","-",colnames(aocs_mean_copynumber_allele),fixed=T)
sum(colnames(aocs_mean_copynumber_allele)==aocs_platinumclass$sample)
#combine data and output:
data_aocs_mean_copynumber_allele=cbind(platinumclass=aocs_platinumclass$platinumclass,as.data.frame(t(aocs_mean_copynumber_allele)))




#find the germline info from sequencing data:
df_probes2keep=data.frame(Name=probes2keep)
aocs_anno=merge(df_probes2keep,omni_1.1)
aocs_anno=aocs_anno[,1:3]
colnames(aocs_anno)=c("probe","chr","start")
aocs_anno$probe=as.character(aocs_anno$probe)
#get the germline data:
sampleid=1
res2=c()
for (sampleid in 1:10)
{
mysample=allsamples[sampleid]

germtable=read.table(file=germlinevcffiles[sampleid],skip=42,header=T,sep="\t",comment.char = "",stringsAsFactors = F)
#overlap by position
data1=germtable[,c(1,2,3,7,10)]
data1=data1[data1$FILTER=="PASS",]
data1=data1[,-4]
colnames(data1)=c("chr","start","probe","control")
data2=aocs_anno
findhtprobesinnormal=function(data1=data1,data2=data2)
{
  tmp=merge(data1,data2,by=c("chr","start"))
  gtflag=sapply(1:nrow(tmp),function(i){
    res1=F
    res=unlist(strsplit(tmp$control[i],":",fixed=T))[1]
    if (res=="0/1" | res=="1/0") res1=T
    return(res1)
  })
  leftprobes=tmp$probe.y[gtflag]
}
leftprobes_normal=findhtprobesinnormal(data1=data1,data2=data2)
leftprobes_tumor=read.table(file=paste0("../result/tmp/",mysample,".allele.txt"),header=T,sep="\t",stringsAsFactors = F) 
leftprobes_tumor=leftprobes_tumor$SNP.Name
leftprobes_all=unique(c(leftprobes_normal,leftprobes_tumor))
res=data.frame(matrix(0,nrow=2,ncol=2))
colnames(res)=c("innormal","notinnormal")
res[1,1]=sum(leftprobes_normal %in% leftprobes_tumor)
res[2,1]=sum(! leftprobes_tumor %in% leftprobes_normal)
res[1,2]=sum(! leftprobes_normal %in% leftprobes_tumor)
res2=rbind(res2,res)
}

