#!/usr/bin/env Rscript
#SBATCH -t 1-10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org


#salloc -t 1-5 -n 61 mpirun -n 1 R --interactive
##salloc -n 40 -t 23:0:0 mpirun -n 1 Rscript --no-save --no-restore find_AUS_germlinemutaions.R
# args <- commandArgs(trailingOnly = TRUE)
# jobn=as.integer(args[1])
library(data.table)
library(DNAcopy)
library(PSCBS)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/aocs_data.RData")



#find germline vcf files from sequencing data
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
#for tumor data files
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
    names(myfile)=mysample
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

germlinevcffiles=germlinevcffiles[!is.na(germlinevcffiles)]
filecon=file("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/allsamples_aus_loh_pscbs.txt","w")
for (i in 1:length(germlinevcffiles))
{
  writeLines(names(germlinevcffiles[i]),filecon)
}
close(filecon)
#work on jobnth germlinevcffiles
PSCBS_normal_allele=function(jobn)
{
  
  if (!is.na(germlinevcffiles[jobn]))
  {
    mysample=names(germlinevcffiles[jobn])
    idx1=which(names(primarysamplefiles)==names(germlinevcffiles[jobn]))
    tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",primarysamplefiles[idx1]),sep="\t",header=T,skip=10,
                       stringsAsFactors=F)
    #germline mutation table
    normaltable=read.table(file=germlinevcffiles[jobn],skip=42,header=T,sep="\t",comment.char = "",stringsAsFactors = F)
    #to extract genotype info from controls
    tmp=normaltable$CONTROL
    genotypes=sapply(1:length(tmp),function(i){
      res2=NA
      res1=unlist(strsplit(tmp[i],":"))[1]
      if (res1=="0/1" | res1=="1/0")
      {
        res2=0.5
      }else
      {
        if (res1=="0/0")
        {
          res2=0
        }else
        {
          if (res1=="1/1")
          {
            res2=1
          }
        }
      }
      return(res2)
    })
    normaltable$genotype=genotypes
    colnames(normaltable)[1:2]=c("chr","start")
    leftprobes_normal=merge(normaltable,aocs_anno,by=c("chr","start"))
    alltable=merge(leftprobes_normal,tumortable,by.y="SNP.Name",by.x="probe")
    alltable=alltable[alltable$FILTER=="PASS",]
    alltable=alltable[!is.na(alltable$B.Allele.Freq),]
    # #for those probes not called by germline
    # uncolled_tumortable=tumortable[! tumortable$SNP.Name %in% leftprobes_normal$probe,]
    # uncolled_tumortable=uncolled_tumortable[!is.na(uncolled_tumortable$B.Allele.Freq),]
    # uncolled_tumortable=uncolled_tumortable[uncolled_tumortable$B.Allele.Freq<0.05 | uncolled_tumortable$B.Allele.Freq>0.95,]
    # betaN=sapply(1:nrow(alltable),function(i){
    #   tmp=unlist(strsplit(alltable$INFO[i],";"))
    #   tmp1=gsub("AF=","",tmp[3])
    #   tmp2=as.numeric(unlist(strsplit(tmp1,",",fixe=T)))
    #   if ((alltable$B.Allele.Freq[i]-0.5)*(tmp2[2]-0.5)>=0)
    #   {
    #     res1=tmp2[1]
    #   }else
    #   {
    #     res1=1-tmp2[1]
    #   }
    #   return(res1)
    # })
    
    PSdata=data.frame(chromosome=alltable$chr,x=alltable$start,CT=2*2^(alltable$Log.R.Ratio),betaT=alltable$B.Allele.Freq,muN=alltable$genotype)
    # #add those not called as germlines
    # PSdata1=data.frame(chromosome=uncolled_tumortable$Chr,x=uncolled_tumortable$Position,CT=2*2^(uncolled_tumortable$Log.R.Ratio),betaT=uncolled_tumortable$B.Allele.Freq,muN=rep(0,nrow(uncolled_tumortable)))
    # PSdata=rbind(PSdata,PSdata1)
    chrs = c(1:22,"X","Y")
    PSdata <- PSdata[!is.na(PSdata$x) & PSdata$chromosome %in% chrs,]
    PSdata$chromosome=gsub("X",23,PSdata$chromosome)
    PSdata$chromosome=gsub("Y",24,PSdata$chromosome)
    PSdata[,1] <- as.integer(as.character(PSdata[,1]))
    PSdata[,2] <- as.integer(as.character(PSdata[,2]))
    PSdata$betaT[!is.na(PSdata$betaT) & PSdata$betaT==-Inf] <- 0
    PSdata$betaT[!is.na(PSdata$betaT) & PSdata$betaT==+Inf] <- 1
    PSdata=PSdata[order(PSdata[,1],PSdata[,2]),]
    # PSdata$betaN[!is.na(PSdata$betaN) & PSdata$betaN==-Inf] <- 0
    # PSdata$betaN[!is.na(PSdata$betaN) & PSdata$betaN==+Inf] <- 1
    #PSdata22=PSdata[PSdata$chromosome==22,]
    gaps<- findLargeGaps(PSdata, minLength = 1e+06)
    knownSegments <- gapsToSegments(gaps)
    fit <- segmentByPairedPSCBS(PSdata, knownSegments = knownSegments,preserveScale = FALSE, seed = 48879, verbose = F,tbn=F)
    deltaAB <- estimateDeltaAB(fit, scale=1)
    fit1 <- callAB(fit, verbose = F)
    deltaLOH <- estimateDeltaLOH(fit1, midpoint=1/2)
    fit1 <- callLOH(fit1, delta=deltaLOH)
    outseg <- fit1$output
    outseg$chromosome=gsub(23,"X",outseg$chromosome)
    outseg$chromosome=gsub(24,"Y",outseg$chromosome)
    outseg <- data.frame(cbind(mysample,outseg))
    names(outseg)[1] <- "ID"
    
    output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.PSCBS.segment.txt")
    write.table(outseg,file=output,row.names=F,col.names=T,sep="\t",quote=F)
  }
  return(jobn)
}

# Sys.getenv(c("SLURM_SUBMIT_DIR"))
# Sys.getenv(c("HOST", "SLURM_JOB_ID", "SLURM_NODELIST","SLURM_NNODES", "SLURM_NTASKS",
#              "SLUR M_CPUS_PER_TASK","SLURM_CPUS_ON_NODE","SLURM_NTASKS_PER_NODE",
#              "SLURM_TASK_PID", "SLURM_ PARTITION"))
# 
opt=0
if (opt==0)
{
  library("Rmpi")
  njobs=mpi.universe.size() - 1
  mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
  .Last <- function(){
    if (is.loaded("mpi_initialize")){
      if (mpi.comm.size(1) > 0){
        print("Please use mpi.close.Rslaves() to close slaves.")
        mpi.close.Rslaves()
      }
      print("Please use mpi.quit() to quit R")
      .Call("mpi_finalize")
    }
  }
  
  mpi.bcast.Robj2slave(PSCBS_normal_allele)
  mpi.bcast.Robj2slave(germlinevcffiles)
  mpi.bcast.Robj2slave(primarysamplefiles)
  mpi.bcast.Robj2slave(aocs_anno)
  mpi.bcast.cmd(library("PSCBS"))
  mpi.bcast.cmd(library("DNAcopy"))
  res=mpi.parSapply(X=1:length(germlinevcffiles),FUN=PSCBS_normal_allele,job.num=length(germlinevcffiles))
  mpi.close.Rslaves()
  mpi.quit()
}
if (opt==1)
{
  print("start pscbs...")
  for (jobn in 1:length(germlinevcffiles))
  {
    cat(jobn,"..")
    res=PSCBS_normal_allele(jobn)
  }
}
#     




