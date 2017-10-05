#!/usr/bin/env Rscript
#for TCGA missing normal data
rm(list=ls())
library(data.table)
filelist=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/filelist.txt",header=T,sep="\t",stringsAsFactors=F)

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
allsamples=rownames(data_copynumber_filtered)
idxtumor=sapply(allsamples,function(mysample){
  res=which(grepl(paste0(mysample,"-01"),filelist$sampleid)==T & grepl("byallele",filelist$filename)==T)
  #for multiple files
  if (mysample=="TCGA-23-1023") res=res[2]
  if (mysample=="TCGA-29-2414") res=res[1]
  if (mysample=="TCGA-13-1817") res=res[2]
  if (mysample=="TCGA-13-1819") res=res[1]
  return(res)
})

idxnormal=sapply(allsamples,function(mysample){
  res=which(grepl(paste0(mysample,"-1"),filelist$sampleid)==T & grepl("byallele",filelist$filename)==T)
  return(res)  
})

idxnormal1=sapply(1:length(allsamples),function(i){
  if (length(idxnormal[[i]])==0)
  {
    res=NA
  }else
  {
    res=idxnormal[[i]]
  }
  names(res)=names(idxnormal)[i]
  return(res)
})

snp_anno=fread(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno=as.data.frame(snp_anno)
snp_anno=snp_anno[,c(1:4,9,10)]
colnames(snp_anno)[1]="Composite.Element.REF"
mutationfiles=list.files(path="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/wustl_protectedmutation/genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.1.0/",patter="genome.*.snv.*.vcf")
missingsamples=names(idxnormal1[which(is.na(idxnormal1))])
write.table(missingsamples,file="../result//missingnormalsamples.txt",row.names = F,col.names = F,quote = F)

#to find the heterozygous sites defined in normal
getgermlinemutations=function(samplename)
{
  errmsg=""
  filenames=mutationfiles[which(grepl(samplename,mutationfiles))]
  if (length(filenames)==0) errmsg="no mutation files found"
  if (errmsg=="")
  {
    filesizes=file.size(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/wustl_protectedmutation/genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.1.0/",filenames))
    filename=filenames[which.max(filesizes)]
  }
  if (errmsg=="")
  {
    if (length(filename)>0)
    {
      #the mutation calls
      mutations=fread(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/wustl_protectedmutation/genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.1.0/",filename),
                      skip=57,sep="\t",stringsAsFactors = F)
      mutations=as.data.frame(mutations)
      #pick germline mutations
      idxss=sapply(mutations$FORMAT,function(x){
        tmp=unlist(strsplit(x,":",fixed=T))
        res1=which(tmp=="SS")
      })
      sum(is.na(idxss))
      ##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">
      SS=sapply(1:nrow(mutations),function(i){
        res1=NA
        if (!is.na(idxss[i]))
        {
          tmp=unlist(strsplit(mutations[i,10],":",fixed=T))
          res1=tmp[idxss[i]]
        }
        res1
      })
      idxgermline=which(SS=="1")
      res1=data.frame(chr=mutations$`#CHROM`[idxgermline],pos=mutations$POS[idxgermline],ref=mutations$REF[idxgermline],alt=mutations$ALT[idxgermline])
      # write.table(res1,file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",samplename,".germline.vcf"),
      #             row.names = F,sep="\t",quote=F)
      idx1=which(names(idxtumor)==samplename)
      tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[idx1]],"/",filelist$filename[idxtumor[idx1]]),
                            skip = 1, sep="\t",header=T,stringsAsFactors = F)
      tumortable1=merge(tumortable,snp_anno)
      tmp=merge(res1,tumortable1,by.x=c("chr","pos"),by.y=c("Chromosome","Physical Position"))
      tmp=tmp[,-c(3,4,8)]
      tmp=tmp[,c(3,4,5,1,2)]
      write.table(tmp,file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",samplename,".heterog.txt"),
             row.names = F,sep="\t",quote=F)
    }
  }else
  {
    print(paste0(samplename,": ",errmsg))
  }
}


for (i in 4:length(missingsamples))
{
  cat(i,'..')
  getgermlinemutations(samplename = missingsamples[i])
}

library(DNAcopy)
mpi.bcast.cmd(library(DNAcopy))
CBS_allele_dh=function(samplename)
{
  
  mysample=samplename
  filenames=mutationfiles[which(grepl(samplename,mutationfiles))]
  if (length(filenames)>0)
  {
    idx1=which(names(idxtumor)==samplename)
    tumortable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[idx1]],"/",filelist$filename[idxtumor[idx1]]),
                          skip = 1, sep="\t",header=T,stringsAsFactors = F)
    colnames(tumortable)[1]="Probe.Set.ID"
    tumortable$Signal_A[which(tumortable$Signal_A<0)]=0
    tumortable$Signal_B[which(tumortable$Signal_B<0)]=0
    #pick the heterozygous sites defined in normal
    normaltable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".heterog.txt"),header=T,sep="\t",stringsAsFactors = F)
    colnames(normaltable)[1]="Probe.Set.ID"
    tumortable=merge(tumortable,normaltable,by.x="Probe.Set.ID",by.y="Probe.Set.ID")
    tumortable=tumortable[,c(1:3,6:7)] 
    colnames(tumortable)[c(2:5)]=c("Signal_A","Signal_B","Chromosome","Physical.Position")
    dh=abs(tumortable$Signal_A/(tumortable$Signal_A+tumortable$Signal_B)-0.5)
    
    tumortable=cbind.data.frame(tumortable,dh=dh)
    tumortable$Physical.Position=as.integer(tumortable$Physical.Position)
    
    chrs = paste0('chr',c(1:22,'X','Y'))
    CNA.obj = CNA(tumortable$dh, 
                  tumortable$Chromosome, 
                  tumortable$Physical.Position, 
                  data.type = "logratio",
                  sampleid=mysample)
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
        segs <- segment(CNAchr.smoothed,alpha=0.1)
        segsall=rbind(segsall,segs$output)
      }
    }
    #write segment results to tmp folder
    output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.cbs.segment.txt")
    write.table(segsall,file=output,row.names=F,col.names=T,sep="\t",quote=F)
    dh_probe_file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.txt")
    write.table(tumortable,file=dh_probe_file,row.names=F,col.names=T,sep="\t",quote=F)
  }
}

for (mysample in missingsamples)
{
  print(mysample)
  CBS_allele_dh(mysample)
}