#!/usr/bin/env Rscript

#form input data-------
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
sum(is.na(idxnormal1))
#[1] 12 samples without normal allele files
idxnormal=idxnormal1



snp_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno <- snp_anno[,1:4]

ascat_form_input=function()
{
  if (exists("tumorbaf")) rm(tumorbaf)
  if (exists("normalbaf")) rm(normalbaf)
  if (exists("tumorlogr")) rm(tumorlogr)
  if (exists("normallogr")) rm(normallogr)
  
  #first only consider samples having normals
  idxs=which(!is.na(idxnormal))
  #pick a few samples
  #idxs=idxs[idxs<10]
  for (i in idxs)
  {
    mysample=allsamples[i]
    tumortable=fread(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[i]],"/",filelist$filename[idxtumor[i]]),
                          skip = 1, sep="\t",header=T,stringsAsFactors = F)
    
    tumortable$Signal_A[which(!is.na(tumortable$Signal_A) & tumortable$Signal_A<0)] <- NA
    tumortable$Signal_B[which(!is.na(tumortable$Signal_A) & tumortable$Signal_B<0)] <- NA
    
    ##tumortable$Signal_A[which(tumortable$Signal_A<0)]=0
    ##tumortable$Signal_B[which(tumortable$Signal_B<0)]=0
    
    colnames(tumortable)[1]="Probe.Set.ID"
    
    
    #if it has normal allele files, find heterogeneity probes first
    if (!is.na(idxnormal[i]))
    {
      normaltable=fread(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxnormal[i]],"/",filelist$filename[idxnormal[i]]),
                             skip = 1, sep="\t",header=T,stringsAsFactors = F)
      
      normaltable$Signal_A[which(!is.na(normaltable$Signal_A) & normaltable$Signal_A<0)] <- NA
      normaltable$Signal_B[which(!is.na(normaltable$Signal_A) & normaltable$Signal_B<0)] <- NA
      
      colnames(normaltable)[1]="Probe.Set.ID"
      
      tumortable=merge(tumortable,snp_anno,by="Probe.Set.ID")
      normaltable=merge(normaltable,snp_anno,by="Probe.Set.ID")
      
      tumortable=tumortable[order(tumortable$Chromosome,tumortable$Physical.Position),]
      normaltable=normaltable[order(normaltable$Chromosome,normaltable$Physical.Position),]
      
      
      tumortable$Physical.Position=as.integer(tumortable$Physical.Position)
      normaltable$Physical.Position=as.integer(normaltable$Physical.Position)
      
      tumortable$LogR=log((tumortable$Signal_A+tumortable$Signal_B)/2)
      normaltable$LogR=log((normaltable$Signal_A+normaltable$Signal_B)/2)
      tumortable$BAF=tumortable$Signal_B/(tumortable$Signal_A+tumortable$Signal_B)
      normaltable$BAF=normaltable$Signal_B/(normaltable$Signal_A+normaltable$Signal_B)
      
      if (!exists("tumorbaf"))
      {
        tumorbaf=data.frame(matrix(NA,nrow=nrow(tumortable),ncol=3))
        colnames(tumorbaf)=c("chrs","pos",mysample)
        rownames(tumorbaf)=tumortable$Probe.Set.ID
        tumorbaf$chrs=tumortable$Chromosome
        tumorbaf$chrs=gsub("---",NA,tumorbaf$chrs)
        tumorbaf$pos=tumortable$Physical.Position
        tumorbaf[,ncol(tumorbaf)]=tumortable$BAF
      }else
      {
        tumorbaf=cbind(tumorbaf,tumortable$BAF)
        colnames(tumorbaf)[ncol(tumorbaf)]=mysample
      }
      
      if (!exists("normalbaf"))
      {
        normalbaf=data.frame(matrix(NA,nrow=nrow(normaltable),ncol=3))
        colnames(normalbaf)=c("chrs","pos",mysample)
        rownames(normalbaf)=normaltable$Probe.Set.ID
        normalbaf$chrs=normaltable$Chromosome
        normalbaf$chrs=gsub("---",NA,normalbaf$chrs)
        normalbaf$pos=normaltable$Physical.Position
        normalbaf[,ncol(normalbaf)]=normaltable$BAF
      }else
      {
        normalbaf=cbind(normalbaf,normaltable$BAF)
        colnames(normalbaf)[ncol(normalbaf)]=mysample
      }
      
      if (!exists("tumorlogr"))
      {
        tumorlogr=data.frame(matrix(NA,nrow=nrow(tumortable),ncol=3))
        colnames(tumorlogr)=c("chrs","pos",mysample)
        rownames(tumorlogr)=tumortable$Probe.Set.ID
        tumorlogr$chrs=tumortable$Chromosome
        tumorlogr$chrs=gsub("---",NA,tumorlogr$chrs)
        tumorlogr$pos=tumortable$Physical.Position
        tumorlogr[,ncol(tumorlogr)]=tumortable$LogR
      }else
      {
        tumorlogr=cbind(tumorlogr,tumortable$LogR)
        colnames(tumorlogr)[ncol(tumorlogr)]=mysample
      }
      
      if (!exists("normallogr"))
      {
        normallogr=data.frame(matrix(NA,nrow=nrow(normaltable),ncol=3))
        colnames(normallogr)=c("chrs","pos",mysample)
        rownames(normallogr)=normaltable$Probe.Set.ID
        normallogr$chrs=normaltable$Chromosome
        normallogr$chrs=gsub("---",NA,normallogr$chrs)
        normallogr$pos=normaltable$Physical.Position
        normallogr[,ncol(normallogr)]=normaltable$LogR
      }else
      {
        normallogr=cbind(normallogr,normaltable$LogR)
        colnames(normallogr)[ncol(normallogr)]=mysample
      }
    }
    if (i %%10==0) cat(i,"..")
  }
  tumorbaf=tumorbaf[tumorbaf$chrs!="Y",]
  normalbaf=normalbaf[normalbaf$chrs!="Y",]
  tumorlogr=tumorlogr[tumorlogr$chrs!="Y",]
  normallogr=normallogr[normallogr$chrs!="Y",]
  
  # tumorbaf1=complete.cases(tumorbaf)
  # normalbaf1=complete.cases(normalbaf)
  # tumorlogr1=complete.cases(tumorlogr)
  # normallogr1=complete.cases(normallogr)
  idx=is.na(tumorbaf$pos)
  tumorbaf=tumorbaf[!idx,]
  normalbaf=normalbaf[!idx,]
  tumorlogr=tumorlogr[!idx,]
  normallogr=normallogr[!idx,]
  write.table(tumorbaf,file="../result/tcga_tumorbaf.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  write.table(normalbaf,file="../result/tcga_normalbaf.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  write.table(tumorlogr,file="../result/tcga_tumorlogr.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  write.table(normallogr,file="../result/tcga_normallogr.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  save(tumorbaf,tumorlogr,normalbaf,normallogr,file="../result/tcga_ascat.RData")
  filecon=file("../result/ascat_samplelist.txt","w")
  for (i in 1:(ncol(tumorbaf)-2))
  {
    writeLines(colnames(tumorbaf)[i+2],filecon)
    Tumor_LogR=tumorlogr[,c(1:2,2+i)]
    Tumor_BAF = tumorbaf[,c(1:2,2+i)]
    Germline_LogR = normallogr[,c(1:2,2+i)]
    Germline_BAF = normalbaf[,c(1:2,2+i)]
    write.table(Tumor_LogR,file=paste0("../result/ascat_input/",colnames(tumorbaf)[i+2],"_tumorlogr.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
    write.table(Tumor_BAF,file=paste0("../result/ascat_input/",colnames(tumorbaf)[i+2],"_tumorbaf.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
    write.table(Germline_LogR,file=paste0("../result/ascat_input/",colnames(tumorbaf)[i+2],"_normallogr.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
    write.table(Germline_BAF,file=paste0("../result/ascat_input/",colnames(tumorbaf)[i+2],"_normalbaf.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
    if (i %% 100==0) cat(i,"..")
  }
  close(filecon)
  # write.table(tumorbaf[,1:7],file="../result/tcga1_tumorbaf.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  # write.table(normalbaf[,1:7],file="../result/tcga1_normalbaf.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  # write.table(tumorlogr[,1:7],file="../result/tcga1_tumorlogr.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  # write.table(normallogr[,1:7],file="../result/tcga1_normallogr.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  # write.table(tumorbaf[,c(1:2,5)],file="../result/tcga2_tumorbaf.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  # write.table(normalbaf[,c(1:2,5)],file="../result/tcga2_normalbaf.txt",row.names = T,col.names = NA,sep="\t",quote=F)
  
  #tumor only samples
  
  idxs1=which(is.na(idxnormal))
  filecon=file("../result/ascat_tumoronlysamplelist.txt","w")
  for (i in idxs1)
  {
    mysample=allsamples[i]
    writeLines(mysample,filecon)
    tumortable=fread(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/",filelist$fileid[idxtumor[i]],"/",filelist$filename[idxtumor[i]]),
                     skip = 1, sep="\t",header=T,stringsAsFactors = F)
    
    tumortable$Signal_A[which(!is.na(tumortable$Signal_A) & tumortable$Signal_A<0)] <- NA
    tumortable$Signal_B[which(!is.na(tumortable$Signal_A) & tumortable$Signal_B<0)] <- NA
    
    ##tumortable$Signal_A[which(tumortable$Signal_A<0)]=0
    ##tumortable$Signal_B[which(tumortable$Signal_B<0)]=0
    
    colnames(tumortable)[1]="Probe.Set.ID"
    tumortable=merge(tumortable,snp_anno,by="Probe.Set.ID")
    tumortable=tumortable[order(tumortable$Chromosome,tumortable$Physical.Position),]
    tumortable$Physical.Position=as.integer(tumortable$Physical.Position)
    tumortable$LogR=log((tumortable$Signal_A+tumortable$Signal_B)/2)
    tumortable$BAF=tumortable$Signal_B/(tumortable$Signal_A+tumortable$Signal_B)
    rownames(tumortable)=tumortable$Probe.Set.ID
    Tumor_LogR=tumortable[,c("Chromosome","Physical.Position","LogR")]
    Tumor_BAF =tumortable[,c("Chromosome","Physical.Position","BAF")]
    rownames(Tumor_LogR)=rownames(Tumor_BAF)=tumortable$Probe.Set.ID
    colnames(Tumor_LogR)=colnames(Tumor_BAF)=c("chrs","pos",mysample)
    Tumor_LogR$chrs=gsub("---",NA,Tumor_LogR$chrs)
    Tumor_BAF$chrs=gsub("---",NA,Tumor_BAF$chrs)
    idx=which(Tumor_LogR$chrs!="Y" | !is.na(Tumor_BAF$pos))
    
    Tumor_LogR1=Tumor_LogR[idx,]
    Tumor_BAF1=Tumor_BAF[idx,]
    rownames(Tumor_BAF1)=rownames(Tumor_LogR1)=rownames(Tumor_BAF)[idx]

    write.table(Tumor_LogR1,file=paste0("../result/ascat_input/",mysample,"_tumorlogr.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
    write.table(Tumor_BAF1,file=paste0("../result/ascat_input/",mysample,"_tumorbaf.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
  }
  close(filecon)
}

#Ran ascat
library(ASCAT)
Sys.time()
ascat.bc = ascat.loadData("../result/tcga1_tumorlogr.txt","../result/tcga1_tumorbaf.txt","../result/tcga1_normallogr.txt", "../result/tcga1_normalbaf.txt")
Sys.time()
ascat.bc = ascat.GCcorrect(ascat.bc,"/fh/fast/dai_j/CancerGenomics/Tools/database/other/GC_AffySNP6_102015.txt")
Sys.time()
ascat.plotRawData(ascat.bc)
Sys.time()
ascat.bc = ascat.aspcf(ascat.bc)
Sys.time()
ascat.plotSegmentedData(ascat.bc)
Sys.time()
ascat.output = ascat.runAscat(ascat.bc)
Sys.time()

add_pos_nA=function(nA=ascat.output$nA,snppos=ascat.bc$SNPpos)
{
  idx=match(rownames(nA),rownames(snppos))
  snppos=snppos[idx,]
  res=cbind(snppos,nA)
}
nA=add_pos_nA()
nB=add_pos_nA(nA=ascat.output$nB)

#without normals
Sys.time()
ascat1.bc = ascat.loadData("../result/tcga1_tumorlogr.txt","../result/tcga1_tumorbaf.txt")
Sys.time()
ascat1.bc = ascat.GCcorrect(ascat1.bc,"/fh/fast/dai_j/CancerGenomics/Tools/database/other/GC_AffySNP6_102015.txt")
ascat1.gg = ascat.predictGermlineGenotypes(ascat1.bc, "AffySNP6") 
Sys.time()
#ascat.plotRawData(ascat1.bc)
ascat1.bc = ascat.aspcf(ascat1.bc,ascat.gg=ascat1.gg) 
Sys.time()
ascat.plotSegmentedData(ascat1.bc)
Sys.time()
ascat1.output = ascat.runAscat(ascat1.bc)
Sys.time()

#check a failed one
ascat2.bc = ascat.loadData("../result/tcga2_tumorlogr.txt","../result/tcga2_tumorbaf.txt")
ascat2.bc = ascat.GCcorrect(ascat2.bc,"/fh/fast/dai_j/CancerGenomics/Tools/database/other/GC_AffySNP6_102015.txt")
ascat2.gg = ascat.predictGermlineGenotypes(ascat2.bc, "AffySNP6") 
Sys.time()
#ascat.plotRawData(ascat1.bc)
ascat2.bc = ascat.aspcf(ascat2.bc,ascat.gg=ascat2.gg) 
Sys.time()
ascat.plotSegmentedData(ascat2.bc)
Sys.time()
ascat2.output = ascat.runAscat(ascat2.bc)
Sys.time()


#generate inputs for AUS samples
ascat_form_input_aus=function()
{
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
  filecon=file("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_AUS_samplelist.txt","w")
  for (i in 1:length(primarysamplefiles))
  {
    writeLines(names(primarysamplefiles)[i],filecon)
  }
  close(filecon)
  #omni_1.1=fread(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni25M-8v1-1_B.annotated.txt",header=T,sep="\t",fill=T,stringsAsFactors=F)
  load("/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni25M-8v1-1_B.annotated.RData")
  #omni_1.1=as.data.frame(omni_1.1)
  load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/probes2keep.RData")
  df_probes2keep=data.frame(Name=probes2keep)
  aocs_anno=merge(df_probes2keep,omni_1.1)
  aocs_anno=aocs_anno[,1:3]
  colnames(aocs_anno)=c("probe","chr","start")
  aocs_anno$probe=as.character(aocs_anno$probe)
  
  chrs=c(1:22,"X")
  for (i in 1:length(primarysamplefiles))
  {
    if (i %% 10==0) cat(i,"..")
    mysample=names(primarysamplefiles)[i]
    tumortable=fread(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",primarysamplefiles[i]),sep="\t",header=T,skip=10,
                          stringsAsFactors=F)
    tumortable=as.data.frame(tumortable)
    tumortable1=tumortable[,c(1,31,32)]
    colnames(tumortable1)=c("probe","baf","logr")
    tumortable1=merge(aocs_anno,tumortable1)
    Tumor_LogR=tumortable1[,c("chr","start","logr")]
    Tumor_BAF =tumortable1[,c("chr","start","baf")]
    rownames(Tumor_LogR)=rownames(Tumor_BAF)=tumortable1$probe
    colnames(Tumor_LogR)=colnames(Tumor_BAF)=c("chrs","pos",mysample)
    idx=which(Tumor_LogR$chrs %in% chrs & !is.na(Tumor_BAF$pos))
    
    Tumor_LogR1=Tumor_LogR[idx,]
    Tumor_BAF1=Tumor_BAF[idx,]
    rownames(Tumor_BAF1)=rownames(Tumor_LogR1)=rownames(Tumor_BAF)[idx]
    
    write.table(Tumor_LogR1,file=paste0("../result/ascat_input/",mysample,"_tumorlogr.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
    write.table(Tumor_BAF1,file=paste0("../result/ascat_input/",mysample,"_tumorbaf.txt"),row.names = T,col.names = NA,sep="\t",quote=F)
  }
}
  


