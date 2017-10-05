#!/usr/bin/env Rscript
#http://firebrowse.org/?cohort=OV#

transformcolnames=function(data=median_mrna)
{
  tmp=data
  tmp=as.matrix(tmp)
  tmp=as.numeric(tmp)
  tmp1=matrix(tmp,ncol=ncol(data),byrow=FALSE)
  colnames(tmp1)=colnames(data)
  rownames(tmp1)=rownames(data)
  data=tmp1
  tmp=sapply(1:ncol(data),function(x){
    tmp1=unlist(strsplit(colnames(data)[x],".",fixed=TRUE))
    return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
  })
  length(unique(tmp))
  colnames(data)=tmp
  #check tumor part name
  tmp=sapply(1:ncol(data),function(x){
    tmp1=unlist(strsplit(colnames(data)[x],"-",fixed=TRUE))
    return(substr(tmp1[4],1,2))
  })
  print(table(tmp))
  return(data)
}

extracttumordata=function(data,offset=0)
{
  tmp=sapply((offset+1):ncol(data),function(x){
    tmp1=unlist(strsplit(colnames(data)[x],"-",fixed=TRUE))
    return(substr(tmp1[4],1,2))
  })
  table(tmp)
  idx_tumor=which(tmp %in% c("01"))
  data_all=data
  data=data_all[,offset+idx_tumor]
  return(data)
}

remove_tumorid=function(data=mrna_exon)
{
  if (grepl("TCGA.",colnames(data)[1],fixed=TRUE)==TRUE)
  {
    tmp=sapply(1:ncol(data),function(x){
      tmp1=unlist(strsplit(colnames(data)[x],".",fixed=TRUE))
      return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
    })
    length(unique(tmp))
    colnames(data)=tmp
  }
  
  tmp=unlist(strsplit(colnames(data)[1],"-",fixed=TRUE))
  #if include tumor id 01,02,10,11...
  if (length(tmp>=4))
  {
    tmp1=sapply(colnames(data),function(x){
      tmp=unlist(strsplit(x,"-",fixed=TRUE))
      res=paste0(tmp[1:3],collapse="-")
    })
  }
  if (ncol(data)>length(unique(tmp1)))
  {
    warning("There are multiple samples for a same patient...")
  }
  colnames(data)=tmp1  
  return(data)
}

load("/fh/fast/dai_j/CancerGenomics/Ovarian/mrna_copynumber_methylation_mutation.RData")
#to get mrna---------------------
median_mrna=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.mRNA_Preprocess_Median.Level_3.2016012800.0.0/OV.medianexp.txt",
                       sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
median_mrna=median_mrna[2:nrow(median_mrna),]
colnames(copynumber)[1]
rownames(median_mrna)=as.character(median_mrna[,1])
median_mrna=median_mrna[,2:ncol(median_mrna)]
dim(median_mrna)
#[1] 18632   594
head(colnames(median_mrna))
#[1] "TCGA.01.0628.11A.01R.0361.03" "TCGA.01.0630.11A.01R.0361.03" "TCGA.01.0631.11A.01R.0361.03"

tmp=median_mrna
tmp=as.matrix(tmp)
tmp=as.numeric(tmp)
tmp1=matrix(tmp,ncol=ncol(median_mrna),byrow=FALSE)
colnames(tmp1)=colnames(median_mrna)
rownames(tmp1)=rownames(median_mrna)
median_mrna=tmp1
tmp=sapply(1:ncol(median_mrna),function(x){
  tmp1=unlist(strsplit(colnames(median_mrna)[x],".",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
})
length(unique(tmp))
colnames(median_mrna)=tmp
#check tumor part name
tmp=sapply(1:ncol(median_mrna),function(x){
  tmp1=unlist(strsplit(colnames(median_mrna)[x],"-",fixed=TRUE))
  return(substr(tmp1[4],1,2))
})
table(tmp)
#01  02  11 
#569  17   8
#01: primary solid tumor; 02: recurrent solid tumor
#idx_tumor=which(tmp %in% c("01","02"))
#TCGA analysis only used 01
idx_tumor=which(tmp %in% c("01"))
median_mrna_all=median_mrna
median_mrna=median_mrna_all[,idx_tumor]
dim(median_mrna)
#[1] 18632   569
# #check patient part name, useful if 01 02 both are included
# tmp=sapply(1:ncol(median_mrna),function(x){
#   tmp1=unlist(strsplit(colnames(median_mrna)[x],"-",fixed=TRUE))
#   return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
# })
# length(tmp)
# #[1] 586
# length(unique(tmp))
# #[1] 571
# # 586 tumor samples, 571 unique tumor samples
#---------------------------

# #check other mrna data:
mrna_u133=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.Level_3.2016012800.0.0/OV.transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.data.txt",
                     sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
mrna_u133=mrna_u133[2:nrow(mrna_u133),]
rownames(mrna_u133)=as.character(mrna_u133[,1])
mrna_u133=mrna_u133[,2:ncol(mrna_u133)]
dim(mrna_u133)
# [1] 12042 543
mrna_u133=transformcolnames(mrna_u133)
mrna_u133=extracttumordata(data=mrna_u133)
mrna_u133=remove_tumorid(mrna_u133)
dim(mrna_u133)
#[1] 12042   520

mrna_4502a=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3.2016012800.0.0/OV.transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt",
                      sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
mrna_4502a=mrna_4502a[2:nrow(mrna_4502a),]
rownames(mrna_4502a)=as.character(mrna_4502a[,1])
mrna_4502a=mrna_4502a[,2:ncol(mrna_4502a)]
dim(mrna_4502a)
# [1] 17814   561
mrna_4502a=transformcolnames(mrna_4502a)
mrna_4502a=extracttumordata(data=mrna_4502a)
mrna_4502a=remove_tumorid(mrna_4502a)
dim(mrna_4502a)
#[1] 17814   541

mrna_4502a_72=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3.2016012800.0.0/OV.transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt",
                         sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
mrna_4502a_72=mrna_4502a_72[2:nrow(mrna_4502a_72),]
rownames(mrna_4502a_72)=as.character(mrna_4502a_72[,1])
mrna_4502a_72=mrna_4502a_72[,2:ncol(mrna_4502a_72)]
dim(mrna_4502a_72)
#[1] 17814    36
mrna_4502a_72=transformcolnames(mrna_4502a_72)
mrna_4502a_72=extracttumordata(data=mrna_4502a_72)
mrna_4502a_72=remove_tumorid(mrna_4502a_72)
sum(colnames(mrna_4502a_72) %in% colnames(mrna_4502a))
#[1] 0 The two datasets don't overlap
mrna_4502a=cbind(mrna_4502a,mrna_4502a_72)
dim(mrna_4502a)
#[1] 17814   572

# 
# mrna_exon=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_exon__huex_1_0_st_v2__lbl_gov__Level_3__quantile_normalization_gene__data.Level_3.2016012800.0.0/OV.exon__huex_1_0_st_v2__lbl_gov__Level_3__quantile_normalization_gene__data.data.txt",
#                        sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
# #it shows that mrna_median is mrna_exon

#methylation
# mean_methylation=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Methylation_Preprocess.Level_3.2016012800.0.0/OV.meth.by_mean.data.txt",
#                        sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
# mean_methylation=mean_methylation[2:nrow(mean_methylation),]
# rownames(mean_methylation)=as.character(mean_methylation[,1])
# mean_methylation=mean_methylation[,2:ncol(mean_methylation)]
# dim(mean_methylation)
# #[1] 13772   612
# head(colnames(mean_methylation))
# #[1] "TCGA.01.0628.11" "TCGA.01.0630.11" "TCGA.01.0631.11" "TCGA.01.0633.11" "TCGA.01.0636.11" "TCGA.01.0637.11"
# 
# #make name of patient in mrna consistent with that in methylation
# sample_median_mrna=sapply(colnames(median_mrna),function(x){
#   strings=unlist(strsplit(x,".",fixed=TRUE))
#   tumorid=substr(strings[4],1,2)
#   res=paste0(strings[1],".",strings[2],".",strings[3],".",tumorid)
# })
# 
# #find common set of genes in mrna and methylation
# colnames(median_mrna)=sample_median_mrna
# tmp=rownames(median_mrna) %in% rownames(mean_methylation)
# medianmrna_meanmethylation=median_mrna[tmp,]
# tmp=rownames(mean_methylation) %in% rownames(median_mrna)
# meanmethylation_medianmrna=mean_methylation[tmp,]
# tmp=sapply(1:nrow(medianmrna_meanmethylation),function(x){
#   res=which(rownames(meanmethylation_medianmrna)==rownames(medianmrna_meanmethylation)[x])
# })
# meanmethylation_medianmrna=meanmethylation_medianmrna[tmp,]
# 
# 
# #find common set of samples
# tmp=colnames(median_mrna) %in% colnames(mean_methylation)
# medianmrna_meanmethylation=medianmrna_meanmethylation[,tmp]
# tmp=colnames(mean_methylation) %in% colnames(median_mrna)
# meanmethylation_medianmrna=meanmethylation_medianmrna[,tmp]
# tmp=sapply(1:ncol(medianmrna_meanmethylation),function(x){
#   res=which(colnames(meanmethylation_medianmrna)==colnames(medianmrna_meanmethylation)[x])
# })
# meanmethylation_medianmrna=meanmethylation_medianmrna[,tmp]
# 
# #compute correlation:
# cor_meanmethylation_medianmrna=sapply(1:nrow(meanmethylation_medianmrna),function(x){
#   res=cor(as.numeric(meanmethylation_medianmrna[x,]),as.numeric(medianmrna_meanmethylation[x,]),use="complete")
# })
# cor_meanmethylation_medianmrna_order=cor_meanmethylation_medianmrna[order(cor_meanmethylation_medianmrna,decreasing=TRUE)]
# names(cor_meanmethylation_medianmrna_order)=rownames(meanmethylation_medianmrna)[order(cor_meanmethylation_medianmrna,decreasing=TRUE)]

#to get methylation---------
mincor_methylation=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Methylation_Preprocess.Level_3.2016012800.0.0/OV.meth.by_min_expr_corr.data.txt",
                              sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
mincor_methylation=mincor_methylation[2:nrow(mincor_methylation),]
rownames(mincor_methylation)=as.character(mincor_methylation[,1])
mincor_methylation=mincor_methylation[,5:ncol(mincor_methylation)]
tmp=mincor_methylation
tmp=as.matrix(tmp)
tmp=as.numeric(tmp)
tmp1=matrix(tmp,ncol=ncol(mincor_methylation),byrow=FALSE)
colnames(tmp1)=colnames(mincor_methylation)
rownames(tmp1)=rownames(mincor_methylation)
mincor_methylation=tmp1

tmp=sapply(1:ncol(mincor_methylation),function(x){
  tmp1=unlist(strsplit(colnames(mincor_methylation)[x],".",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
})
length(unique(tmp))
colnames(mincor_methylation)=tmp
#check tumor part name
tmp=sapply(1:ncol(mincor_methylation),function(x){
  tmp1=unlist(strsplit(colnames(mincor_methylation)[x],"-",fixed=TRUE))
  return(substr(tmp1[4],1,2))
})
table(tmp)
# 01  02  11 
#582  18  12 
#idx_tumor=which(tmp %in% c("01","02"))
idx_tumor=which(tmp %in% c("01"))
mincor_methylation_all=mincor_methylation
mincor_methylation=mincor_methylation_all[,idx_tumor]
dim(mincor_methylation)
#[1] 12087   582
# #check patient part name
# tmp=sapply(1:ncol(mincor_methylation),function(x){
#   tmp1=unlist(strsplit(colnames(mincor_methylation)[x],"-",fixed=TRUE))
#   return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
# })
# length(tmp)
# #[1] 600
# length(unique(tmp))
# #[1] 584
# # 600 tumor samples, 584 unique tumor samples.
#------------------------

# #find common set of genes in mrna and methylation
# #colnames(median_mrna)=sample_median_mrna
# tmp=rownames(median_mrna) %in% rownames(mincor_methylation)
# medianmrna_mincormethylation=median_mrna[tmp,]
# tmp=rownames(mincor_methylation) %in% rownames(median_mrna)
# mincormethylation_medianmrna=mincor_methylation[tmp,]
# tmp=sapply(1:nrow(medianmrna_mincormethylation),function(x){
#   res=which(rownames(mincormethylation_medianmrna)==rownames(medianmrna_mincormethylation)[x])
# })
# mincormethylation_medianmrna=mincormethylation_medianmrna[tmp,]
# #find common set of samples
# tmp=colnames(median_mrna) %in% colnames(mincor_methylation)
# medianmrna_mincormethylation=medianmrna_mincormethylation[,tmp]
# tmp=colnames(mincor_methylation) %in% colnames(median_mrna)
# mincormethylation_medianmrna=mincormethylation_medianmrna[,tmp]
# tmp=sapply(1:ncol(medianmrna_mincormethylation),function(x){
#   res=which(colnames(mincormethylation_medianmrna)==colnames(medianmrna_mincormethylation)[x])
# })
# mincormethylation_medianmrna=mincormethylation_medianmrna[,tmp]
# #compute correlation:
# cor_mincormethylation_medianmrna=sapply(1:nrow(mincormethylation_medianmrna),function(x){
#   res=cor(as.numeric(mincormethylation_medianmrna[x,]),as.numeric(medianmrna_mincormethylation[x,]),use="complete")
# })
# cor_mincormethylation_medianmrna_order=cor_mincormethylation_medianmrna[order(cor_mincormethylation_medianmrna,decreasing=TRUE)]
# names(cor_mincormethylation_medianmrna_order)=rownames(mincormethylation_medianmrna)[order(cor_mincormethylation_medianmrna,decreasing=TRUE)]
# tail(cor_mincormethylation_medianmrna_order)
# 
# methylation_mincor_common=mincormethylation_medianmrna
# mrna_median_common=medianmrna_mincormethylation
# save(median_mrna,mrna_median_common,mincor_methylation,methylation_mincor_common,file="../data/firehose.RData")



#mRNAseq
mrnaseq=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.mRNAseq_Preprocess.Level_3.2016012800.0.0/OV.uncv2.mRNAseq_RSEM_normalized_log2.txt",
                   sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
rownames(mrnaseq)=as.character(mrnaseq[,1])
mrnaseq=mrnaseq[,2:ncol(mrnaseq)]
tmp=sapply(1:ncol(mrnaseq),function(x){
  tmp1=unlist(strsplit(colnames(mrnaseq)[x],".",fixed=TRUE))
  return(substr(tmp1[4],1,2))
})
table(tmp)
idx_tumor=which(tmp %in% c("01"))
mrnaseq=mrnaseq[,idx_tumor]
mrnaseq=remove_tumorid(mrnaseq)


#copy number
library(CNTools)
#to get copy number-----------------------
#SNP6
snp6copynumberhg18=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg.seg.txt",
                              sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
colnames(snp6copynumberhg18)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
dim(snp6copynumberhg18)
length(unique(snp6copynumberhg18$Sample))
cnseg_snp6hg18=CNSeg(snp6copynumberhg18)
# geneInfo was built for human genes based on build 36.
data("geneInfo")
snp6hg18_median <- getRS(cnseg_snp6hg18, by = "gene", imput = FALSE, XY = TRUE, geneMap = geneInfo, what = "median")
snp6hg18_median=rs(snp6hg18_median)
tmp=sapply(6:ncol(snp6hg18_median),function(x){
  tmp1=unlist(strsplit(colnames(snp6hg18_median)[x],"-",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
})
length(unique(tmp))
colnames(snp6hg18_median)[6:ncol(snp6hg18_median)]=tmp
#check tumor part name
tmp=sapply(6:ncol(snp6hg18_median),function(x){
  tmp1=unlist(strsplit(colnames(snp6hg18_median)[x],"-",fixed=TRUE))
  return(substr(tmp1[4],1,2))
})
table(tmp)
#01  02  10  11 
#583  18 441 130 
idx_tumor=which(tmp %in% c("01","02"))
snp6hg18_median_all=snp6hg18_median
snp6hg18_median=snp6hg18_median_all[,5+idx_tumor]
#there a few duplicated gene names
idx_dup=which(snp6hg18_median_all$genename %in% c("CRSP8","DUB3","PRG2","SKIP"))
cbind(snp6hg18_median_all[idx_dup,1:5],idx_dup)
snp6hg18_median_all$genename[21250]="CRSP8P"

snp6hg18_median_all[19724,1:5] #Gene ID: 402164, discontinued on 4-Apr-2009,remove it
snp6hg18_median_all$genename[19724]=NA
snp6hg18_median_all[12823,1:5]
snp6hg18_median_all$genename[12823]="PLPPR3"
snp6hg18_median_all[16110,1:5]
snp6hg18_median_all$genename[16110]="SPHKAP"
snp6hg18_median_all=snp6hg18_median_all[!is.na(snp6hg18_median_all$genename),]
snp6hg18_median=snp6hg18_median_all[,5+idx_tumor]
rownames(snp6hg18_median)=snp6hg18_median_all$genename
#check patient part name
tmp=sapply(1:ncol(snp6hg18_median),function(x){
  tmp1=unlist(strsplit(colnames(snp6hg18_median)[x],"-",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
})
length(tmp)
#[1] 601
length(unique(tmp))
#[1] 586
# 601 tumor samples, 586 unique tumor samples
#----------------------------------------
# tmp=sapply(6:ncol(snp6hg18_median_all),function(x){
#   tmp1=unlist(strsplit(colnames(snp6hg18_median_all)[x],"-",fixed=TRUE))
#   return(substr(tmp1[4],1,2))
# })
# idx_normal=which(tmp %in% c("10","11"))
# snp6hg18_median_normal=snp6hg18_median_all[,5+idx_normal]
# 
# 
# tmp1=sapply(6:ncol(snp6hg18_median_all),function(x){
#   tmp1=unlist(strsplit(colnames(snp6hg18_median_all)[x],"-",fixed=TRUE))
#   return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
# })
# length(unique(tmp1))

#for TCGA segments, include tumor and normal segments
readcnasegs=function(hg19=T,segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt")
{
  library(CNTools)
  segments=read.table(file=segfile,sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
  colnames(segments)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  dim(segments)
  length(unique(segments$ID))
  segments=segments[!is.na(segments$seg.mean),]
  cnseg_snp6=CNSeg(segments)
  #data("geneInfo"): it includes duplicate genes, need to be removed
  #load 
  load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")
  hg19_firehose=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
  idxdup=duplicated(hg19_firehose$gene)
  hg19_firehose=hg19_firehose[! idxdup,]
  #change the column name of chr
  colnames(hg19_firehose)[2]="chrom"
  if (hg19==T)
  {
    nstart=5
    snp6_max<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19_firehose, what = "max")
    snp6_max=rs(snp6_max)
    snp6_min<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19_firehose, what = "min")
    snp6_min=rs(snp6_min)
    snp6=snp6_max
  }else
  {
    nstart=6
    snp6_max<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = geneInfo, what = "max")
    snp6_max=rs(snp6_max)
    snp6_min<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = geneInfo, what = "min")
    snp6_min=rs(snp6_min)
    snp6=snp6_max
  }
  
  #pick the extreme one
  for (i in nstart:ncol(snp6))
  {
    pickmin=which(abs(snp6_min[,i])-abs(snp6_max[,i])>0)
    if (length(pickmin)>0) snp6[pickmin,i]=snp6_min[pickmin,i]
  }
  
  tmp=sapply(nstart:ncol(snp6),function(x){
    tmp1=unlist(strsplit(colnames(snp6)[x],"-",fixed=TRUE))
    return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
  })
  length(unique(tmp))
  colnames(snp6)[nstart:ncol(snp6)]=tmp
  #check tumor part name
  tmp=sapply(nstart:ncol(snp6),function(x){
    tmp1=unlist(strsplit(colnames(snp6)[x],"-",fixed=TRUE))
    return(substr(tmp1[4],1,2))
  })
  print(table(tmp))
  #01  02  10  11 
  #583  18 441 130 
  idx_tumor=which(tmp %in% c("01"))
  snp6_all=snp6
  snp6=snp6_all[,nstart-1+idx_tumor]
  rownames(snp6)=snp6_all$gene
  #check patient part name
  tmp=sapply(1:ncol(snp6),function(x){
    tmp1=unlist(strsplit(colnames(snp6)[x],"-",fixed=TRUE))
    return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
  })
  print(length(tmp))
  #[1] 583
  print(length(unique(tmp)))
  if (sum(duplicated(tmp)==0))
  {
    colnames(snp6)=tmp
  }
  
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
snp6cna_hg19=readcnasegs()
snp6cna_germline_hg19=readcnasegs(segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")

copynumber_CGH1M=readcnasegs(hg19=F,segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.Level_3.2016012800.0.0/OV.cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.seg.txt")
copynumber_CGH415K=readcnasegs(hg19=F,segfile='../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg.Level_3.2016012800.0.0/OV.cna__hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg.seg.txt')
copynumber_CGH244A=readcnasegs(hg19=F,segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__hg_cgh_244a__hms_harvard_edu__Level_3__segmentation__seg.Level_3.2016012800.0.0/OV.cna__hg_cgh_244a__hms_harvard_edu__Level_3__segmentation__seg.seg.txt")
copynumber_1MDUO=readcnasegs(hg19=F,segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__human1mduo__hudsonalpha_org__Level_3__segmented_cna__seg.Level_3.2016012800.0.0/OV.snp__human1mduo__hudsonalpha_org__Level_3__segmented_cna__seg.seg.txt")


#the data from gistic result----------------------
#scna_minus_germline_cnv was used as the input segment file
gisticdata=read.table(file="../data/firehose_20160128/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt",
                      header=T,sep="\t",quote="")
rownames(gisticdata)=gisticdata$Gene.Symbol
gisticdata=gisticdata[,4:ncol(gisticdata)]
tmp=sapply(1:ncol(gisticdata),function(x){
  tmp1=unlist(strsplit(colnames(gisticdata)[x],".",fixed=TRUE))
  return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
})
length(unique(tmp))
colnames(gisticdata)=tmp
#check tumor part name
tmp=sapply(1:ncol(gisticdata),function(x){
  tmp1=unlist(strsplit(colnames(gisticdata)[x],"-",fixed=TRUE))
  return(substr(tmp1[4],1,2))
})
table(tmp)
#01
#579
sum(!colnames(gisticdata) %in% colnames(snp6hg18_median))
#[1] 0 all gistic samples were included in snp6 data
dim(gisticdata)
#[1] 24776   579
#---------------------------------------



#verify copynumber data using firehose's correlation result
#find overlap of 2 platforms,order genes/samples in the same way (as in platform 1)
#in firehose's mRNA Analyses, cor(copynumber vs mRNA), gisticdata was used as copynumber,median_mrna was used as mRNA

# 
# #old function to extract primary tumors (01) copy number gene matrix, use median
# extractcopynumber=function(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg.Level_3.2016012800.0.0/OV.cna__hg_cgh_415k_g4124a__hms_harvard_edu__Level_3__segmentation__seg.seg.txt")
# {
#   cghcopynumber=read.table(file=file,
#                                sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
#   colnames(cghcopynumber)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
#   cghcopynumber=cghcopynumber[!is.na(cghcopynumber$seg.mean),]
#   cnseg_cgh=CNSeg(cghcopynumber)
#   cgh_median <- getRS(cnseg_cgh, by = "gene", imput = FALSE, XY = TRUE, geneMap = geneInfo, what = "median")
#   cgh_median=rs(cgh_median)
#   tmp=sapply(6:ncol(cgh_median),function(x){
#     tmp1=unlist(strsplit(colnames(cgh_median)[x],"-",fixed=TRUE))
#     return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3],"-",substr(tmp1[4],1,2)))
#   })
#   length(unique(tmp))
#   colnames(cgh_median)[6:ncol(cgh_median)]=tmp
#   #check tumor part name
#   tmp=sapply(6:ncol(cgh_median),function(x){
#     tmp1=unlist(strsplit(colnames(cgh_median)[x],"-",fixed=TRUE))
#     return(substr(tmp1[4],1,2))
#   })
#   table(tmp)
#   idx_tumor=which(tmp %in% c("01"))
#   cgh_median_all=cgh_median
#   cgh_median=cgh_median_all[,5+idx_tumor]
#   #there a few duplicated gene names
#   idx_dup=which(cgh_median_all$genename %in% c("CRSP8","DUB3","PRG2","SKIP"))
#   cbind(cgh_median_all[idx_dup,1:5],idx_dup)
#   cgh_median_all$genename[21250]="CRSP8P"
#   
#   cgh_median_all[19724,1:5] #Gene ID: 402164, discontinued on 4-Apr-2009,remove it
#   cgh_median_all$genename[19724]=NA
#   cgh_median_all[12823,1:5]
#   cgh_median_all$genename[12823]="PLPPR3"
#   cgh_median_all[16110,1:5]
#   cgh_median_all$genename[16110]="SPHKAP"
#   cgh_median_all=cgh_median_all[!is.na(cgh_median_all$genename),]
#   cgh_median=cgh_median_all[,5+idx_tumor]
#   rownames(cgh_median)=cgh_median_all$genename
#   #check patient part name
#   tmp=sapply(1:ncol(cgh_median),function(x){
#     tmp1=unlist(strsplit(colnames(cgh_median)[x],"-",fixed=TRUE))
#     return(paste0(tmp1[1],"-",tmp1[2],"-",tmp1[3]))
#   })
#   length(tmp)
#   #[1] 601
#   length(unique(tmp))
#   if (length(tmp)==length(unique(tmp)))
#   {
#     cgh_median=remove_tumorid(cgh_median)
#   }
#   return(cgh_median)
# }
# copynumber_CGH1M=extractcopynumber(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.Level_3.2016012800.0.0/OV.cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.seg.txt")
# copynumber_CGH415K=extractcopynumber()
# copynumber_CGH244A=extractcopynumber(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_cna__hg_cgh_244a__hms_harvard_edu__Level_3__segmentation__seg.Level_3.2016012800.0.0/OV.cna__hg_cgh_244a__hms_harvard_edu__Level_3__segmentation__seg.seg.txt")
# copynumber_1MDUO=extractcopynumber(file="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__human1mduo__hudsonalpha_org__Level_3__segmented_cna__seg.Level_3.2016012800.0.0/OV.snp__human1mduo__hudsonalpha_org__Level_3__segmented_cna__seg.seg.txt")



overlap2platforms=function(platform1=mrna_4502a,platform2=gisticdata)
{
  comm_genes=rownames(platform1)[rownames(platform1) %in% rownames(platform2)]
  comm_samples=colnames(platform1)[colnames(platform1) %in% colnames(platform2)]
  res1=platform1[rownames(platform1) %in% comm_genes,]
  res1=res1[,colnames(platform1) %in% comm_samples]
  idx_row=sapply(comm_genes,function(x){
    res=which(rownames(platform2)==x)
  })
  idx_col=sapply(comm_samples,function(x){
    res=which(colnames(platform2)==x)
  })
  res2=platform2[idx_row,]
  res2=res2[,idx_col]
#   cor_res1_res2=sapply(1:nrow(res1),function(x){
#     res=cor(as.numeric(res1[x,]),as.numeric(res2[x,]),use="complete")
#   })
#   cor_res1_res2_order=cor_res1_res2[order(cor_res1_res2,decreasing=TRUE)]
#   names(cor_res1_res2_order)=rownames(res1)[order(cor_res1_res2,decreasing=TRUE)]
  result=list(platform1=res1,platform2=res2)
}

#only consider overlap of genes
overlap2platforms1=function(platform1=mrna_4502a,platform2=gisticdata)
{
  comm_genes=rownames(platform1)[rownames(platform1) %in% rownames(platform2)]
  #comm_samples=colnames(platform1)[colnames(platform1) %in% colnames(platform2)]
  res1=platform1[rownames(platform1) %in% comm_genes,]
  #res1=res1[,colnames(platform1) %in% comm_samples]
  idx_row=sapply(comm_genes,function(x){
    res=which(rownames(platform2)==x)
  })
  #idx_col=sapply(comm_samples,function(x){
  #  res=which(colnames(platform2)==x)
  #})
  res2=platform2[idx_row,]
  #res2=res2[,idx_col]
  #   cor_res1_res2=sapply(1:nrow(res1),function(x){
  #     res=cor(as.numeric(res1[x,]),as.numeric(res2[x,]),use="complete")
  #   })
  #   cor_res1_res2_order=cor_res1_res2[order(cor_res1_res2,decreasing=TRUE)]
  #   names(cor_res1_res2_order)=rownames(res1)[order(cor_res1_res2,decreasing=TRUE)]
  result=list(platform1=res1,platform2=res2)
}

overlap3platforms=function(platform1=median_mrna,platform2=gisticdata,platform3=mincor_methylation)
{
  tmp=overlap2platforms(platform1,platform2)
  tmp1=overlap2platforms(tmp$platform1,platform3)
  res1=tmp1$platform1
  res3=tmp1$platform2
  tmp1=overlap2platforms(res1,platform2)
  res2=tmp1$platform2
  result=list(platform1=res1,platform2=res2,platform3=res3)
}

# #save the data


#overlap of three platforms
# tmp=overlap3platforms()
# mrna_fh_com=remove_tumorid(tmp$platform1)
# copynumber_com=remove_tumorid(tmp$platform2)
# methylation_com=remove_tumorid(tmp$platform3)
# mrna_4502a=remove_tumorid(median_mrna)
# copynumber_snp6=remove_tumorid(gisticdata)
# methylation_27k=remove_tumorid(mincor_methylation)
#overlap of two platforms:
tmp=overlap2platforms(copynumber_snp6,mrna_4502a)
copynumber_comcm=tmp$platform1
mrna_comcm=tmp$platform2

tmp=overlap2platforms1(copynumber_snp6,mrna_4502a)
copynumber_comcm1=tmp$platform1
mrna_comcm1=tmp$platform2
copynumber=copynumber_snp6
mrna=mrna_4502a
methylation=methylation_27k
# save(mrna,copynumber,methylation,copynumber_comcm,mrna_comcm,file="../data/firehose.RData")
#save data in other platforms:
#save(mrna_u133,mrna_exon,copynumber_CGH1M,copynumber_1MDUO,mrnaseq,snp6cna_hg19,snp6cna_germline_hg19,file="../data/otherdata.RData")

#for oncotated mutation folder:
readmutationdata=function(mutationfolder)
{
  mutfiles=list.files(mutationfolder)
  mutfiles=mutfiles[which(grepl("TCGA-",mutfiles,fixed=T)==T)]
  mutation=data.frame(matrix(NA,ncol=10,nrow=0))
  colnames(mutation)=c("sample","gene","chr","start","end","variant_classification","variant_type","ref_allele","tumor_allele1","tumor_allele2")
  count=1
  for (i in 1:length(mutfiles))
  {
    tmptable=read.table(file=paste0(mutationfolder,"/",mutfiles[i]),header=T,sep="\t",fill=T,quote="",stringsAsFactors=F)
    for (j in 1:nrow(tmptable))
    {
      mutation[count,1]=paste0(unlist(strsplit(mutfiles[i],"-",fixed=T))[1:3],collapse="-")
      mutation[count,2:ncol(mutation)]=tmptable[j,c(1,5,6,7,9,10,11,12,13)]
      count=count+1
    }
  }
  uniq_genes=unique(mutation$gene)
  uniq_samples=unique(mutation$sample)
  mutationtable=data.frame(matrix(0,nrow=length(mutfiles),ncol=length(uniq_genes)))
  rownames(mutationtable)=uniq_samples
  colnames(mutationtable)=uniq_genes
  for (i in 1:nrow(mutation))
  {
    rowid=which(uniq_samples==mutation$sample[i])
    colid=which(uniq_genes==mutation$gene[i])
    mutationtable[rowid,colid]=mutationtable[rowid,colid]+1
  }
  return(mutationtable)
}
#save mutation 
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")

#mutation data, from raw mutation
mutationfolder="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV.Mutation_Packager_Oncotated_Raw_Calls.Level_3.2016012800.0.0/"
rawmutation=readmutationdata(mutationfolder)

mutationfolder="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/"
mutation=readmutationdata(mutationfolder)

hrgenes=c("ATM","BARD1","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FAM175A","NBN","PALB2","RAD51C","RAD51D","MRE11A")
hrgenes=hrgenes[hrgenes %in% colnames(rawmutation)]
rawmutation_hrgenes=data.frame(matrix(nrow=nrow(rawmutation),ncol=1))
rownames(rawmutation_hrgenes)=rownames(rawmutation)
colnames(rawmutation_hrgenes)="hrgenes"
idx=colnames(rawmutation) %in% hrgenes
for (i in 1:nrow(rawmutation))
{
  rawmutation_hrgenes[i,1]=sum(rawmutation[i,idx])
}
save(mrna,copynumber,methylation,copynumber_comcm,mrna_comcm,mutation,rawmutation,rawmutation_hrgenes,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")


#regress mrna on copynumber
# #salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
# njob=100
# library(Rmpi)
# #mpi.spawn.Rslaves(nslaves=4)
# #mpi.spawn.Rslaves()
# mpi.spawn.Rslaves(needlog = FALSE)
# # In case R exits unexpectedly, have it automatically clean up # resources taken up by Rmpi (slaves, memory, etc...)
# .Last <- function()
# { if (is.loaded("mpi_initialize")){
#   if (mpi.comm.size(1) > 0){
#     print("Please use mpi.close.Rslaves() to close slaves.")
#     mpi.close.Rslaves()
#   }
#   print("Please use mpi.quit() to quit R")
#   .Call("mpi_finalize")
# }
# }
# mpi.bcast.Robj2slave(mrna_comcm)
# mpi.bcast.Robj2slave(copynumber_comcm)
fitted_copynumber=residuals_copynumber=data.frame(matrix(NA,nrow=nrow(copynumber_comcm),ncol=ncol(copynumber_comcm)))
rownames(fitted_copynumber)=rownames(residuals_copynumber)=rownames((copynumber_comcm))
colnames(fitted_copynumber)=colnames(residuals_copynumber)=colnames((copynumber_comcm))
for (i in 1:nrow(copynumber_comcm))
{
  mod=lm(as.numeric(copynumber_comcm[i,])~as.matrix(mrna_comcm[i,]))
  fitted_copynumber[i,]=mod$fitted.values
  residuals_copynumber[i,]=mod$residuals
}
#save(fitted_copynumber,residuals_copynumber,file="../data/fittedcopynumber.RData")

#compare mrna,methylation,and copynumber
compute_corr=function(platform1=comp_mrna$platform1,platform2=comp_mrna$platform2)
{
  corr=sapply(1:nrow(platform1),function(x){
    res=NA
    if (sum(is.na(platform1[x,]))<ncol(platform1)-2 & sum(is.na(platform2[x,]))<ncol(platform2)-2 )
    {
      res=cor(as.numeric(platform1[x,]),as.numeric(platform2[x,]),use="complete")
    }
    return(res)
  })
  names(corr)=rownames(platform1)
  return(corr)
}

comp_mrna=overlap2platforms(platform1=mrna_fh,platform2=mrna)
corr_mrna=compute_corr(comp_mrna$platform1,comp_mrna$platform2)
hist(corr_mrna,breaks=20)
dim(comp_mrna$platform1)

comp_methylation=overlap2platforms(platform1=methylation_fh,platform2=methylation)
corr_methylation=compute_corr(comp_methylation$platform1,comp_methylation$platform2)
hist(corr_methylation,breaks=20)
dim(comp_methylation$platform1)

comp_copynumber=overlap2platforms(platform1=copynumber_snp6,platform2=copynumber)
corr_copynumber=compute_corr(comp_copynumber$platform1,comp_copynumber$platform2)
hist(corr_copynumber,breaks=20)
dim(comp_copynumber$platform1)

cmp_mrnas=overlap3platforms(mrna_fh,mrna_u133,mrna_4502a)
mrna_comgenes=rownames(cmp_mrnas$platform1)
sum(colnames(cmp_mrnas$platform1) %in% colnames(mrna))
sum(rownames(cmp_mrnas$platform1) %in% rownames(mrna))
test=data.frame(m1=cmp_mrnas$platform1[1,],m2=cmp_mrnas$platform2[1,],m3=cmp_mrnas$platform3[1,])
test1=getGeneFA(test)

getGeneFA = function (geneExprs) {
  factA = factanal(geneExprs, factors=1, scores="regression", nstart=2);  
  return (
    list(beta = factA$loading,
         Psi = factA$uniquenesses,
         corMx = factA$correlation,
         unifiedExprs = factA$scores)
  );
}

imputematbyrow=function(data)
{
  for (i in 1:nrow(data))
  {
    data[i,] <- ifelse(is.na(data[i,]),mean(data[i,!is.na(data[i,])]),data[i,])
  }
  return(data)
}
unifyexp=function(genes=rownames(mrna),samples=colnames(mrna),mrna1=mrna_fh,mrna2=mrna_u133,mrna3=mrna_4502a)
{
  mrna1=imputematbyrow(mrna1)
  mrna2=imputematbyrow(mrna2)
  mrna3=imputematbyrow(mrna3)
  #standardize
  mrna1=scale(mrna1)
  mrna2=scale(mrna2)
  mrna3=scale(mrna3)

  res=data.frame(matrix(NA,nrow=length(genes),ncol=length(samples)))
  colnames(res)=samples
  rownames(res)=genes
  #first work on samples on all 3 platforms
  tmp=samples[samples %in% colnames(mrna1)]
  tmp=tmp[tmp %in% colnames(mrna2)]
  samples_3platforms=tmp[tmp %in% colnames(mrna3)]
  idx_col=which(samples %in% samples_3platforms)
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idxcol1=which(colnames(mrna1) %in% samples_3platforms)
    idxrow1=which(rownames(mrna1)==gene)
    tmp1=mrna1[idxrow1,idxcol1]
    idxcol2=which(colnames(mrna2) %in% samples_3platforms)
    idxrow2=which(rownames(mrna2)==gene)
    tmp2=mrna2[idxrow2,idxcol2]
    idxcol3=which(colnames(mrna3) %in% samples_3platforms)
    idxrow3=which(rownames(mrna3)==gene)
    tmp3=mrna3[idxrow3,idxcol3]
    extmat=data.frame(m1=tmp1,m2=tmp2,m3=tmp3)
    tmpres=getGeneFA(extmat)
    res[i,idx_col]=tmpres$unifiedExprs
  }
  #work on samples on platforms 1 and 2
  left_samples=samples[!samples %in% samples_3platforms]
  tmp=left_samples[left_samples %in% colnames(mrna1)]
  samples_1_2=tmp[tmp %in% colnames(mrna2)]
  idx_col=which(samples %in% samples_1_2)
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idxcol1=which(colnames(mrna1) %in% samples_1_2)
    idxrow1=which(rownames(mrna1)==gene)
    tmp1=mrna1[idxrow1,idxcol1]
    idxcol2=which(colnames(mrna2) %in% samples_1_2)
    idxrow2=which(rownames(mrna2)==gene)
    tmp2=mrna2[idxrow2,idxcol2]
    extmat=data.frame(m1=tmp1,m2=tmp2)
    tmpres=getGeneFA(extmat)
    res[i,idx_col]=tmpres$unifiedExprs
  }
  left_samples1=left_samples[!left_samples %in% samples_1_2]
}

getProbeFA = function (probeExprs, nstart=5) {
  fa = try(factanal(probeExprs), factors=1, scores="regression",
           nstart=nstart);

if (class(fa) == "factanal") {
  ue = fa$scores;
  beta = fa$loadings;
  psi = diag(fa$uniquenesses)
  
  return ( list(ue=ue,
                beta=beta,
                psi=psi) )
} else {
  return (NULL)
}
}
