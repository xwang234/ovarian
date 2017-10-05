#!/usr/bin/env Rscript

source("functions.R")
#first re-compute all the pvalues
#for aocs:

load("../data/aocs_data.RData")

#for only platinumclass added data (1st column is platinumclass)
compute_pvalues1=function(data)
{
  pvalues=sapply(2:ncol(data),function(i){
    fm=glm(as.factor(data[,"platinumclass"])~as.numeric(as.character(data[,i])),family = "binomial")
    res=NA
    if (nrow(coef(summary(fm)))>1)
    {
      res=coef(summary(fm))[2,4]
    }
    names(res)=colnames(data)[i-1]
    return(res)
  })
  return(pvalues)
}

pvalues_aocs_copynumber=compute_pvalues1(data=data_aocs_copynumber)
pvalues_aocs_copynumber_filtered=compute_pvalues1(data=data_aocs_copynumber_filtered)
pvalues_aocs_mean_copynumber=compute_pvalues1(data=data_aocs_mean_copynumber)
pvalues_aocs_mean_copynumber_allele=compute_pvalues1(data=data_aocs_mean_copynumber_allele)
pvalues_aocs_mean_copynumber_allele_rmcnv=compute_pvalues1(data=data_aocs_mean_copynumber_allele_rmcnv)

# save(allsamples,primarysamplefiles,aocs_platinumclass,omni_1.0,omni_1.1,probes2keep,probes2keep_rmcnv,
# GSM_AOCS_ID,AOCS_GEO_ID,aocs_copynumber,aocs_mean_copynumber,aocs_mean_copynumber_rmcnv,aocs_mean_copynumber_allele,aocs_mean_copynumber_allele_rmcnv,
# data_aocs_copynumber,data_aocs_copynumber_filtered,data_aocs_mean_copynumber,data_aocs_mean_copynumber_allele,data_aocs_mean_copynumber_allele_rmcnv,
# pvalues_aocs_2df,pvalues_aocs_mean_2df,pvalues_aocs_mean_copynumber_rmcnv_2df,
# pvalues_aocs_copynumber,pvalues_aocs_copynumber_filtered,pvalues_aocs_mean_copynumber,pvalues_aocs_mean_copynumber_allele,pvalues_aocs_mean_copynumber_allele_rmcnv,
# file="../data/aocs_data.RData")

#for TCGA data
load("../data/tcga_data.RData")
pvalues_copynumber=compute_pvalues1(data=data_copynumber)
pvalues_copynumber_filtered=compute_pvalues1(data=data_copynumber_filtered)
pvalues_copynumber_tangent=compute_pvalues1(data=data_copynumber_tangent)
pvalues_copynumber_tangent_filtered=compute_pvalues1(data=data_copynumber_tangent_filtered)
pvalues_copynumber_tangent_gistic=compute_pvalues1(data=data_copynumber_tangent_gistic)
pvalues_copynumber_tangent_gistic_filtered=compute_pvalues1(data=data_copynumber_tangent_gistic_filtered)

# save(filelist,allsamples,tcga_platinumclass,anno1,anno1_rmcnv,copynumber,copynumber_tangent,copynumber_tangent_gistic,copynumber_firehosesegment,copynumber_firehosesegment_gistic,
#      copynumber_tangent_rmcnv,copynumber_allele,data_copynumber_allele,pvalues_copynumber_allele,
#      data_copynumber,data_copynumber_filtered,data_copynumber_firehosesegment,data_copynumber_firehosesegment_gistic,data_copynumber_tangent,
#      data_copynumber_tangent_filtered,data_copynumber_tangent_gistic,data_copynumber_tangent_gistic_filtered,data_copynumber_tangent_rmcnv,
#      pvalues_2df,pvalues_2df_tangent,pvalues_2df_tangent_gistic,pvalues_2df_tangent_rmcnv,
#      pvalues_copynumber,pvalues_copynumber_filtered,pvalues_copynumber_tangent,pvalues_copynumber_tangent_filtered,pvalues_copynumber_tangent_gistic,pvalues_copynumber_tangent_gistic_filtered,
#      file="../data/tcga_data.RData")


#load("../result/tmp_pvalues.RData")

dataframe2vector=function(data)
{
  if (class(data)=="data.frame")
  {
    res=data[,1]
    names(res)=rownames(data)
  }else
  {
    res=data
  }
  return(res)
}

geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)

#for aocs:
df_probes2keep=data.frame(Name=probes2keep)
aocs_anno=merge(df_probes2keep,omni_1.1)
aocs_anno=aocs_anno[,1:3]
colnames(aocs_anno)=c("probe","chr","start")
myfile=primarysamplefiles[1]
aocs_probedata=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/GSE65821/",myfile),sep="\t",header=T,skip=10,
                     stringsAsFactors=F)
aocs_probedata=aocs_probedata[,c("SNP.Name","Log.R.Ratio")]
colnames(aocs_probedata)=c("probe","value")
aocs_genes=colnames(data_aocs_mean_copynumber)[2:ncol(data_aocs_mean_copynumber)]
probes_in_genes=function(probedata=aocs_probedata,anno=aocs_anno,genes=aocs_genes)
{
  probedata=probedata[probedata$probe %in% anno$probe,]
  probedata=merge(anno,probedata)
  
  
  
}
  
  

    
comparepvalues_2data=function(pvalues1=pvalues_2df,pvalues2=pvalues_aocs_2df,cutoff=0.01)
{
  pvalues1=dataframe2vector(pvalues1)
  pvalues2=dataframe2vector(pvalues2)
  commongenes=intersect(names(pvalues1),names(pvalues2))
  df_genes=data.frame(gene=commongenes)
  df_pvalues1=data.frame(gene=names(pvalues1),pvalue1=pvalues1)
  df_pvalues1=merge(df_genes,df_pvalues1)
  df_pvalues2=data.frame(gene=names(pvalues2),pvalue2=pvalues2)
  df_pvalues=merge(df_pvalues1,df_pvalues2)
  idx=which(df_pvalues$pvalue1<=cutoff & df_pvalues$pvalue2<=cutoff)
  res=NA
  if (length(idx)>0)
  {
    res=df_pvalues[idx,]
  }
  return(res)
}

comp_pvalues_tcga_aocs=comparepvalues_2data(pvalues1=pvalues_copynumber,pvalues2=pvalues_aocs_copynumber)

