#load data and result:
load("../data/tcga_data.RData")
load("../data/aocs_data.RData")
load("../result/pvalues_level2data.Rdata") #probe info
load("../data/probe2gene.RData") #generated gene data using the mean of probes
source("functions.R") #keep functions compute_pvalues1 and two_df_pvalues_data

smallfdr=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t",stringsAsFactors = F)

#pvalues were calculated before, to recalculate them set recal=T
recal=F
if (recal==T)
{
  pvalues_copynumber_tangent=compute_pvalues1(data=data_copynumber_tangent) #created using tangent data. data were generated in process_TCGA_copynumber_leve2.R
  pvalues_copynumber_allele=compute_pvalues1(data=data_copynumber_allele) #created using tangent data, data were generated in process_TCGA_copynumber_BAF.R
  pvalues_2df_tangent=two_df_pvalues_data(data1=data_copynumber_tangent)
  pvalues_aocs_mean_copynumber=compute_pvalues1(data=data_aocs_mean_copynumber)
  pvalues_aocs_mean_copynumber_allele=compute_pvalues1(data=data_aocs_mean_copynumber_allele) #allele data were genereated in process_Aus_copynumber_level2.R
  pvalues_aocs_mean_2df=two_df_pvalues_data(data1=data_aocs_mean_copynumber,data2=data_aocs_mean_copynumber_allele)
}

#compute correlation between tcn and dh
cor_tcn_dh_tcga=cor_tcn_dh_data()
cor_tcn_dh_aocs=cor_tcn_dh_data(data1=data_aocs_mean_copynumber,data2=data_aocs_mean_copynumber_allele)

#count number of probes for DH data
geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
library(GenomicRanges)
countdhprobeforgenes=function(genes,type="tcga")
{
  
  if (type=="tcga")
  {
    allsamples=rownames(data_copynumber)
  }else
  {
    allsamples=rownames(data_aocs_copynumber)
  }
  
  res=data.frame(matrix(0,nrow=length(allsamples),ncol=length(genes)))
  rownames(res)=allsamples
  colnames(res)=genes
  for (i in 1:length(allsamples))
  {
    mysample=allsamples[i]
    dhprobedata=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.txt"),header=T,sep="\t",stringsAsFactors=F)
    if ("Chromosome" %in% colnames(dhprobedata))
    {
      position=as.numeric(dhprobedata$Physical.Position)
      dhprobedata=dhprobedata[!is.na(position),]
      position=position[!is.na(position)]
      gr_dhprobedata=GRanges(seqnames=dhprobedata$Chromosome,ranges=IRanges(start=position,width=1),snp=dhprobedata$Probe.Set.ID)
    }else
    {
      position=as.numeric(dhprobedata$Position)
      dhprobedata=dhprobedata[!is.na(position),]
      position=position[!is.na(position)]
      gr_dhprobedata=GRanges(seqnames=dhprobedata$Chr,ranges=IRanges(start=position,width=1),snp=dhprobedata$SNP.Name)
    }
    for (j in 1:length(genes))
    {
      gene=genes[j]
      idx=which(geneposition$gene==gene)
      gr_gene=GRanges(seqnames=geneposition$chr[idx],ranges=IRanges(start=geneposition$start[idx],end=geneposition$end[idx]))
      olap=subsetByOverlaps(gr_dhprobedata,gr_gene)
      res[i,j]=length(olap)
    }
  }

  return(res)
}

smallfdr1=cbind.data.frame(smallfdr,num_tcgaprobes=rep(0,nrow(smallfdr)),num_aocsprobes=rep(0,nrow(smallfdr)),num_tcga_dhprobes=rep(0,nrow(smallfdr)),num_aocs_dhprobes=rep(0,nrow(smallfdr)),
                           mean_cp_aocs_sense=rep(NA,nrow(smallfdr)),mean_cp_aocs_resi=rep(NA,nrow(smallfdr)),
                           pvalues_tcga_cn=rep(NA,nrow(smallfdr)),pvalues_tcga_dh=rep(NA,nrow(smallfdr)),cor_tcga_cn_dh=rep(NA,nrow(smallfdr)),pvalues_tcga_2df=rep(NA,nrow(smallfdr)),
                           pvalues_aocs_cn=rep(NA,nrow(smallfdr)),pvalues_aocs_dh=rep(NA,nrow(smallfdr)),cor_aocs_cn_dh=rep(NA,nrow(smallfdr)),pvalues_aocs_2df=rep(NA,nrow(smallfdr)),
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
  idx=which(names(cor_tcn_dh_tcga)==gene)
  if (length(idx)>0) smallfdr1$cor_tcga_cn_dh[i]=cor_tcn_dh_tcga[idx]
  idx=which(names(pvalues_2df_tangent)==gene)
  if (length(idx)>0) smallfdr1$pvalues_tcga_2df[i]=pvalues_2df_tangent[idx]
  
  idx=which(names(pvalues_aocs_mean_copynumber)==gene)
  if (length(idx)>0) smallfdr1$pvalues_aocs_cn[i]=pvalues_aocs_mean_copynumber[idx]
  idx=which(names(pvalues_aocs_mean_copynumber_allele)==gene)
  if (length(idx)>0) smallfdr1$pvalues_aocs_dh[i]=pvalues_aocs_mean_copynumber_allele[idx]
  idx=which(names(cor_tcn_dh_aocs)==gene)
  if (length(idx)>0) smallfdr1$cor_aocs_cn_dh[i]=cor_tcn_dh_aocs[idx]
  idx=which(names(pvalues_aocs_mean_2df)==gene)
  if (length(idx)>0) smallfdr1$pvalues_aocs_2df[i]=pvalues_aocs_mean_2df[idx]
  
   idx=which(names(pvalues_tcga_probe2gene)==gene)
   if (length(idx)>0) smallfdr1$pvalue_tcga_probe[i]=pvalues_tcga_probe2gene[idx]
   idx=which(names(pvalues_aocs_probe2gene)==gene)
   if (length(idx)>0) smallfdr1$pvalue_aocs_probe[i]=pvalues_aocs_probe2gene[idx]
  
}
#convert the value to copynumber
smallfdr1$cp_mean_resi=2+smallfdr1$cp_mean_resi
smallfdr1$cp_mean_sens=2+smallfdr1$cp_mean_sens

numdhprobes=countdhprobeforgenes(genes=smallfdr$gene,type="tcga")
smallfdr1$num_tcga_dhprobes=colSums(numdhprobes,na.rm=T)/nrow(numdhprobes)
numdhprobes1=countdhprobeforgenes(genes=smallfdr$gene,type="aocs")
smallfdr1$num_aocs_dhprobes=colSums(numdhprobes1,na.rm=T)/nrow(numdhprobes1)


#write.table(smallfdr1,file="../result/copynumber_result_fdrlessthan0.05_update1.txt",col.names = T,row.names = F,sep="\t",quote=F)

