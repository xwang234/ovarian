#!/usr/bin/env Rscript
#sample sizes of 3 copy number platforms:
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")
samples_snp6=c(rownames(data_copynumber$data1),rownames(data_copynumber$data2))
samples_CGH1M=c(rownames(data_copynumber_CGH1M$data1),rownames(data_copynumber_CGH1M$data2))
samples_1MDUO=c(rownames(data_copynumber_1MDUO$data1),rownames(data_copynumber_1MDUO$data2))
area1=length(samples_snp6)
#[1] 381
area2=length(samples_CGH1M)
#[1] 385
area3=length(samples_1MDUO)
#[1] 354
n12=sum(samples_snp6 %in% samples_CGH1M)
#[1] 379
n13=sum(samples_snp6 %in% samples_1MDUO)
#[1] 350
n23=sum(samples_CGH1M %in% samples_1MDUO)
#[1] 353
n123=sum(samples_snp6[samples_snp6 %in% samples_CGH1M] %in% samples_1MDUO)
#[1] 349
library(VennDiagram)
grid.newpage();
venn.plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("SNP6", "CGH1M", "1MDUO"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 3,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
);
grid.draw(venn.plot);

source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
generatesmallfdrtable=function(fdrs,data=data_copynumber,threshold=0.05)
{
  genenames=rownames(fdrs)
  fdrs=fdrs[,2] #fdr column
  names(fdrs)=genenames
  res=smallfdrtable(fdrs=fdrs,data=data,threshold=threshold,includecor=F)
}
fdrs1000p_copynumber=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_1000permutation.txt")
smallfdr1000p_copynumber=generatesmallfdrtable(fdrs=fdrs1000p_copynumber)
fdrs1000p_copynumber_CGH1M=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_CGH1M_1000permutation.txt")
smallfdr1000p_copynumber_CGH1M=generatesmallfdrtable(fdrs=fdrs1000p_copynumber_CGH1M,data=data_copynumber_CGH1M)
fdrs1000p_copynumber_1MDUO=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_1MDUO_1000permutation.txt")
smallfdr1000p_copynumber_1MDUO=generatesmallfdrtable(fdrs=fdrs1000p_copynumber_1MDUO,data=data_copynumber_1MDUO)

#no fdr are < 0.05 in 2 other platforms, check the pvalues
genes=as.character(smallfdr1000p_copynumber$gene)

extract_gene_pvalues=function(pvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_CGH1M_platinum.txt",genes)
{
  pvalues=read.table(file=pvaluefile)
  res=data.frame(matrix(NA,nrow=length(genes),ncol=2))
  colnames(res)=c("pvalue","rank")
  pvalues_order=pvalues[order(pvalues[,1]),1]
  names(pvalues_order)=rownames(pvalues)[order(pvalues[,1])]
  for (i in 1:length(genes))
  {
    idx=which(rownames(pvalues)==genes[i])
    if (length(idx)>0)
    {
      res[i,1]=pvalues[idx[1],1]
      theorder=which(names(pvalues_order)==rownames(pvalues)[idx[1]])
      res[i,2]=min(theorder)
      rownames(res)[i]=rownames(pvalues)[idx[1]]
    }
  }
  return(res)
}

pvaluestable=extract_gene_pvalues(pvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_platinum.txt",genes)
tmp=extract_gene_pvalues(pvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_CGH1M_platinum.txt",genes)
pvaluestable=cbind(pvaluestable,tmp)
tmp=extract_gene_pvalues(pvaluefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1MDUO_platinum.txt",genes)
pvaluestable=cbind(pvaluestable,tmp)
tmp=c(rep("snp6_",2),rep("CGH1M_",2),rep("1MDUO_",2))
colnames(pvaluestable)=paste0(tmp,colnames(pvaluestable))
write.table(pvaluestable,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/smallfdrs_snp6_pvalues_3platforms.txt",row.names=T,col.names = T,sep="\t",quote=F)

#check relationship between cna and mrna
copynumber_comm=rbind(data_copynumber_commc$data1,data_copynumber_commc$data2)
mrna_comm=rbind(data_mrna_commc$data1,data_mrna_commc$data2)
#select mrna cutoff based on all genes
tmp=quantile(unlist(mrna_comm[,14:ncol(mrna_comm)]),c(0.25,0.75))
mrna_l=tmp[1]
mrna_h=tmp[2]
# tmp=quantile(unlist(copynumber_comm),c(0.25,0.75))
# cp_l=tmp[1]
# cp_h=tmp[2]
cp_l=-0.5
cp_h=0.5

check_mrna_copynumber=function(gene,mrna_h,mrna_l,cp_hi,cp_l)
{
  idx=which(colnames(mrna_comm)==gene)
  res=data.frame(matrix(NA,ncol=2,nrow=3))
  colnames(res)=paste0("mrna_",c("high","low"))
  rownames(res)=paste0(gene,"_",c("del","interm","amp"))
  if (length(idx)>0)
  {
    res[1,1]=sum(mrna_comm[,idx]>=mrna_h & copynumber_comm[,idx]<=cp_l)
    res[1,2]=sum(mrna_comm[,idx]<=mrna_l & copynumber_comm[,idx]<=cp_l)
    res[2,1]=sum(mrna_comm[,idx]>=mrna_h & copynumber_comm[,idx]<cp_h & copynumber_comm[,idx]>cp_l)
    res[2,2]=sum(mrna_comm[,idx]<=mrna_l & copynumber_comm[,idx]<cp_h & copynumber_comm[,idx]>cp_l)
    res[3,1]=sum(mrna_comm[,idx]>=mrna_h & copynumber_comm[,idx]>=cp_h)
    res[3,2]=sum(mrna_comm[,idx]<=mrna_l & copynumber_comm[,idx]>=cp_h)
  }
  return(res)
}

mrna_vs_copynumber=NULL
for (gene in genes)
{
  mrna_vs_copynumber=rbind(mrna_vs_copynumber,check_mrna_copynumber(gene,mrna_h=mrna_h,mrna_l=mrna_l,cp_hi=cp_hi,cp_l=cp_l))
}

#use mrna_l mrna_h specified for each gene
mrna_vs_copynumber1=NULL
for (gene in genes)
{
  idx=which(colnames(mrna_comm)==gene)
  if (length(idx)>0)
  {
    tmp=quantile(mrna_comm[,idx],c(0.25,0.75))
    mrna_l1=tmp[1]
    mrna_h1=tmp[2]
  }else
  {
    mrna_l1=mrna_l
    mrna_h1=mrna_h
  }
  
  mrna_vs_copynumber1=rbind(mrna_vs_copynumber1,check_mrna_copynumber(gene,mrna_h=mrna_h1,mrna_l=mrna_l1,cp_hi=cp_hi,cp_l=cp_l))
}

#differential expression
mrna_difexp=data.frame(matrix(NA,ncol=4,nrow=length(genes)))
colnames(mrna_difexp)=c("fdr","sens_mean","resi_mean","pvalue")
mrna_difexp$fdr=smallfdr1000p_copynumber$fdr
rownames(mrna_difexp)=genes
for (i in 1:length(genes))
{
  gene=genes[i]
  idx=which(colnames(mrna_comm)==gene)
  if (length(idx)>0)
  {
    sens=mrna_comm[which(mrna_comm[,"platinumclass"]=="Sensitive"),idx]
    resi=mrna_comm[which(mrna_comm[,"platinumclass"]=="Resistant"),idx]
    mrna_difexp[i,2:4]=c(mean(sens,na.rm=T),mean(resi,na.rm=T),t.test(sens,resi)$p.value)
    
  }
}
smallfdr1000p_copynumber=cbind(smallfdr1000p_copynumber,mrna_difexp[,2:4])
colnames(smallfdr1000p_copynumber)[9:11]=paste0("cp_",colnames(smallfdr1000p_copynumber)[9:11])
colnames(smallfdr1000p_copynumber)[13:15]=paste0("mrna_",colnames(smallfdr1000p_copynumber)[13:15])
#include mutation data
compare_sens_resi=function(data=data_mutation,genes=as.character(smallfdr1000p_copynumber$gene),opt="sum")
{
  alldata=formwholedataframe(data)
  #remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
  alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
  alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
  alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  idx_y=which(colnames(alldata)=="platinumclass")
  resistantdata=alldata[alldata[,"platinumclass"]=="Resistant",(idx_y+3):ncol(alldata)]
  sensitivedata=alldata[alldata[,"platinumclass"]=="Sensitive",(idx_y+3):ncol(alldata)]
  mean_res=data.frame(mean_resi=rep(NA,length(genes)),mean_sens=rep(NA,length(genes)),pvalue=rep(NA,length(genes)))
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idx=which(colnames(resistantdata)==gene)
    if (length(idx)>0)
    {
      if (opt=="sum")
      {
        mean_res[i,1]=sum(resistantdata[,idx]>0)
        mean_res[i,2]=sum(sensitivedata[,idx]>0)
        tmp=c(sum(resistantdata[,idx]>0),nrow(resistantdata)-sum(resistantdata[,idx]>0))
        tmp=rbind(tmp,c(sum(sensitivedata[,idx]>0),nrow(sensitivedata)-sum(sensitivedata[,idx]>0)))
        #mean_res[i,3]=prop.test(tmp)$p.value
        mean_res[i,3]=fisher.test(tmp)$p.value
      }else # for mean
      {
        mean_res[i,1]=mean(resistantdata[,idx])
        mean_res[i,2]=mean(sensitivedata[,idx])
        mean_res[i,3]=t.test(resistantdata[,idx],sensitivedata[,idx])$p.value
      }
    }
  }
  return(mean_res)
}
smallfdr_mutation=compare_sens_resi()
smallfdr1000p_copynumber=cbind(smallfdr1000p_copynumber,smallfdr_mutation)
colnames(smallfdr1000p_copynumber)[16:18]=paste0("mutation_",colnames(smallfdr1000p_copynumber)[16:18])

write.table(smallfdr1000p_copynumber,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_result_fdrlessthan0.05.txt",row.names = F,col.names = T,sep="\t",quote=F)
