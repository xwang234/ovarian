#!/usr/bin/env Rscript

library(gdata)
clinical_4classification=read.xls("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_4classification_table_jl1.xlsx",1,sep=",")
clinical_4classification$corrected_classification=as.character(clinical_4classification$corrected_classification)
clinical_classification=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/clinical_vairables.txt",header=T,sep="\t",stringsAsFactors=F)
clinical_classification$platinumclass=as.character(clinical_classification$platinumclass)

idx=sapply(clinical_classification$sample,function(mysample){
  which(clinical_4classification$sample==mysample)
})
clinical_4classification=clinical_4classification[idx,]
#check the two-class def are identical in the two tables, yes.
sum(clinical_classification$platinumclass==clinical_4classification$corrected_classification,na.rm = T)
idx=NULL
for (i in 1:nrow(clinical_classification))
{
  if (! is.na(clinical_classification$platinumclass[i]))
  if (clinical_classification$platinumclass[i] != clinical_4classification$corrected_classification[i])
  {
    idx=c(idx,i)
  }
}

#check the result if use the first treatment date
res=data.frame(matrix(NA,nrow=nrow(clinical_4classification),ncol=3))
colnames(res)=c("new_druginterval","new_platinumclass","change")
