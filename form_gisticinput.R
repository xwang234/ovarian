#!/usr/bin/env Rscript

#update cilinical items of "clinical stage" and "residual disease"
updateclinicalitem=function(data=data_mrna)
{
  if (is.data.frame(data)==F)
  {
    alldata=NULL
    for (i in 1:length(data))
    {
      alldata=rbind(alldata,data[[i]])
    }
  }else
  {
    alldata=data
  }
  
  #process clinical stage, merge A,B,C
  alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IA","Stage IB","Stage IC"))
  alldata[tmp,"clinical_stage"]="Stage1"
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IIA","Stage IIB","Stage IIC"))
  alldata[tmp,"clinical_stage"]="Stage2"
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IIIA","Stage IIIB","Stage IIIC"))
  alldata[tmp,"clinical_stage"]="Stage3"
  tmp=which(alldata[,"clinical_stage"] %in% c("Stage IV"))
  alldata[tmp,"clinical_stage"]="Stage4"
  table(alldata[,"clinical_stage"])
  #Stage1 Stage2 Stage3 Stage4 
  #8     18    300     54 
  #process residual names
  alldata[,"residual_disease_largest_nodule"]=as.character(alldata[,"residual_disease_largest_nodule"])
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "No Macroscopic disease")
  alldata[tmp,"residual_disease_largest_nodule"]="No_Macroscopic_disease"
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "1-10 mm")
  alldata[tmp,"residual_disease_largest_nodule"]="1_10mm"
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% "11-20 mm")
  alldata[tmp,"residual_disease_largest_nodule"]="11_20mm"
  tmp=which(alldata[,"residual_disease_largest_nodule"] %in% ">20 mm")
  alldata[tmp,"residual_disease_largest_nodule"]="20mm_"
  alldata[,"residual_disease_largest_nodule"]=as.factor(alldata[,"residual_disease_largest_nodule"])
  return(alldata)
}

#segtable=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.txt",header=T,sep="\t")
segtable=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",header=T,sep="\t")

tmp=as.character(segtable$Sample)
tmp1=sapply(1:length(tmp),function(i){
  tmp2=unlist(strsplit(tmp[i],"-",fixed=T))
  res=paste0(tmp2[4])
})
tmp2=grepl("01",tmp1)
segtable1=segtable[tmp2,]
tmp=as.character(segtable1$Sample)
tmp1=sapply(1:length(tmp),function(i){
  tmp2=unlist(strsplit(tmp[i],"-",fixed=T))
  res=paste0(tmp2[4])
})
idx_01A=which(tmp1=="01A")
idx_01B=which(tmp1=="01B")
View(segtable1[idx_01A,])
View(segtable1[idx_01B,])
tmp=as.character(segtable1$Sample)
tmp1=sapply(1:length(tmp),function(i){
  tmp2=unlist(strsplit(tmp[i],"-",fixed=T))
  res=paste0(tmp2[1:3],collapse = "-")
})
segtable2=segtable1
segtable2$Sample=tmp1

alldata=rbind(data_mrna$data1,data_mrna$data2)
alldata=updateclinicalitem(alldata)
#remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
alldata[,"tumor_grade"]=as.character(alldata[,"tumor_grade"])
alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]
alldata[,"clinical_stage"]=as.character(alldata[,"clinical_stage"])

sensitivesamples=rownames(alldata)[which(alldata[,"platinumclass"]=="Sensitive")]
resistantsamples=rownames(alldata)[which(alldata[,"platinumclass"]=="Resistant")]
sensitiveseg=segtable2[which(segtable2$Sample %in% sensitivesamples),]
resistantseg=segtable2[which(segtable2$Sample %in% resistantsamples),]
write.table(sensitiveseg,file="../data/snp6copynumber_sensitiveseg.txt",sep="\t",col.names = F,row.names = F,quote=F)
write.table(resistantseg,file="../data/snp6copynumber_resistantseg.txt",sep="\t",col.names = F,row.names = F,quote=F)

#for uspersensitive and refractory
supersensitiveseg=segtable2[which(segtable2$Sample %in% supersensitivesamples),]
refractoryseg=segtable2[which(segtable2$Sample %in% refractorysamples),]
write.table(supersensitiveseg,file="../data/snp6copynumber_supersensitiveseg.txt",sep="\t",col.names = F,row.names = F,quote=F)
write.table(refractoryseg,file="../data/snp6copynumber_refractoryseg.txt",sep="\t",col.names = F,row.names = F,quote=F)

#work on supersensitive/refractory data, #find original cnv segment in certain region
library(GenomicRanges)
findregion=function(segtable=supersensitiveseg,region)
{
  segtable$Sample=as.character(segtable$Sample)
  tmp=unlist(strsplit(region,":",fixed=T))
  tmp1=substr(tmp[1],4,nchar(tmp[1]))
  tmp2=unlist(strsplit(tmp[2],"-",fixed=T))
  gr_region=GRanges(seqnames = tmp1,ranges=IRanges(start=as.integer(tmp2[1]),end=as.integer(tmp2[2])))
  samples=unique(as.character(segtable$Sample))
  res=data.frame(matrix(NA,ncol=3,nrow=length(samples)))
  colnames(res)=c("sample","count","amplitude")
  res$sample=samples
  for (i in 1:length(samples))
  {
    subsegtable=segtable[segtable$Sample==samples[i],]
    gr_subsegtable=GRanges(seqnames = subsegtable$Chromosome,ranges=IRanges(start=subsegtable$Start,end=subsegtable$End),amplitude=subsegtable$Segment_Mean)
    olap=subsetByOverlaps(gr_subsegtable,gr_region)
    if (length(olap)>0)
    {
      res$count[i]=length(olap)
      res$amplitude[i]=max(mcols(olap)$amplitude)
    }
  }
  return(res)
}

res_super1=findregion(segtable=supersensitiveseg,region="chr5:14968148-15000267")
res_refra1=findregion(segtable=refractoryseg,region="chr5:14968148-15000267")


#work on supersensitive/refractory data, use cnv called by gistic
supersensitivecall=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/supersensitive_rx0_conf0.95_armpeel0_brlen0.98_broad1/all_thresholded.by_genes.txt",header=T,sep="\t")
refractorycall=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/refractory_rx0_conf0.95_armpeel0_brlen0.98_broad1/all_thresholded.by_genes.txt",header=T,sep="\t")
findregion1=function(cnvcall=supersensitivecall,cytoband="5p15.2")
{
  cnvcall$Cytoband=as.character(cnvcall$Cytoband)
  subtable=cnvcall[cnvcall$Cytoband==cytoband,]
  numsamples=ncol(cnvcall)-3
  res=list(numsamples=numsamples,numamp=0,numdel=0)
  if (nrow(subtable)>0)
  {
    res$numamp=sum(subtable[1,4:ncol(subtable)]>=1)
    res$numdel=sum(subtable[1,4:ncol(subtable)]<0)
  }
  return(res)
}
res_super_5p15.2=findregion1(cnvcall=supersensitivecall,cytoband="5p15.2")
res_refra_5p15.2=findregion1(cnvcall=refractorycall,cytoband="5p15.2")

#should work on focal data, the above used all_data
supersensitivecall=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/supersensitive_rx0_conf0.95_armpeel0_brlen0.98_broad1/focal_data_by_genes_sups.txt",header=T,sep="\t")
refractorycall=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/refractory_rx0_conf0.95_armpeel0_brlen0.98_broad1/focal_data_by_genes_refr.txt",header=T,sep="\t")
findregion2=function(cnvcall=supersensitivecall,cytoband="5p15.2")
{
  cnvcall$Cytoband=as.character(cnvcall$Cytoband)
  subtable=cnvcall[cnvcall$Cytoband==cytoband,]
  numsamples=ncol(cnvcall)-3
  res=list(numsamples=numsamples,numamp=0,numdel=0)
  if (nrow(subtable)>0)
  {
    res$numamp=sum(subtable[1,4:ncol(subtable)]>=0.1)
    res$numdel=sum(subtable[1,4:ncol(subtable)]<=-0.1)
  }
  return(res)
}
res_super_5p15.2=findregion2(cnvcall=supersensitivecall,cytoband="5p15.2")
res_refra_5p15.2=findregion2(cnvcall=refractorycall,cytoband="5p15.2")

res_super_7q34=findregion2(cnvcall=supersensitivecall,cytoband="7q34")
res_refra_7q34=findregion2(cnvcall=refractorycall,cytoband="7q34")

res_super_7q36.3=findregion2(cnvcall=supersensitivecall,cytoband="7q36.3")
res_refra_7q36.3=findregion2(cnvcall=refractorycall,cytoband="7q36.3")

res_super_10q22.2=findregion2(cnvcall=supersensitivecall,cytoband="10q22.2")
res_refra_10q22.2=findregion2(cnvcall=refractorycall,cytoband="10q22.2")

library(MASS)
proptest=function(res1=res_super_5p15.2,res2=res_refra_5p15.2)
{
  tmp=matrix(c(res1$numamp,res1$numsamples-res1$numamp,res2$numamp,res2$numsamples-res2$numamp),nrow=2,byrow = T)
  tmp1=prop.test(tmp)
  res=tmp1$p.value
}
pvalue_5p15.2=proptest()
pvalue_7q34=proptest(res1=res_super_7q34,res2=res_refra_7q34)
pvalue_7q36.3=proptest(res1=res_super_7q36.3,res2=res_refra_7q36.3)
pvalue_10q22.2=proptest(res1=res_super_10q22.2,res2=res_refra_10q22.2)
# #use markerfile genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt
# #form marker file
# sensitivemarker1=sensitivemarker2=data.frame(matrix(NA,nrow=nrow(sensitiveseg),ncol=3))
# colnames(sensitivemarker1)=colnames(sensitivemarker2)=c("unitName","chromosome","position")
# tmp=sapply(1:nrow(sensitiveseg),function(i){
#   tmp1=paste0("chr",paste0(sensitiveseg[i,2:3],collapse ="-"))
# })
# sensitivemarker1$unitName=tmp
# sensitivemarker1$chromosome=sensitiveseg$Chromosome
# sensitivemarker1$position=sensitiveseg$Start
# 
# tmp=sapply(1:nrow(sensitiveseg),function(i){
#   tmp1=paste0("chr",paste0(sensitiveseg[i,c(2,4)],collapse ="-"))
# })
# sensitivemarker2$unitName=tmp
# sensitivemarker2$chromosome=sensitiveseg$Chromosome
# sensitivemarker2$position=sensitiveseg$End
# sensitivemarker=rbind(sensitivemarker1,sensitivemarker2)
# #remove duplicated markers
# tmp=duplicated(sensitivemarker$unitName)
# sensitivemarker=sensitivemarker[!tmp,]
# 
# write.table(sensitivemarker,file="../data/snp6copynumber_sensitiveseg_marker.txt",sep="\t",row.names = F,col.names = T,quote=F)
# 
# resistantmarker1=resistantmarker2=data.frame(matrix(NA,nrow=nrow(resistantseg),ncol=3))
# colnames(resistantmarker1)=colnames(resistantmarker2)=c("unitName","chromosome","position")
# tmp=sapply(1:nrow(resistantseg),function(i){
#   tmp1=paste0("chr",paste0(resistantseg[i,2:3],collapse ="-"))
# })
# resistantmarker1$unitName=tmp
# resistantmarker1$chromosome=resistantseg$Chromosome
# resistantmarker1$position=resistantseg$Start
# 
# tmp=sapply(1:nrow(resistantseg),function(i){
#   tmp1=paste0("chr",paste0(resistantseg[i,c(2,4)],collapse ="-"))
# })
# resistantmarker2$unitName=tmp
# resistantmarker2$chromosome=resistantseg$Chromosome
# resistantmarker2$position=resistantseg$End
# resistantmarker=rbind(resistantmarker1,resistantmarker2)
# #remove duplicated markers
# tmp=duplicated(resistantmarker$unitName)
# resistantmarker=resistantmarker[!tmp,]
# 
# write.table(resistantmarker,file="../data/snp6copynumber_resistantseg_marker.txt",sep="\t",row.names = F,col.names = T,quote=F)
