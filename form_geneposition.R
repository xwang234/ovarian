#!/usr/bin/env Rscript
library(CNTools)
data("geneInfo")
#process duplicate genes
keepflag=rep(T,nrow(geneInfo))
tmp=duplicated(geneInfo$genename)
which(tmp==T)
#[1] 12823 16110 25212 27126
which(geneInfo$genename==geneInfo$genename[12823])
#[1]  4671 12823
geneInfo[c(4671,12823),]
#      chrom    start      end geneid genename
#141253    11 56911410 56914706   5553     PRG2
#380262    19   763518   772952  79948     PRG2
#manually set
keepflag[12823]=F
which(geneInfo$genename==geneInfo$genename[16110])
#[1] 10931 16110
geneInfo[c(10931,16110),]
#       chrom     start       end geneid genename
#318525    17   1344621   1366704  51763     SKIP
#480721     2 228552919 228754586  80309     SKIP
keepflag[16110]=F
which(geneInfo$genename==geneInfo$genename[25212])
#[1] 19724 25212
geneInfo[c(19724,25212),]
#       chrom    start      end geneid genename
#589695     4  8969207  8970799 402164     DUB3
#765943     8 12032086 12033678 377630     DUB3
keepflag[19724]=F
which(geneInfo$genename==geneInfo$genename[27126])
#[1] 21250 27126
geneInfo[c(21250,27126),]
#chrom     start       end geneid genename
#631614     5  79682181  79683542 441089    CRSP8
#823142     9 133725320 133945074   9442    CRSP8
keepflag[21250]=F
geneInfo=geneInfo[keepflag,]
geneInfo.bed=data.frame(matrix(NA,ncol=4,nrow=nrow(geneInfo)))
colnames(geneInfo.bed)=c("chrom","chromStart","chromEnd","name")
geneInfo.bed[,2]=geneInfo$start
geneInfo.bed[,3]=geneInfo$end
geneInfo.bed[,4]=geneInfo$genename
geneInfo.bed[,1]=paste0("chr",geneInfo$chrom)
#write.table(geneInfo.bed,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneinfohg18.bed.txt",row.names=F,col.names=F,sep="\t",quote=F)
#run liftover and get /home/xwang234/Downloads/geneInfohg19_liftover.bed
geneInfo_hg19=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg19_liftover.bed",header=F,sep="\t")
colnames(geneInfo_hg19)=colnames(geneInfo)[c(1,2,3,5)]
geneInfo_hg19$chrom=gsub("chr","",geneInfo_hg19$chrom)
#check duplicates in geneInfo genes
tmp=duplicated(geneInfo_hg19$genename)
which(tmp==T)
#save(geneInfo,geneInfo_hg19,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")


platforms=c("copynumber","copynumber_CGH1M","copynumber_1MDUO","mrna","mrna_u133","mrna_exon")
addtable=function(res,tmptable)
{
  if (nrow(tmptable)>0)
  {
    idxnew=which(! tmptable$gene %in% res$gene)
    if (length(idxnew)>0)
    {
      res=rbind(res,tmptable[idxnew,])
    }
  }
  return(res)
}
generate_geneposition=function(platforms)
{
  res=NULL
  for (platform in platforms)
  {
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_biomart.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    if (platform == platforms[1])
    {
      res=tmptable[idxnoNA,]
    }else
    {
      tmptable=tmptable[idxnoNA,]
      res=addtable(res,tmptable)
    }
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_gisticknowntable.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    tmptable=tmptable[idxnoNA,]
    res=addtable(res,tmptable)
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_ucscknowntable.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    tmptable=tmptable[idxnoNA,]
    res=addtable(res,tmptable)
  }
  return(res)
}

geneposition=generate_geneposition(platforms)

#add the genes from geneInfo
tmp=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg19_liftover.bed",header=F,sep="\t")
geneInfo_hg19=tmp[,c(4,1,2,3)]
colnames(geneInfo_hg19)=c("gene","chr","start","end")
geneInfo_hg19$chr=gsub("chr","",geneInfo_hg19$chr)


geneposition=addtable(geneInfo_hg19,geneposition)


#process alias genes, map previous alias genes to current gene symbol
mapaliasgenes=function(genes)
{
  aliastable=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/genes_multiplesymbols.txt",sep="\t",fill=T,quote="",header=T,stringsAsFactors=F )
  multiplegenes=paste(aliastable$symbol,aliastable$alias_symbol,sep="|")
  multiplegenes=gsub("\\|$","",multiplegenes,perl=T)
  multiplegenes=gsub("\"","",multiplegenes,perl=T)
  multiplegenes=gsub("\'","",multiplegenes,perl=T)
  
  genes=as.character(genes)
  res=genes
  
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idx=which(grepl(paste0("\\<",gene,"\\>"),multiplegenes)==T)
    if (length(idx)==0)
    {
      idx=which(grepl(paste0("\\<",tolower(gene),"\\>"),tolower(multiplegenes))==T)
    }
    if (length(idx)>0)
    {
      res[i]=unlist(strsplit(multiplegenes[idx[1]],"|",fixed=T))[1]
    }
  }
  if (sum(duplicated(res))>0)
  {
    warning("there are duplicated gene names after transformation")
  }
  return(res)
}
#geneposition has 39407 genes, 36938 are unique genes (after removing alias genes). Genes with same letters can be different genes. RGR and Rgr (alias of RGL4) are different genes.


#form a aliastable gene-symbol. symbol is the current gene symbol name
aliastable=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/genes_multiplesymbols.txt",sep="\t",fill=T,quote="",header=T,stringsAsFactors=F )
multiplegenes=paste(aliastable$symbol,aliastable$alias_symbol,sep="|")
multiplegenes=gsub("\\|$","",multiplegenes,perl=T)
multiplegenes=gsub("\"","",multiplegenes,perl=T)
multiplegenes=gsub("\'","",multiplegenes,perl=T)
aliasgenetable=matrix(NA,ncol=2,nrow=0)
for (i in 1:length(multiplegenes))
{
  tmp=unlist(strsplit(multiplegenes[i],"|",fixed=T))
  for (j in 1:length(tmp))
  {
    aliasgenetable=rbind(aliasgenetable,c(tmp[j],tmp[1]))
  }
}
colnames(aliasgenetable)=c("gene","symbol")
aliasgenetable=as.data.frame(aliasgenetable)
aliasgenetable[,1]=as.character(aliasgenetable[,1])
aliasgenetable[,2]=as.character(aliasgenetable[,2])

#add symbol into geneposition table
genesymbolposition=data.frame(matrix(NA,ncol=5,nrow=nrow(geneposition)))
genesymbolposition[,c(1,3,4,5)]=geneposition
colnames(genesymbolposition)[c(1,3,4,5)]=colnames(geneposition)
colnames(genesymbolposition)[2]="symbol"
genesymbolposition[,2]=geneposition[,1]
genesymbolposition[,1]=as.character(genesymbolposition[,1])
genesymbolposition[,2]=as.character(genesymbolposition[,2])
for (i in 1:nrow(genesymbolposition))
{
  idx=which(aliasgenetable$gene==geneposition$gene[i])
  if (length(idx)>0)
  {
    genesymbolposition[i,2]=aliasgenetable$symbol[idx[1]]
  }
}

genesymbolposition$chr=gsub("chr","",genesymbolposition$chr)
genesymbolposition=sortgenetable(genesymbolposition)

save(aliasgenetable,genesymbolposition,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genesymbolposition.RData")
# #check width
# tmp=geneposition$end-geneposition$start
# tmp1=which(tmp>1000000)
# 
# 
# #check duplicated genes
# keepflag=rep(T,nrow(geneposition))
# tmp=duplicated(geneposition$gene)
# which(tmp==T)
# #[1] 12810 16095 25154 27066
# which(geneposition$gene==geneposition$gene[12810])
# #[1]  4665 12810
# geneposition[c(4665,12810),]
# #      gene chr    start      end
# #4665  PRG2  11 57154834 57158130
# #12810 PRG2  19   812518   821952
# keepflag[12810]=F
# which(geneposition$gene==geneposition$gene[16095])
# #[1] 10924 16095
# geneposition[c(10924,16095),]
# #      gene chr     start       end
# #10924 SKIP  17   1397871   1419954
# #16095 SKIP   2 228844675 229046342
# keepflag[16095]=F
# which(geneposition$gene==geneposition$gene[25154])
# #[1] 19709 25154
# geneposition[c(19709,25154),]
# #      gene chr    start      end
# #19709 DUB3   4  9360109  9361701
# #25154 DUB3   8 11994677 11996269
# keepflag[19709]=F
# which(geneposition$gene==geneposition$gene[27066])
# #[1] 21225 27066
# geneposition[c(21225,27066),]
# #      gene chr     start       end
# #21225 CRSP8   5  79646425  79647786
# #27066 CRSP8   9 134735499 134955253
# keepflag[21225]=F
# geneposition=geneposition[keepflag,]

# geneposition$chr=gsub("chr","",geneposition$chr)
# geneposition=sortgenetable(geneposition)
# write.table(geneposition,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",row.names=F,col.names=T,sep="\t",quote=F)
