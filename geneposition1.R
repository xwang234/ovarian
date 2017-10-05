#!/usr/bin/env Rscript

platforms=c("copynumber","copynumber_CGH1M","copynumber_1MDUO","mrna","mrna_exon","mrna_u133")
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
    #tmptable=tmptable[idxnoNA,]
    res=addtable(res,tmptable)
    tmptable=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_ucscknowntable.txt"),header=T)
    idxnoNA=which(!is.na(tmptable$chr))
    #tmptable=tmptable[idxnoNA,]
    res=addtable(res,tmptable)
  }
  return(res)
}
genealiasposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/aliasgeneposition_biomart.txt",header=T)
filteraliasgene=function(geneposition,genealiasposition)
{
  geneposition$gene=as.character(geneposition$gene)
  genealiasposition$gene=as.character(genealiasposition$gene)
  geneposition$chr=as.character(geneposition$chr)
  genealiasposition$chr=as.character(genealiasposition$chr)
  #for alias gene aaa|bbb|ccc, change bbb/ccc to aaa
  for (i in 1:nrow(geneposition))
  {
    tmp=which(grepl(paste0("\\<",geneposition$gene[i],"\\>"),genealiasposition$gene)==T)
    if (length(tmp)>0)
    {
      geneposition$gene[i]=unlist(strsplit(genealiasposition$gene[tmp[1]],"|",fixed=T))[1]
    }
  }
  uniqgenes=unique(geneposition$gene)
  uniqgenetable=data.frame(matrix(NA,nrow=length(uniqgenes),ncol=ncol(geneposition)))
  colnames(uniqgenetable)=colnames(geneposition)
  uniqgenetable$gene=uniqgenes
  for (i in 1:length(uniqgenes))
  {
    tmptable=geneposition[which(geneposition$gene==uniqgenes[i]),]
    idx=which(!is.na(tmptable$chr)==T)
    if (length(idx)>0)
    {
      uniqgenetable[i,]=tmptable[idx[1],]
    }
  }
  return(uniqgenetable)
}

geneposition=generate_geneposition(platforms)
geneposition=filteraliasgene()
searchaliastableforgenepostion=function(geneposition,genealiasposition)
{
  eneposition$gene=as.character(geneposition$gene)
  genealiasposition$gene=as.character(genealiasposition$gene)
  geneposition$chr=as.character(geneposition$chr)
  genealiasposition$chr=as.character(genealiasposition$chr)
  idx=which(is.na(geneposition$chr)==T)
  for (i in idx)
  {
    tmp=which(grepl(paste0("\\<",geneposition$gene[i],"\\>"),genealiasposition$gene)==T)
    if (length(tmp)>0)
    {
      geneposition[i,2:ncol(geneposition)]=genealiasposition[tmp[1],2:ncol(genealiasposition)]
    }
  }
  return(geneposition)
}
geneposition=searchaliastableforgenepostion()
searchgenecardforgeneposition=function(geneposition)
{
  cmd="curl -v -c cookies.txt -b cookies.txt http://www.genecards.org/cgi-bin/carddisp.pl?gene=TP53"
  system(cmd)
  cmd="wget -U Mozilla/5.0 (Windows; U; Windows NT 5.1; de; rv:1.9.2.3) Gecko/20100401 Firefox/3.6.3 http://www.genecards.org/cgi-bin/carddisp.pl?gene=TP53 -O test.txt"
  system(cmd)
}
#geneposition=sortgenetable(geneposition)
write.table(geneposition,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",row.names=F,col.names=T,sep="\t",quote=F)

#instead use refgenes from firehose, hg19_GENCODE_v18_20140127.mat
tmp=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/refgenefiles/hg19_GENECODE_v18_20140127.txt",sep="\t")
geneposition1=tmp[,c(4,1,2,3)]
colnames(geneposition1)=c('gene','chr','start','end')
geneposition1$chr=gsub(23,"X",geneposition1$chr)
geneposition1$chr=gsub(24,"Y",geneposition1$chr)
#geneposition1 has duplicated gene def
idxkeep=rep(T,nrow(geneposition1))
idxdup=duplicated(geneposition1$gene)
idxkeep[idxdup]=F
geneposition1=geneposition1[idxkeep,]
#deals with overlap gene defs
overlaptable=data.frame(matrix(NA,nrow=nrow(geneposition1),ncol=10))
library(GenomicRanges)
gr_geneposition1=GRanges(seqnames = geneposition1$chr,ranges=IRanges(start=geneposition1$start,end=geneposition1$end),gene=geneposition1$gene)
for (i in 1:nrow(geneposition1))
{
  gr_gene=GRanges(seqnames=geneposition1$chr[i],ranges=IRanges(start=geneposition1$start[i],end=geneposition1$end[i]),gene=geneposition1$gene[i])
  olap=subsetByOverlaps(gr_geneposition1,gr_gene)
  if (length(olap)>1)
  {
    count=1
    for (j in 1:length(olap))
    {
      if (mcols(olap)$gene[j] != geneposition1$gene[i])
      {
        overlaptable[i,count]=as.character(mcols(olap)$gene[j])
        count=count+1
      }
    }
  }
}

write.table(geneposition1,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",row.names=F,col.names=T,sep="\t",quote=F)
