#!/usr/bin/env Rscript
#find gene postion.
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
#salloc -t 2-1 -n 100 mpirun -n 1 ./mpi_findgeneposition.R
masterhost=system('uname -n',intern=T)
print(paste0("the masterhost: ",masterhost))
.libPaths(c("/home/xwang234/R/x86_64-unknown-linux-gnu-library/3.0","/app/R/3.2.2/lib/R/library"))
njob=100
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classification_otherdata_stringent.RData")

library("biomaRt")
#mart=useMart("ENSEMBL_MART_ENSEMBL", host="may2009.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #NCBI36, use listDatasets(mart)
mart=useMart("ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #GRCh37,NCBI37
#use biomart
findgeneposition=function(gene)
{
  filter="hgnc_symbol"
  #filter="hgnc_id"
  #filter="entrezgene"
  attributes=c("chromosome_name","start_position","end_position")
  res <- getBM(attributes=attributes, filters=filter, values=gene, mart=mart)
  res1=rep(NA,3)
  if (nrow(res)>0)
  {
    res1=c(res[1,1],res[1,2],res[1,3])
  }
  return(res1)
}

#use knowntable
findgeneposition_knowntable0=function(gene)
{
  allchrs=paste0("chr",c(1:22,"X","Y"))
  res=rep(NA,3)
  idx=which(tolower(knowntable$symbol)==tolower(gene) & knowntable$chr %in% allchrs)
  if (length(idx)>0)
  {
    allchr=as.character(knowntable[idx,"chr"])
    if (length(unique(allchr))>0) #if the gene has multiple position on different chromosomes
    {
      width=rep(0,length(unique(allchr)))
      start=rep(0,length(unique(allchr)))
      end=rep(0,length(unique(allchr)))
      for (j in 1:length(unique(allchr)))
      {
        chr=unique(allchr)[j]
        idx1=idx[which(knowntable$chr[idx]==chr)]
        start[j]=min(knowntable[idx1,"start"],na.rm=TRUE)
        end[j]=max(knowntable[idx1,"end"],na.rm=TRUE)
        width[j]=end[j]-start[j]
      }
      #pick the widest one
      idx2=which.max(width)
      chr=unique(allchr)[idx2]
      start=start[idx2]
      end=end[idx2]
    }else
    {
      chr=unlist(strsplit(chr,"_",fixed=TRUE))
      chr=chr[1]
      chr=gsub("chr","",chr)
      start=min(knowntable[idx,"start"],na.rm=TRUE)
      end=max(knowntable[idx,"end"],na.rm=TRUE)
    }
    res=c(chr,start,end)
  }
  return(res)
}

#only use one item in the table
findgeneposition_knowntable=function(gene)
{
  allchrs=paste0("chr",c(1:22,"X","Y"))
  res=rep(NA,3)
  idx=which(tolower(knowntable$symbol)==tolower(gene) & knowntable$chr %in% allchrs)
  if (length(idx)>0)
  {
    width=knowntable$end[idx]-knowntable$start[idx]
    #pick the widest one
    idx1=which.max(width)
    chr=knowntable$chr[idx][idx1]
    start=knowntable$start[idx][idx1]
    end=knowntable$end[idx][idx1]
    res=c(chr,start,end)
  }
  return(res)
}

mpi_findgeneposition_knowntable=function(genes,outputfile="output.txt",knowntable=knowntable)
{
  
  mpi.bcast.Robj2slave(genes)
  mpi.bcast.Robj2slave(findgeneposition_knowntable)
  mpi.bcast.Robj2slave(knowntable)
  #mpi.remote.exec(system('uname -n',intern=T))
  res1=data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(res1)=c("gene","chr","start","end")
  nrun <- ceiling(length(genes)/1000)
  print("start2")
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):length(genes)
    z=genes[cseq]
    res=mpi.parSapply(X=z,FUN=findgeneposition_knowntable,job.num=njob)
    idx=seq(1,length(res),3)
    res2=data.frame(gene=genes[cseq],chr=res[idx],start=res[idx+1],end=res[idx+2])
    res1=rbind(res1,res2)
    write.table(res1,file=outputfile,col.names=T,row.names=F,sep="\t",quote=F)
  }

#   print("done")
#   mpi.close.Rslaves()
#   mpi.quit()
}

mpi_findgeneposition=function(genes,outputfile="output.txt")
{
#   require(Rmpi)
#   mpi.spawn.Rslaves(needlog = FALSE)
#   .Last <- function()
#   { if (is.loaded("mpi_initialize")){
#     if (mpi.comm.size(1) > 0){
#       print("Please use mpi.close.Rslaves() to close slaves.")
#       mpi.close.Rslaves()
#     }
#     print("Please use mpi.quit() to quit R")
#     .Call("mpi_finalize")
#   }
#   } 
  mpi.bcast.Robj2slave(genes)
  mpi.bcast.Robj2slave(findgeneposition)
  mpi.bcast.Robj2slave(mart)
  mpi.bcast.cmd(library(biomaRt))
  mpi.remote.exec(system('uname -n',intern=T))
  
  res1=data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(res1)=c("gene","chr","start","end")
  nrun <- ceiling(length(genes)/1000)
  print("start2")
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):length(genes)
    z=genes[cseq]
    res=mpi.parSapply(X=z,FUN=findgeneposition,job.num=njob)
    idx=seq(1,length(res),3)
    res2=data.frame(gene=genes[cseq],chr=res[idx],start=res[idx+1],end=res[idx+2])
    res1=rbind(res1,res2)
    write.table(res1,file=outputfile,col.names=T,row.names=F,sep="\t",quote=F)
  }
}
knowntable_ucsc=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGene1.txt",header=TRUE,sep="\t",comment.char="")
knowntable_ucsc$hg19.knownGene.chrom=as.character(knowntable_ucsc$hg19.knownGene.chrom)
knowntable_ucsc$hg19.knownGene.txStart=as.integer(knowntable_ucsc$hg19.knownGene.txStart)
knowntable_ucsc$hg19.knownGene.txEnd=as.integer(knowntable_ucsc$hg19.knownGene.txEnd)
knowntable_ucsc$hg19.kgXref.geneSymbol=as.character(knowntable_ucsc$hg19.kgXref.geneSymbol)
knowntable_ucsc=data.frame(chr=knowntable_ucsc$hg19.knownGene.chrom,start=knowntable_ucsc$hg19.knownGene.txStart,
                      end=knowntable_ucsc$hg19.knownGene.txEnd,symbol=knowntable_ucsc$hg19.kgXref.geneSymbol,stringsAsFactors=F)

knowntable_gistic=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/copynumbergenedefinition_fromgisticmatlabdata.txt",header=T,sep=",")
knowntable_gistic=data.frame(chr=knowntable_gistic$chrn,start=knowntable_gistic$start,end=knowntable_gistic$end,symbol=knowntable_gistic$symb)
knowntable_gistic$chr=as.character(knowntable_gistic$chr)
knowntable_gistic$chr=gsub("23","X",knowntable_gistic$chr)
knowntable_gistic$chr=gsub("24","Y",knowntable_gistic$chr)
if (grepl("chr",knowntable_gistic$chr[1])==F)
{
  knowntable_gistic$chr=paste0("chr",knowntable_gistic$chr)
}

require(Rmpi)
mpi.spawn.Rslaves(needlog = FALSE)
.Last <- function()
{ if (is.loaded("mpi_initialize")){
  if (mpi.comm.size(1) > 0){
    print("Please use mpi.close.Rslaves() to close slaves.")
    mpi.close.Rslaves()
  }
  print("Please use mpi.quit() to quit R")
  .Call("mpi_finalize")
}
} 
platforms=c("copynumber","copynumber_CGH1M","copynumber_1MDUO","mrna","mrna_exon","mrna_u133")
#platforms=c("mrna","mrna_exon","mrna_u133")
for (platform in platforms)
{
  print(platform)
  if (platform=="copynumber")
  {
    data=data_copynumber
  }
  if (platform=="copynumber_CGH1M")
  {
    data=data_copynumber_CGH1M
  }
  if (platform=="copynumber_1MDUO")
  {
    data=data_copynumber_1MDUO
  }
  if (platform=="mrna")
  {
    data=data_mrna
  }
  if (platform=="mrna_exon")
  {
    data=data_mrna_exon
  }
  if (platform=="mrna_u133")
  {
    data=data_mrna_u133
  }
  genes=colnames(data$data1[,14:ncol(data$data1)])
  #modify gene names to run biomart
  genes=gsub("'","",genes,fixed=T)
  genes=gsub("\"","",genes,fixed=T)
#   tmp=sapply(genes,function(i){
#     tmp1=unlist(strsplit(i,"|",fixed=T))
#     res=unname(tmp1[1])
#     return(res)
#   })
#   genes=unname(tmp)
  mpi_findgeneposition_knowntable(genes,
                                  outputfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_ucscknowntable.txt"),
                                  knowntable=knowntable_ucsc)
  mpi_findgeneposition_knowntable(genes,
                                  outputfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_gisticknowntable.txt"),
                                  knowntable=knowntable_gistic)
  #mpi_findgeneposition(genes,outputfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/",platform,"_geneposition_biomart.txt"))
  
}

print("done")
mpi.close.Rslaves()
mpi.quit()
