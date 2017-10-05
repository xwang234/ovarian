#!/usr/bin/env Rscript
njob=100
library("biomaRt")
#mart=useMart("ENSEMBL_MART_ENSEMBL", host="may2009.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #NCBI36, use listDatasets(mart)
mart=useMart("ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #GRCh37,NCBI37
#use biomart

#alias gene is in the format of aaa|bbb|xxx
findgeneposition=function(gene)
{
  res1=rep(NA,3)
  allchrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","23","24")
  filter="hgnc_symbol"
  #filter="hgnc_id"
  #filter="entrezgene"
  attributes=c("chromosome_name","start_position","end_position")
  mygenes=unlist(strsplit(gene,"|",fixed=T))
  for (i in 1:length(mygenes))
  {
    res <- getBM(attributes=attributes, filters=filter, values=mygenes[i], mart=mart)
    if (nrow(res)>0)
    {
      if (sum(as.character(res$chromosome_name) %in% allchrs)>0)
      {
        idx=which(as.character(res$chromosome_name) %in% allchrs)
        res1=c(res[idx[1],1],res[idx[1],2],res[idx[1],3])
        break
      }
    }
  }
  return(res1)
}

mpi_findgeneposition=function(genes,outputfile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/aliasgeneposition_biomart.txt")
{

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
  return(res1)
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

aliastable=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/genes_multiplesymbols.txt",sep="\t",fill=T,quote="",header=T,stringsAsFactors=F )
multiplegenes=paste(aliastable$symbol,aliastable$alias_symbol,sep="|")
multiplegenes=gsub("\\|$","",multiplegenes,perl=T)
multiplegenes=gsub("\"","",multiplegenes,perl=T)
multiplegenes=gsub("\'","",multiplegenes,perl=T)

aliasgenepos=mpi_findgeneposition(genes=multiplegenes,outputfile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/aliasgeneposition_biomart.txt")
mpi.close.Rslaves()
mpi.quit()
