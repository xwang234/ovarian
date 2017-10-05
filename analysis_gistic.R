
load("../data/firehose.RData")
load('../data/platinum_classificationdata_stringent.RData')
readgisticcna=function(genefile,cnfile=NULL)
{
  gistictable=read.table(genefile,header=F,sep="\t")
  idxNA=is.na(gistictable[1,])
  gistictable=gistictable[,!idxNA]
  res=data.frame(matrix(NA,nrow=ncol(gistictable)-1,ncol=8))
  colnames(res)=c('chr','start','end','cytoband','qvalue','residualq','numgenes','genes')
  for (i in 2:ncol(gistictable))
  {
    tmp=as.character(gistictable[4,i])
    tmp1=unlist(strsplit(tmp,":"))
    tmp2=unlist(strsplit(tmp1[2],"-"))
    res[i-1,1]=tmp1[1]
    res[i-1,2]=tmp2[1]
    res[i-1,3]=tmp2[2]
    res[i-1,4]=as.character(gistictable[1,i])
    res[i-1,5]=as.character(gistictable[2,i])
    res[i-1,6]=as.character(gistictable[3,i])
    genes=as.character(gistictable[5,i])
    n=1
    for (j in 6:nrow(gistictable))
    {
      gene=as.character(gistictable[j,i])
      if (gene!="")
      {
        n=n+1
        genes=paste0(genes,',',gene)
      }
    }
    res[i-1,7]=n
    res[i-1,8]=genes
  }
  #output=paste0(gisticdir,'/',prefix,'.gistic.cna.txt')
  #write.table(res,file=cnafile,row.names=F,col.names=T,sep="\t",quote=F)
  return(res)
}

ampres=readgisticcna("../data/firehose_20160128/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt")
delres=readgisticcna("../data/firehose_20160128/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt")
allgenes=c(ampres$genes,delres$genes)
allgenes=gsub("[","",allgenes,fixed = T)
allgenes=gsub("]","",allgenes,fixed = T)
allgenes1=NULL
for (i in 1:length(allgenes))
{
  allgenes1=c(allgenes1,unlist(strsplit(allgenes[i],",",fixed=T)))
}

allgenes1=unique(allgenes1)

length(allgenes1)
sum(allgenes1 %in% rownames(copynumber))
alldata=rbind(data_copynumber$data1,data_copynumber$data2)
sum(allgenes1 %in% colnames(alldata))
index=which(colnames(alldata) %in% allgenes1)
data_gistic=alldata[,c(11,index)]
pvalues_gistic=rep(0,length(index))
names(pvalues_gistic)=colnames(alldata)[index]
for (i in 2:ncol(data_gistic))
{
  pvalues_gistic[i-1]=calpvalues(data_gistic[,i],data_gistic[,1])
}

pvalues_gistic=pvalues_gistic[order(pvalues_gistic)]
smallfdr=read.table(file="../result/copynumber_result_fdrlessthan0.05.txt",header=T,sep="\t")
sum(allgenes1 %in% smallfdr$gene)
#[1] 17
draw_manhattan(pvalues=pvalues_gistic)

#Aus data:
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/data_aocs_copynumber.RData")

#find genes in copynumber with small qvalues
geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t")
copynumbergenes=rownames(copynumber)
df_copynumbergenes=data.frame(gene=copynumbergenes)
tmp=merge(df_copynumbergenes,geneposition)
tmp=sortgenetable(tmp)
gisticscore=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose_20160128/gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/scores.gistic",header=T,sep="\t")
ampgisticscore=gisticscore[gisticscore$Type=="Amp",]
delgisticscore=gisticscore[gisticscore$Type=="Del",]
idxkeep=rep(F,nrow(tmp))
ampqvalues=delqvalues=rep(NA,nrow(tmp))
library(GenomicRanges)
gr_ampgisticscore=GRanges(seqnames=ampgisticscore$Chromosome,ranges=IRanges(start=ampgisticscore$Start,end=ampgisticscore$End),qvalue=10^(-ampgisticscore$X.log10.q.value.))
gr_delgisticscore=GRanges(seqnames=delgisticscore$Chromosome,ranges=IRanges(start=delgisticscore$Start,end=delgisticscore$End),qvalue=10^(-delgisticscore$X.log10.q.value.))

for (i in 1:nrow(tmp))
{
  gr_gene=GRanges(seqnames=tmp$chr[i],ranges=IRanges(start=tmp$start[i],end=tmp$end[i]))
  olpa=subsetByOverlaps(gr_ampgisticscore,gr_gene)
  if (length(olpa)>0)
  {
    qvalue=min(mcols(olpa)$qvalue)
    ampqvalues[i]=qvalue
    if (qvalue<0.1)
    {
      idxkeep[i]=T
    }
  }
  
  olpd=subsetByOverlaps(gr_delgisticscore,gr_gene)
  if (length(olpd)>0)
  {
    qvalue=min(mcols(olpd)$qvalue)
    delqvalues[i]=qvalue
    if (qvalue<0.1)
    {
      idxkeep[i]=T
    }
  }
}
copynumbergenes_qvalue=data.frame(gene=tmp$gene,chr=tmp$chr,start=tmp$start,end=tmp$end,ampqvalue=ampqvalues,delqvalue=delqvalues)
save(copynumbergenes_qvalue,file="../result/copynumbergeneswithgisticqvalues.RData")

cutoff=0.00001
selgenes=copynumbergenes_qvalue$gene[which(copynumbergenes_qvalue$ampqvalue<cutoff | copynumbergenes_qvalue$delqvalue<cutoff)]
length(selgenes)
sum(smallfdr$gene %in% selgenes)

