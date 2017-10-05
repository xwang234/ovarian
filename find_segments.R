#! /usr/bin/env Rscript

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData")

library(GenomicRanges)
geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
#only keep genes in copynumber
geneposition=geneposition[geneposition$gene %in% rownames(copynumber),]
gr_geneposition = GRanges(seqnames=geneposition$chr,ranges=IRanges(start=geneposition$start,end=geneposition$end),gene=geneposition$gene)

segfile="../data/firehose_20160128/gdac.broadinstitute.org_OV.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"

segs=read.table(file=segfile,sep="\t",quote="",header=TRUE,stringsAsFactors=FALSE)
allsamples=unique(segs$Sample)
#just for check
test=segs[segs$Sample=="TCGA-04-1335-11A-01D-0428-01",]
test1=segs[segs$Sample=="TCGA-04-1335-01A-01D-0428-01",]
plot(x=c(1,247813706),y=c(0,0),ylim=c(-2,2))
for (i in 1:11)
{
  lines(c(test$Start[i],test$End[i]),rep(test$Segment_Mean[i],2),col="blue")
}
for (i in 1:10)
{
  lines(c(test1$Start[i],test1$End[i]),rep(test1$Segment_Mean[i],2),col="red")
}
legend("topleft",legend=c("normal","tumor"),col=c("blue","red"),lty=c(1,1))

sampletypes=sapply(allsamples,function(asample){
  res=unlist(strsplit(asample,"-",fixed=T))
  res=res[4]
})

samplenames=sapply(allsamples,function(asample){
  res=unlist(strsplit(asample,"-",fixed=T))
  res=paste0(res[1:3],collapse="-")
})

tumoridx=grepl("01",sampletypes)==T
#gr_segs=GRanges(seqnames=segs$Chromosome,ranges=IRanges(start=segs$Start,end=segs$End),Segment_Mean=segs$Segment_Mean)
# gr_segs1=sort(gr_segs)
# gr_segs1=unique(gr_segs1)
# gr_segs2=NULL

tumorsamples=allsamples[tumoridx]
tumorsamplenames=samplenames[tumoridx]
length(tumorsamplenames)
length(unique(tumorsamplenames))
if (length(tumorsamplenames) != length(unique(tumorsamplenames))) warning("there are duplicated tumor samples")

#only keep tumor segs
segs=segs[segs$Sample %in% tumorsamples,]
segs$Chromosome=gsub(23,"X",segs$Chromosome)
segs$Chromosome=gsub(24,"Y",segs$Chromosome)
chrs=unique(segs$Chromosome)
if ("Y" %in% chrs) chrs=chrs[-which(chrs=="Y")]
gr_segs=GRanges(seqnames=segs$Chromosome,ranges=IRanges(start=segs$Start,end=segs$End),cn=segs$Segment_Mean,samplename=segs$Sample )
gr_segs=sort(gr_segs)
upthreshold=0.2
lowthreshold = -0.2
gr_finalsegs=GRanges()
for (chr in chrs)
{
  cat(chr,'...')
  tmp=segs[segs$Chromosome==chr,]
  breaks=unique(c(tmp$Start,tmp$End))
  breaks=breaks[order(breaks)]
  #start idx
  idx=1:(length(breaks)-1)
  gr_breaks=GRanges(seqnames=chr,ranges=IRanges(start=breaks[idx],end=breaks[idx+1]),freq=rep(0,length(idx)))
  for (i in 1: length(gr_breaks))
  {
    olp=subsetByOverlaps(gr_segs,gr_breaks[i,],minoverlap=10)
    #the ones have cna greater than thresholds
    idx1=mcols(olp)$cn >= upthreshold | mcols(olp)$cn <= lowthreshold
    olp=olp[idx1,]
    mcols(gr_breaks)$freq[i]=length(unique(mcols(olp)$samplename))/length(unique(segs$Sample))
  }
  gr_finalsegs=append(gr_finalsegs,gr_breaks)
}

save(gr_finalsegs,file="test.RData")
gr_finalsegs=gr_finalsegs[mcols(gr_finalsegs)$freq>0,]

idx_freq=mcols(gr_finalsegs)$freq >0.05
width_segs=width(gr_finalsegs)

#plot segments
gr_segs1=gr_segs[mcols(gr_segs)$cn>0.2 | mcols(gr_segs)$cn < -0.2,]
n=10938
plot(c(min(start(gr_segs1)[1:n]),max(end(gr_segs1)[1:n])),c(0,0),ylim=c(0,1.1))
for (i in 1:n)
{
  if (mcols(gr_segs1)$cn[i]>0)
  {
    mycol="red"
  }else
  {
    mycol="blue"
  }
  lines(c(start(gr_segs1)[i],end(gr_segs1)[i]),c(i/n,i/n),col=mycol)
}


finalsegs=data.frame(matrix(NA,nrow=length(gr_finalsegs),ncol=5))
colnames(finalsegs)=c("chr","start","end","freq","genes")
finalsegs$chr=seqnames(gr_finalsegs)
finalsegs$start=start(gr_finalsegs)
finalsegs$end=end(gr_finalsegs)
finalsegs$freq=mcols(gr_finalsegs)$freq
for (i in 1:length(gr_finalsegs))
{
  olp=subsetByOverlaps(gr_geneposition,gr_finalsegs[i,],minoverlap=10)
  if (length(olp)>0)
  {
    mygenes=mcols(olp)$gene
    if (length(mygenes)>10)
    {
      mygenes=paste0(mygenes[1:10],collapse="__")
    }else
    {
      mygenes=paste0(mygenes,collapse="__")
    }
    finalsegs$genes[i]=mygenes
  }
}

finalsegs1=finalsegs[!is.na(finalsegs$genes),]
finalsegs1$chr=as.character(finalsegs1$chr)
rownames(finalsegs1)=1:nrow(finalsegs1)
length(unique(finalsegs1$genes))
finalsegs2=data.frame(matrix(NA,nrow=nrow(finalsegs1),ncol=4))
colnames(finalsegs2)=c("chr","start","end","genes")
count=1
oldgenes=finalsegs1$genes[1]
oldgenes=unlist(strsplit(oldgenes,"__",fixed=T))
i=2
startline=1
endline=1
while(i<=nrow(finalsegs1))
{
  genes=finalsegs1$genes[i]
  genes=unlist(strsplit(genes,"__",fixed=T))
  #new segments
  if (sum(genes %in% oldgenes)==0)
  {
    endline=i-1
    finalsegs2$chr[count]=as.character(finalsegs1$chr[i-1])
    finalsegs2$start[count]=min(finalsegs1$start[startline:endline])
    finalsegs2$end[count]=max(finalsegs1$end[startline:endline])
    startline=i
    alloldgenes=paste0(oldgenes,collapse="__")
    finalsegs2$genes[count]=alloldgenes
    oldgenes=genes
    count=count+1
  }else
  {
    oldgenes=c(oldgenes,genes)
    oldgenes=unique(oldgenes)
  }
  i=i+1
}
#to add the last row
endline=i-1
finalsegs2$chr[count]=as.character(finalsegs1$chr[i-1])
finalsegs2$start[count]=min(finalsegs1$start[startline:endline])
finalsegs2$end[count]=max(finalsegs1$end[startline:endline])
alloldgenes=paste0(oldgenes,collapse="__")
finalsegs2$genes[count]=alloldgenes
finalsegs2=finalsegs2[!is.na(finalsegs2$chr),]

