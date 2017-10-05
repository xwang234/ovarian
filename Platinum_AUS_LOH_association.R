
rm(list=ls())

#get primary tumor samples
samplefile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/sample.OV-AU.1474317952273.tsv"
allsamples=read.table(file=samplefile,header=T,sep="\t",stringsAsFactors=F)
primarysamples=allsamples[allsamples$specimen_type=="Primary tumour - solid tissue",]
primarysamples=primarysamples[primarysamples$sequencing_strategy=="non-NGS",]
primarysamples=primarysamples[primarysamples$repository=="GEO",]
uniq_primarysamples=unique(primarysamples$submitted_sample_id)

#read resistance def
outcomedeffile="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/nature14410-Supplementary Table 5.xlsx"
library(gdata)
outcomedef=read.xls(outcomedeffile,1,skip=1)
#the last two lines are comments
outcomedef=outcomedef[1:(nrow(outcomedef)-2),]
outcomedef$Chemotherapy.response=as.character(outcomedef$Chemotherapy.response)
outcomedef$Case.ID=as.character(outcomedef$Case.ID)
outcomedef$Case.ID=gsub("_","-",outcomedef$Case.ID,fixed=T)
outcomedef=outcomedef[outcomedef$Sample.type=="7PrimaryTumour",]
sample_platinumstatus=data.frame(matrix(NA,nrow=nrow(outcomedef),ncol=3))
colnames(sample_platinumstatus)=c("patient","sample","status")
sample_platinumstatus$patient=outcomedef$Case.ID
sample_platinumstatus$status=outcomedef$Chemotherapy.response
for (i in 1:nrow(sample_platinumstatus))
{
  idx=which(grepl(sample_platinumstatus$patient[i],uniq_primarysamples)==T)
  if (length(idx)>0)
  {
    sample_platinumstatus$sample[i]=uniq_primarysamples[idx]
  }
}
resist=sample_platinumstatus[,2:3]
names(resist) <- c("sampleID","resistance")

#find samples with pbscbs results
allsamples=list.files(path="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",pattern ="allele.PSCBS.segment.txt" )
allsamples=allsamples[grepl("^AOCS",allsamples)]
allsamples=gsub(".allele.PSCBS.segment.txt","",allsamples)



## first regress platinum resistance on length of LOH ###

LOHseg <- matrix(NA,length(allsamples),4)

for (i in 1:length(allsamples)){
    mysample=allsamples[i]
    output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.PSCBS.segment.txt")
    dat <- read.table(output,header=T)
    LOHseg[i,1] <- sum(dat$dhEnd[!is.na(dat$lohCall) & dat$lohCall]-dat$dhStart[!is.na(dat$lohCall) & dat$lohCall])/10^6
    LOHseg[i,2] <- length(dat$dhEnd[!is.na(dat$lohCall) & dat$lohCall])
    LOHseg[i,3] <- sum(dat$tcnEnd[!is.na(dat$tcnStart) & dat$tcnStart!= "-Inf"] -dat$tcnStart[!is.na(dat$tcnStart) & dat$tcnStart!= "-Inf"])/10^6
    LOHseg[i,4] <- length(dat$tcnEnd[!is.na(dat$tcnStart) & dat$tcnStart!= "-Inf"])
}

loh <- data.frame(cbind(allsamples,LOHseg))
loh[,1] <- as.character(loh[,1])
loh <- loh[order(loh[,1]),]
names(loh) <- c("sampleID","loh_length","n_lohsegs","tcn_length","n_tcnsegs")


# resist <- data_copynumber_tangent_filtered[,1]
# sampleID <- rownames(data_copynumber_tangent_filtered)
# resist <- data.frame(cbind(sampleID,as.character(resist)))
# names(resist) <- c("sampleID","resistance")

resist <- resist[order(resist$sampleID),]
resist <- merge(resist,loh,by="sampleID")
for (i in 3:6)  resist[,i] <- as.numeric(as.character(resist[,i]))
table(resist$resistance)
save(resist,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aus_pscbs_loh_resist.RData")
summary(glm(I(resistance=="sensitive")~loh_length,family=binomial,data=resist))
summary(glm(I(resistance=="sensitive")~n_lohsegs,family=binomial,data=resist))
summary(glm(I(resistance=="sensitive")~tcn_length,family=binomial,data=resist))
summary(glm(I(resistance=="sensitive")~n_tcnsegs,family=binomial,data=resist))



#####################################
## next construct gene-level data ###
#####################################

segsall=c()

for (myfile in allsamples)
{
  tmp=read.table(file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",myfile,".allele.PSCBS.segment.txt"),header=T)
  names(tmp)[1] <- "ID"
  tmp <- tmp[,c(1,2,5,6,7,8,11,12,10,14,34)]
  segsall=rbind(segsall,tmp)
}

#form copynumber data at gene level based on segment result from DNAcopy-------------------

sortgenetable=function(genetable)
{
  genetable$chr=as.character(genetable$chr)
  genetable$chr=gsub("chr","",genetable$chr)
  genetable$chr=gsub("23","X",genetable$chr)
  genetable$chr=gsub("24","Y",genetable$chr)
  if (class(genetable$start[1])=="factor")
  {
    genetable$start=as.numeric(as.character(genetable$start))
  }
  if (class(genetable$start[1])=="character")
  {
    genetable$start=as.numeric(genetable$start)
  }
  
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genetable)))
  for (chr in chrs)
  {
    tmptable=genetable[which(genetable$chr==chr),]
    tmptable=tmptable[order(tmptable$start),]
    res=rbind(res,tmptable)
  }
  return(res)
}



readcnasegs1=function(hg19=T,snp6copynumber,opt=1)
{
  library(CNTools)
  colnames(snp6copynumber)=c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
  dim(snp6copynumber)
  length(unique(snp6copynumber$ID))
  snp6copynumber=snp6copynumber[!is.na(snp6copynumber$seg.mean),]
  cnseg_snp6=CNSeg(snp6copynumber)
  if (opt==1)
  {
    hg19_snp6=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/copynumber_geneposition_biomart.txt",header=T,sep="\t")
    geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
    hg19=merge(hg19_snp6,geneposition,by="gene")
    hg19=hg19[,c(1,5,6,7)]
  }
  if (opt==2) #use the genesets defined by CNTools
  {
    load(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/geneInfohg1819.RData")
    hg19=geneInfo_hg19[,c(4,1,2,3)]
  }
  colnames(hg19)=c("gene","chr","start","end")
  hg19=sortgenetable(hg19)
  
  hg19$chr=gsub("23","X",hg19$chr)
  hg19$chr=gsub("24","Y",hg19$chr)
  idxdup=duplicated(hg19$gene)
  hg19=hg19[! idxdup,]
  #change the column name of chr
  colnames(hg19)[2]="chrom"
  nstart=5
  snp6_max<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "max")
  snp6_max=rs(snp6_max)
  snp6_min<- getRS(cnseg_snp6, by = "gene", imput = FALSE, XY = TRUE, geneMap = hg19, what = "min")
  snp6_min=rs(snp6_min)
  snp6=snp6_max
  
  #pick the extreme one
  for (i in nstart:ncol(snp6))
  {
    pickmin=which(abs(snp6_min[,i])-abs(snp6_max[,i])>0)
    if (length(pickmin)>0) {
      cat(i,"..") 
      snp6[pickmin,i]=snp6_min[pickmin,i]
    }  
  }
  
  rownames(snp6)=snp6$gene
  snp6=snp6[,nstart:ncol(snp6)]
  
  return(snp6)
}

#the data was generated based on mean value of probes
segsall_cn <- segsall[!is.na(segsall$tcnMean),c(1:6)]
segsall_cn$tcnMean <- segsall_cn$tcnMean -2 

segsall_loh <- segsall[segsall$lohCall & !is.na(segsall$dhMean),c(1:2,7:10)]


AUS_gene_copynumber=readcnasegs1(hg19=T,snp6copynumber=segsall_cn)
AUS_gene_copynumber <- AUS_gene_copynumber[,order(names(AUS_gene_copynumber))]
AUS_gene_loh=readcnasegs1(hg19=T,snp6copynumber=segsall_loh)
AUS_gene_loh <- AUS_gene_loh[,order(names(AUS_gene_loh))]


### perform the association test including two-df test ###

outp <- matrix(NA,nrow(AUS_gene_copynumber),6)

for (i in 1:nrow(AUS_gene_copynumber)){
   if ((i %% 1000)==0) cat(i,"..")
   temp <- matrix(NA,ncol(AUS_gene_copynumber),2)
   temp[,1] <- as.vector(names(AUS_gene_copynumber))
   temp[,2] <- as.character(AUS_gene_copynumber[i,])
   
   temp1 <- data.frame(temp)
   names(temp1) <- c("sampleID","tcn")
   
   temp1$tcn <- as.numeric(as.character(temp1$tcn))
   
   if (mean(temp1$tcn!=0)>0.05) {
     resist1 <- merge(resist,temp1,by="sampleID")
     fit1 <- glm(I(resistance=="sensitive")~tcn,family=binomial,data=resist1,x=T,y=T)
     outp[i,1] <- summary(fit1)$coeff[2,4]
   }
   
   if (mean(temp1$tcn>0.5|temp1$tcn<(-0.4))>0.05) {
     resist1 <- merge(resist,temp1,by="sampleID")
     resist1$cn <- 0
     resist1$cn <- ifelse(resist1$tcn>0.5,1,resist1$cn)
     resist1$cn <- ifelse(resist1$tcn<(-0.4),-1,resist1$cn)
     fit3 <- glm(I(resistance=="sensitive")~cn,family=binomial,data=resist1,x=T,y=T)
     outp[i,4] <- summary(fit3)$coeff[2,4]
   }
   
   temp <- matrix(NA,ncol(AUS_gene_loh),2)
   temp[,1] <- as.vector(names(AUS_gene_loh))
   temp[,2] <- as.character(AUS_gene_loh[i,])
   
   temp2 <- data.frame(temp)
   names(temp2) <- c("sampleID","loh")
   temp2$loh <- as.numeric(as.character(temp2$loh))
   
   if (mean(temp2$loh!=0)>0.05) {
     resist1 <- merge(resist,temp2,by="sampleID")
     fit2 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
     outp[i,2] <- summary(fit2)$coeff[2,4]
   }
   
   if (mean(temp2$loh>0.2)>0.05) {
     resist1 <- merge(resist,temp2,by="sampleID")
     resist1$loh <- 1*(resist1$loh>0.2)
     fit4 <- glm(I(resistance=="sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
     outp[i,5] <- summary(fit4)$coeff[2,4]
   }
   
   if (mean(temp2$loh!=0)>0.05  & mean(temp1$tcn!=0)>0.05) {  
     ss1 <- fit1$x*(fit1$y-fit1$fitted)
     ss2 <- fit2$x*(fit2$y-fit2$fitted)
     
     ss <- cbind(ss1,ss2)
     bread <- matrix(0,4,4)
     bread[1:2,1:2] <- t(fit1$x) %*% (fit1$x*(1-fit1$fitted)*fit1$fitted)
     bread[3:4,3:4] <- t(fit2$x) %*% (fit2$x*(1-fit2$fitted)*fit2$fitted)
     
     beef <- t(ss)%*%ss
     
     bcov <- solve(bread) %*% beef %*% solve(bread)
     
     bb <- c(fit1$coef[2],fit2$coef[2])
     outp[i,3] <- 1-pchisq(drop(bb %*% solve(bcov[c(2,4),c(2,4)]) %*% bb),2)
     
   }
   
   if (mean(temp1$tcn>0.5|temp1$tcn<(-0.4))>0.05 & mean(temp2$loh>0.2)>0.05) {  
     ss1 <- fit3$x*(fit3$y-fit3$fitted)
     ss2 <- fit4$x*(fit4$y-fit4$fitted)
   
     ss <- cbind(ss1,ss2)
     bread <- matrix(0,4,4)
     bread[1:2,1:2] <- t(fit3$x) %*% (fit3$x*(1-fit3$fitted)*fit3$fitted)
     bread[3:4,3:4] <- t(fit4$x) %*% (fit4$x*(1-fit4$fitted)*fit4$fitted)
   
     beef <- t(ss)%*%ss
     bcov <- solve(bread) %*% beef %*% solve(bread)
     bb <- c(fit3$coef[2],fit4$coef[2])
     outp[i,6] <- 1-pchisq(drop(bb %*% solve(bcov[c(2,4),c(2,4)]) %*% bb),2)
   }
}

min(outp[!is.na(outp[,1]),1])
length(outp[!is.na(outp[,1]),1])

min(outp[!is.na(outp[,2]),2])
length(outp[!is.na(outp[,2]),2])

min(outp[!is.na(outp[,3]),3])
length(outp[!is.na(outp[,3]),3])

min(outp[!is.na(outp[,4]),4])
length(outp[!is.na(outp[,4]),4])

min(outp[!is.na(outp[,5]),5])
length(outp[!is.na(outp[,5]),5])

min(outp[!is.na(outp[,6]),6])
length(outp[!is.na(outp[,6]),6])

outp <- data.frame(outp)

rownames(outp) <- rownames(AUS_gene_copynumber)
colnames(outp) <- paste0("p",1:6)
outp1=outp

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_pvalues_qvalues.RData")
for (i in 1:6)
{
  print(i)
  idx=which(outp[,6+i]<=0.05)
  if(length(idx)>0)
  {
    print(paste0("number of genes: ",sum(outp1[idx,i]<=0.05,na.rm = T)))
  }
}
outp1=as.data.frame(outp1)
save(outp1,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/AUS_platinum_cnloh.RData")


#draw manhanttan plot:
source("functions.R")
numit=1000
fdrcutoff=0.05
for (i in 1:6)
{
  print(i)
  idx=which(outp[,6+i]<=fdrcutoff)
  print(paste0("number of genes: ",length(idx)))
  if(length(idx)>0)
  {
    print(paste0("number of valid genes: ",sum(outp1[idx,i]<=0.05,na.rm = T)))
    #png(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/manhanttan_tcga_aus",i,".png"),width=1500,height=500)
    postscript(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/manhanttan_tcga_aus",i,".ps"))
    par(mfrow=c(2,1))
    outputprefix=paste0("TCGA",i)
    fdrs=read.table(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"_",numit,"permutations.txt"))
    tcga=draw_manhattan(fdrs=fdrs,fdrthreshold=fdrcutoff,maxy=6,chrs=NULL,keepres=T,logscale=T,main=paste0("TCGA"))
    idx1=which(fdrs$qvalues_permut<=fdrcutoff)
    print(c("chrs (fdr<cutoff) :",unique(tcga$chromosome[idx1])))
    pvalues=outp1[,i]
    names(pvalues)=rownames(outp1)
    idx0=pvalues==0
    if (sum(idx0,na.rm = T)>0) pvalues[idx0]=NA#min(pvalues[!idx0],na.rm=T)
    aus=draw_manhattan(pvalues = pvalues,keepres = T, main=paste0("AUS"))
    idx1=which(outp1[idx,i]<=0.05)
    validgenes=rownames(outp1[idx[idx1],])
    idx1=match(validgenes,rownames(aus))
    #View(aus[idx1,])
    ymax=max(-log10(aus$pvalue),na.rm=T)
    if (length(idx1)>0) #confirmed in AUS
    {
      #idxs=idx[idx1]
      points(aus$allpos[idx1],rep(ymax*1.05,length(idx1)),pch=8,col="red",cex=1.3)
    }
    dev.off()
  }
  
}

