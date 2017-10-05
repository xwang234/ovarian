
rm(list=ls())

filelist=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/filelist.txt",header=T,sep="\t",stringsAsFactors=F)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
allsamples=rownames(data_copynumber_filtered)
idxtumor=sapply(allsamples,function(mysample){
  res=which(grepl(paste0(mysample,"-01"),filelist$sampleid)==T & grepl("byallele",filelist$filename)==T)
  #for multiple files
  if (mysample=="TCGA-23-1023") res=res[2]
  if (mysample=="TCGA-29-2414") res=res[1]
  if (mysample=="TCGA-13-1817") res=res[2]
  if (mysample=="TCGA-13-1819") res=res[1]
  return(res)
})

idxnormal=sapply(allsamples,function(mysample){
  res=which(grepl(paste0(mysample,"-1"),filelist$sampleid)==T & grepl("byallele",filelist$filename)==T)
  return(res)  
})

idxnormal1=sapply(1:length(allsamples),function(i){
  if (length(idxnormal[[i]])==0)
  {
    res=NA
  }else
  {
    res=idxnormal[[i]]
  }
  names(res)=names(idxnormal)[i]
  return(res)
})
sum(is.na(idxnormal1))
#[1] 12 samples without normal allele files
idxnormal=idxnormal1

snp_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp_anno <- snp_anno[,1:4]

## first regress platinum resistance on length of LOH ###

LOHseg <- matrix(NA,length(idxnormal1),4)

for (i in 1:length(idxnormal1)){
mysample=allsamples[i]
if (!i %in% c(122,359,363)) {
    output=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".allele.PSCBS.segment.txt")
    dat <- read.table(output,header=T)
    LOHseg[i,1] <- sum(dat$dhEnd[!is.na(dat$lohCall) & dat$lohCall]-dat$dhStart[!is.na(dat$lohCall) & dat$lohCall])/10^6
    LOHseg[i,2] <- length(dat$dhEnd[!is.na(dat$lohCall) & dat$lohCall])
    LOHseg[i,3] <- sum(dat$tcnEnd[!is.na(dat$tcnStart) & dat$tcnStart!= "-Inf"] -dat$tcnStart[!is.na(dat$tcnStart) & dat$tcnStart!= "-Inf"])/10^6
    LOHseg[i,4] <- length(dat$tcnEnd[!is.na(dat$tcnStart) & dat$tcnStart!= "-Inf"])
    
}    
}

loh <- data.frame(cbind(allsamples,LOHseg))
loh[,1] <- as.character(loh[,1])
loh <- loh[order(loh[,1]),]
names(loh) <- c("sampleID","loh_length","n_lohsegs","tcn_length","n_tcnsegs")

save(loh,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_loh.RData")


resist <- data_copynumber_tangent_filtered[,1]
sampleID <- rownames(data_copynumber_tangent_filtered)
resist <- data.frame(cbind(sampleID,as.character(resist)))
names(resist) <- c("sampleID","resistance")

resist <- resist[order(resist$sampleID),]
resist <- merge(resist,loh,by="sampleID")
for (i in 3:6)  resist[,i] <- as.numeric(as.character(resist[,i]))


summary(glm(I(resistance=="Sensitive")~loh_length,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~n_lohsegs,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~tcn_length,family=binomial,data=resist))
summary(glm(I(resistance=="Sensitive")~n_tcnsegs,family=binomial,data=resist))


library(beeswarm)
postscript("~/SCHARPHOME/ScharpFile/Article/TCGA_ovarian/LOH_platinum.ps",width=8,height=6)
beeswarm(resist$loh_length/3000 ~resist$resistance,pch=19,cex=1.2,xlab="Platinum response", ylab="% genome has LOH")
bxplot(resist$loh_length/3000 ~resist$resistance,probs=0.5,col=2,add=T)
dev.off()

postscript("~/SCHARPHOME/ScharpFile/Article/TCGA_ovarian/LOHsegs_platinum.ps",width=8,height=6)
beeswarm(resist$n_lohsegs ~resist$resistance,pch=19,cex=1.2,xlab="Platinum response", ylab="number of LOH segments",ylim=c(0,300))
bxplot(resist$n_lohsegs ~resist$resistance,probs=0.5,col=2,add=T)
dev.off()



#####################################
## next construct gene-level data ###
#####################################

segsall=c()

for (myfile in allsamples[-c(122,359,363)])
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


TCGA_gene_copynumber=readcnasegs1(hg19=T,snp6copynumber=segsall_cn)
TCGA_gene_copynumber <- TCGA_gene_copynumber[,order(names(TCGA_gene_copynumber))]
TCGA_gene_loh=readcnasegs1(hg19=T,snp6copynumber=segsall_loh)
TCGA_gene_loh <- TCGA_gene_loh[,order(names(TCGA_gene_loh))]

save(TCGA_gene_loh,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_gene_loh.RData")


### perform the association test including two-df test ###

outp <- matrix(NA,nrow(TCGA_gene_copynumber),6)

for (i in 1:nrow(TCGA_gene_copynumber)){
   if ((i %% 1000)==0) cat(i,"..")
   temp <- matrix(NA,ncol(TCGA_gene_copynumber),2)
   temp[,1] <- as.vector(names(TCGA_gene_copynumber))
   temp[,2] <- as.character(TCGA_gene_copynumber[i,])
   
   temp1 <- data.frame(temp)
   names(temp1) <- c("sampleID","tcn")
   
   temp1$tcn <- as.numeric(as.character(temp1$tcn))
   
   if (mean(temp1$tcn!=0)>0.05) {
     resist1 <- merge(resist,temp1,by="sampleID")
     fit1 <- glm(I(resistance=="Sensitive")~tcn,family=binomial,data=resist1,x=T,y=T)
     outp[i,1] <- summary(fit1)$coeff[2,4]
   }
   
   if (mean(temp1$tcn>0.5|temp1$tcn<(-0.4))>0.05) {
     resist1 <- merge(resist,temp1,by="sampleID")
     resist1$cn <- 0
     resist1$cn <- ifelse(resist1$tcn>0.5,1,resist1$cn)
     resist1$cn <- ifelse(resist1$tcn<(-0.4),-1,resist1$cn)
     fit3 <- glm(I(resistance=="Sensitive")~cn,family=binomial,data=resist1,x=T,y=T)
     outp[i,4] <- summary(fit3)$coeff[2,4]
   }
   
   temp <- matrix(NA,ncol(TCGA_gene_loh),2)
   temp[,1] <- as.vector(names(TCGA_gene_loh))
   temp[,2] <- as.character(TCGA_gene_loh[i,])
   
   temp2 <- data.frame(temp)
   names(temp2) <- c("sampleID","loh")
   temp2$loh <- as.numeric(as.character(temp2$loh))
   
   if (mean(temp2$loh!=0)>0.05) {
     resist1 <- merge(resist,temp2,by="sampleID")
     fit2 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
     outp[i,2] <- summary(fit2)$coeff[2,4]
   }
   
   if (mean(temp2$loh>0.2)>0.05) {
     resist1 <- merge(resist,temp2,by="sampleID")
     resist1$loh <- 1*(resist1$loh>0.2)
     fit4 <- glm(I(resistance=="Sensitive")~loh,family=binomial,data=resist1,x=T,y=T)
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

rownames(outp) <- rownames(TCGA_gene_copynumber)


save(outp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_platinum_cnloh.RData")

###################################################################
#### study the joint prediction from LOH segs and gene-level LOH ##
###################################################################

source("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_platinum_cnloh.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_gene_loh.RData")

hg19_snp6=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genepositions/copynumber_geneposition_biomart.txt",header=T,sep="\t")
geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition.txt",header=T,sep="\t",stringsAsFactors=F)
hg19=merge(hg19_snp6,geneposition,by="gene")
hg19=hg19[,c(1,5,6,7)]
colnames(hg19)=c("gene","chr","start","end")
hg19=sortgenetable(hg19)

hg19$chr=gsub("23","X",hg19$chr)
hg19$chr=gsub("24","Y",hg19$chr)
idxdup=duplicated(hg19$gene)
hg19=hg19[! idxdup,]
#change the column name of chr
colnames(hg19)[2]="chrom"

################################################################
### find out gene names that top in LOH association analysis ###
################################################################
topgenes <- NULL
for (i in c(3,4,5,6,10,15,17)){
  kk <- which.min(outp[!is.na(outp[,5]) & hg19$chrom==i,5])
  topgenes<- c(topgenes,as.character(hg19$gene[!is.na(outp[,5]) & hg19$chrom==i][kk]))
}  

topgene_loh <- TCGA_gene_loh[row.names(TCGA_gene_loh) %in% topgenes,]
topgene_loh <- t(topgene_loh)
topgene_loh <- cbind(names(TCGA_gene_loh),topgene_loh)
topgene_loh <- data.frame(topgene_loh)

names(topgene_loh)[1] <- "sampleID"

for (i in 2:8) topgene_loh[,i] <- as.numeric(as.character(topgene_loh[,i]))

resist1 <- merge(resist,topgene_loh,by="sampleID")

summary(glm(I(resistance=="Sensitive")~LRRN1 + PDLIM5+ARHGEF28+MAP3K4+DIP2C+GPR176+OR1E1,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~PDLIM5,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~LRRN1,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~ARHGEF28,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~MAP3K4,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~DIP2C,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~GPR176,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~OR1E1,family=binomial,data=resist1))

cor(resist1[,c(3,7:13)])

summary(glm(I(resistance=="Sensitive")~loh_length+LRRN1 + PDLIM5+ARHGEF28+MAP3K4+DIP2C+GPR176,family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~loh_length+I(LRRN1>0.2) + I(PDLIM5>0.2)+I(ARHGEF28>0.2)+I(MAP3K4>0.2)+I(DIP2C>0.2)+I(GPR176>0.2)+I(OR1E1>0.2),family=binomial,data=resist1))
summary(glm(I(resistance=="Sensitive")~I(LRRN1>0.2) + I(PDLIM5>0.2)+I(ARHGEF28>0.2)+I(MAP3K4>0.2)+I(DIP2C>0.2)+I(GPR176>0.2)+I(OR1E1>0.2),family=binomial,data=resist1))


#plot.roc(I(resist$resistance=="Sensitive"), resist$loh_length,xlim=c(1,0))

fit1 <- glm(I(resistance=="Sensitive")~loh_length+n_lohsegs,family=binomial,data=resist1)
pred1 <- fit1$fitted

roc1 <- roc(I(resistance=="Sensitive")~pred1,data=resist1)
roc2 <- roc(I(resistance=="Sensitive")~I(LRRN1>0.2) + I(PDLIM5>0.2)+I(ARHGEF28>0.2)+I(MAP3K4>0.2)+I(DIP2C>0.2)+I(GPR176>0.2)+I(OR1E1>0.2),data=resist1)


fit2 <- glm(I(resistance=="Sensitive")~I(LRRN1>0.2) + I(PDLIM5>0.2)+I(ARHGEF28>0.2)+I(MAP3K4>0.2)+I(DIP2C>0.2)+I(GPR176>0.2)+I(OR1E1>0.2),family=binomial,data=resist1)
pred2 <- fit2$fitted

roc2 <- roc(I(resistance=="Sensitive")~pred2,data=resist1)

roc.test(roc1, roc2)

##DeLong's test for two correlated ROC curves
##data:  roc1 and roc2
##Z = -4.543, p-value = 5.546e-06
##alternative hypothesis: true difference in AUC is not equal to 0
##sample estimates:
##AUC of roc1 AUC of roc2 
##  0.6341999   0.7490458 


fit3 <- glm(I(resistance=="Sensitive")~loh_length+I(LRRN1>0.2) + I(PDLIM5>0.2)+I(ARHGEF28>0.2)+I(MAP3K4>0.2)+I(DIP2C>0.2)+I(GPR176>0.2),family=binomial,data=resist1)
pred3 <- fit3$fitted


cutpt <- c(pred1,pred2,pred3)
cutpt <- cutpt[order(cutpt)]
TP1 <- rep(0,length(cutpt))
FP1 <- rep(0,length(cutpt))

for (i in 1:length(cutpt)) {
  TP1[i] <- mean(pred1[resist1$resistance=="Sensitive"]>cutpt[i]) 
  FP1[i] <- mean(pred1[resist1$resistance!="Sensitive"]>cutpt[i]) 
}

#cutpt <- pred2[order(pred2)]
TP2 <- rep(0,length(cutpt))
FP2 <- rep(0,length(cutpt))
for (i in 1:length(cutpt)) {
  TP2[i] <- mean(pred2[resist1$resistance=="Sensitive"]>cutpt[i]) 
  FP2[i] <- mean(pred2[resist1$resistance!="Sensitive"]>cutpt[i]) 
}

#cutpt <- pred3[order(pred3)]
TP3 <- rep(0,length(cutpt))
FP3 <- rep(0,length(cutpt))
for (i in 1:length(cutpt)) {
  TP3[i] <- mean(pred3[resist1$resistance=="Sensitive"]>cutpt[i]) 
  FP3[i] <- mean(pred3[resist1$resistance!="Sensitive"]>cutpt[i]) 
}


postscript("~/SCHARPHOME/ScharpFile/Article/TCGA_ovarian/ROC_platinum.ps",width=6,height=6)
plot(FP1,TP1,type="l",xlim=c(0,1),ylim=c(0,1),lwd=3,xlab="1-specificity",ylab="Sensitivity")
lines(FP2,TP2,col=2,lwd=3)
#lines(FP3,TP3,col=3,lwd=3)
segments(0,0,1,1)
legend(0.4,0.2,c("HRD: AUC=0.64","LOH associations: AUC=0.75"),col=1:2,lwd=3,bty="n")
dev.off()

######################
### now try glmnet ###
######################

lohmat <- t(TCGA_gene_loh)

del <- rep(0,ncol(lohmat))
for (l in 1:ncol(lohmat)){
   if (sum(is.na(lohmat[,l]))==368 | sum(lohmat[,l]>0.2) <5) del[l] <- 1
}


gdat <- as.matrix(lohmat)
for (l in 1:ncol(gdat))  gdat[,l] <- I(gdat[,l]>0.2)

outcome <- I(resist1$resistance=="Sensitive")
fit <- glmnet(gdat,outcome,family="binomial",nlambda=100)
cvfit <- cv.glmnet(gdat,outcome,family="binomial",nfolds=10,nlambda=100)
plot(cvfit)
coeff <- as.matrix(coef(fit,s=cvfit$lambda.min))[,1]
selected <- which(as.matrix(coef(fit,s=cvfit$lambda.min))[,1]!=0)
coeff[selected]


  

