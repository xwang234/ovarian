#!/usr/bin/env Rscript

#chemo-resistance defination
clinical_classification=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/clinical_vairables.txt",header=T,sep="\t",stringsAsFactors = F)
load("../data/TCGA_methylation27.RData")
sample1=colnames(methylation)[grepl("-01$",colnames(methylation))]
sample2=colnames(methylation)[grepl("-02$",colnames(methylation))]
sample3=colnames(methylation)[grepl("-11$",colnames(methylation))] #4 normal samples
#samplename
sample11=sapply(sample1,function(x){
  tmp=unlist(strsplit(x,"-",fixed = T))
  paste0(tmp[1:3],collapse = "-")
})
sample22=sapply(sample2,function(x){
  tmp=unlist(strsplit(x,"-",fixed = T))
  paste0(tmp[1:3],collapse = "-")
})

length(intersect(sample11,sample22))
#[1] 16
length(sample22)
#[1] 18
#there are two samples only with -02
idxsample2=which(!sample22 %in% sample11)
methysamples=c(sample1,sample2[idxsample2])
length(methysamples)
#[1] 584 tumor samples having methylation data
methylation=cbind(methylation[,1:4],methylation[,colnames(methylation) %in% methysamples])
colnames(methylation)[5:ncol(methylation)]=sapply(colnames(methylation)[5:ncol(methylation)],function(x){
  tmp=unlist(strsplit(x,"-",fixed = T))
  paste0(tmp[1:3],collapse = "-")
})

#check methylation samples having def
alldefsamples=clinical_classification$sample[clinical_classification$platinumclass %in% c("Sensitive","Resistant")]
methydefsamples=colnames(methylation)[colnames(methylation) %in% alldefsamples]
length(methydefsamples)
#[1] 388 samples have methylation and def

#check copynumber data
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
copynumberdefsamples=rownames(data_copynumber_filtered)
length(intersect(copynumberdefsamples,methydefsamples))
#[1] 369 samples have copynumber,methylation,def

#save data
resistall=clinical_classification[clinical_classification$platinumclass %in%c("Sensitive","Resistant"),c("sample","platinumclass","new_data")]
methylation=cbind(methylation[,1:4],methylation[,colnames(methylation) %in% methydefsamples])
for (i in 5:ncol(methylation))
{
  methylation[,i]=as.numeric(as.character(methylation[,i]))
}
save(methylation,resistall,file="../data/TCGA_methylation.RData")

#try bumphunter----------------------
load("../data/TCGA_methylation.RData")
#add age and residual disease, there is no smoking data
clinical_classification=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/clinical_vairables.txt",header=T,sep="\t",stringsAsFactors = F)
sum(is.na(methylation$Chromosome))
#remove probes having many NAs
methylation1 = subset(methylation,subset = (rowSums(is.na(methylation[,5:ncol(methylation)])) <ncol(methylation)-4))
rownames(methylation1)=methylation1$`Composite Element REF`
methylationdata=methylation1[,5:ncol(methylation1)]
NARows=which(rowSums(is.na(methylationdata))>0)
numNA=rowSums(is.na(methylationdata[NARows,]))
quantile(numNA)
#replace these NAs with mean
for (i in NARows)
{
  idx=which(is.na(methylationdata[i,]))
  methylationdata[i,idx]=mean(as.numeric(methylationdata[i,]),na.rm=T)
}


#PCA
pc=prcomp(t(methylationdata),scale=T,center=T)
#the order of pca is correct
sum(rownames(pc$x)==colnames(methylationdata))


clinicaldata=clinical_classification[,c("sample","age","residual_disease_largest_nodule")]
sum(colnames(methylationdata) %in% clinicaldata$sample)
idx=match(colnames(methylationdata),clinicaldata$sample)
clinicaldata=clinicaldata[idx,]

#add resistance
sum(clinicaldata$sample %in% resistall$sample)
#[1] 388
clinicaldata=merge(clinicaldata,resistall,by="sample")
idx=match(colnames(methylationdata),clinicaldata$sample)
clinicaldata=clinicaldata[idx,]
sum(clinicaldata$sample==rownames(pc$x))
clinicaldata=cbind(clinicaldata,pc$x)
sum(clinicaldata$sample==colnames(methylationdata))
clinicaldata=clinicaldata[,-c(which(colnames(clinicaldata)=="sample"),which(colnames(clinicaldata)=="new_data"))]
colnames(clinicaldata)[which(colnames(clinicaldata)=="residual_disease_largest_nodule")]="residualdisease"
clinicaldata=cbind(clinicaldata[,3],clinicaldata[,-3])
colnames(clinicaldata)[1]="resistance"
clinicaldata$resistance=factor(clinicaldata$resistance,levels = c("Sensitive","Resistant"))



out <- matrix(NA,nrow(methylationdata),5)
colnames(out)=c("p_pc10","coef_pc10","p","coef","var")
for (l in 1:nrow(methylationdata)) {
  yy <- as.numeric(log(methylationdata[l,]/(1-methylationdata[l,]),base=2))
  fit1 <- glm(yy~resistance+age+I(residualdisease=="No Macroscopic disease")+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=clinicaldata)
  out[l,1] <- summary(fit1)$coef[2,4]
  out[l,2] <- summary(fit1)$coef[2,1]
  fit2 <-  glm(yy~resistance+age+I(residualdisease=="No Macroscopic disease"),data=clinicaldata)
  out[l,3] <- summary(fit2)$coef[2,4]
  out[l,4] <- summary(fit2)$coef[2,1]
  out[l,5] <- var(yy)
  if (l %% 1000==0) cat(l,"..")
}
out_with_residualdisease=out
nprobes=nrow(methylationdata)
plot(-log((1:nprobes)/nprobes,base=10),-log(out[order(out[,1]),1],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="with 10pc adjusted")
abline(0,1)
plot(-log((1:nprobes)/nprobes,base=10),-log(out[order(out[,3]),3],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="without pc adjusted")
abline(0,1)
out_with_residualdisease=out

out <- matrix(NA,nrow(methylationdata),5)
colnames(out)=c("p_pc10","coef_pc10","p","coef","var")
for (l in 1:nrow(methylationdata)) {
  yy <- as.numeric(log(methylationdata[l,]/(1-methylationdata[l,]),base=2))
  fit1 <- glm(yy~resistance+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=clinicaldata)
  out[l,1] <- summary(fit1)$coef[2,4]
  out[l,2] <- summary(fit1)$coef[2,1]
  fit2 <-  glm(yy~resistance+age,data=clinicaldata)
  out[l,3] <- summary(fit2)$coef[2,4]
  out[l,4] <- summary(fit2)$coef[2,1]
  out[l,5] <- var(yy)
  if (l %% 1000==0) cat(l,"..")
}

nprobes=nrow(methylationdata)
plot(-log((1:nprobes)/nprobes,base=10),-log(out[order(out[,1]),1],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="with 10pc adjusted")
abline(0,1)
plot(-log((1:nprobes)/nprobes,base=10),-log(out[order(out[,3]),3],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="without pc adjusted")
abline(0,1)

out_20pc=out

out <- matrix(NA,nrow(methylationdata),5)
colnames(out)=c("p_pc10","coef_pc10","p","coef","var")
for (l in 1:nrow(methylationdata)) {
  yy <- as.numeric(log(methylationdata[l,]/(1-methylationdata[l,]),base=2))
  fit1 <- glm(yy~resistance+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=clinicaldata)
  out[l,1] <- summary(fit1)$coef[2,4]
  out[l,2] <- summary(fit1)$coef[2,1]
  fit2 <-  glm(yy~resistance+age,data=clinicaldata)
  out[l,3] <- summary(fit2)$coef[2,4]
  out[l,4] <- summary(fit2)$coef[2,1]
  out[l,5] <- var(yy)
  if (l %% 1000==0) cat(l,"..")
}
out_10pc=out

nprobes=nrow(methylationdata)
plot(-log((1:nprobes)/nprobes,base=10),-log(out[order(out[,1]),1],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="with 20pc adjusted")
abline(0,1)
plot(-log((1:nprobes)/nprobes,base=10),-log(out[order(out[,3]),3],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="without pc adjusted")
abline(0,1)

#draw manhanttan plot
source("functions.R")
pvalues=data.frame(chr=methylation1$Chromosome,start=methylation1$Genomic_Coordinate,pvalue=out_20pc[,1])
pvalues=pvalues[!is.na(pvalues$chr),]
pvalues$chr=as.character(pvalues$chr)
pvalues=sortgenetable(pvalues)
postscript("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_methylation_pvalues.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
draw_manhattan(pvalues = pvalues,maxy=7,keepres=F,logscale=T,ylab="observed p-value (-log base 10)")
dev.off()
save(methylation1,methylationdata,clinicaldata,out_10pc,out_20pc,out_with_residualdisease_10pc,out_with_residualdisease_20pc,file="../result/tmp/methylation_result.RData")
#out_xxx and methylation1 have the same row index


#bumphunter------------------------
library(bumphunter)
methylation1$Genomic_Coordinate=as.integer(as.character(methylation1$Genomic_Coordinate))
methylation1$Gene_Symbol=as.character(methylation1$Gene_Symbol)
designMat <- model.matrix(~resistance+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=clinicaldata)

# library(doParallel)
# registerDoParallel(cores = 1)
methylation1$Chromosome=droplevels(methylation1$Chromosome)
#the data need to be sorted
methylation2=cbind(chr=methylation1$Chromosome,start=methylation1$Genomic_Coordinate,gene=methylation1$Gene_Symbol,methylationdata)
methylation2=sortgenetable(methylation2)
methylation2$gene=as.character(methylation2$gene)
set.seed(1000)
dmrs <- bumphunter(as.matrix(methylation2[,4:ncol(methylation2)]),designMat,chr=methylation2$chr,pos=methylation2$start,cutoff=0.03,nullMethod="bootstrap",B=1000,type="beta")
head(dmrs$table,n=2)
# chr    start      end       value      area cluster indexStart indexEnd     L clusterL      p.value  fwer  p.valueArea fwerArea
# 680   8 22960385 22961088 -0.06676431 0.3338215   16337       1195    22929 21735        5 5.749770e-06 0.004 0.0002673643    0.172
# 507  17 42083193 42083195 -0.10860010 0.2172002    7615       7961    22704 14744        2 1.437443e-05 0.010 0.0011427668    0.565
sum(methylation2$chr==pvalues$chr & methylation2$start==pvalues$start)
out1 <- matrix(0,0,7)
colnames(out1)=c("region","chr","pos","gene","zscore","pvalue","bumphunter")
for (i in 1:2)
{
  idx=which(methylation2$chr==dmrs$table$chr[i] & methylation2$start>=dmrs$table$start[i] & methylation2$start<=dmrs$table$end[i]) #the ones found by bumphunter
  
  if (length(idx)>0)
  {
    idx1=idx[1]
    p=pvalues$pvalue[idx1]
    while(p<=0.05)
    {
      idx1=idx1-1
      p=pvalues$pvalue[idx1]
    }
    idxstart=idx1+1
    idx2=idx[length(idx)]
    p=pvalues$pvalue[idx2]
    while(p<=0.05)
    {
      idx2=idx2+1
      p=pvalues$pvalue[idx2]
    }
    idxend=idx2-1
    idxall=idxstart:idxend
    idxall=unique(c(idx,idxall))
    idxall=idxall[order(idxall)]
    for (j in idxall)
    {
      yy <- as.numeric(log(methylation2[j,4:ncol(methylation2)]/(1-methylation2[j,4:ncol(methylation2)]),base=2))
      fit1 <- glm(yy~resistance+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=clinicaldata)
      inbump=0
      if (j %in% idx)
      {
        inbump=1
      }
      out1=rbind(out1,c(i,dmrs$table$chr[i],methylation2$start[j],methylation2$gene[j],summary(fit1)$coef[2,3],summary(fit1)$coef[2,4],inbump))
    }
  }
}

for (i in 1:2)
{
  idx=which(out1[,1]==i)
  plot(out1[idx,3],out1[idx,5],xlab="position",ylab="z-score",ylim=c(-4,0),type="n",cex.axis=0.75)
  pchs=rep(1,length(idx))
  pchs[out1[idx,7]==1]=4
  cols=rep(1,length(idx))
  cols[out1[idx,7]==1]=2
  points(out1[idx,3],out1[idx,5],pch=pchs,col=cols,cex=1.5)

  genes=unique(out1[idx,4])
  genes1=NULL
  for (j in 1:length(genes))
  {
    tmp=unlist(strsplit(genes[j],";"))
    genes1=c(genes1,tmp)
    genes1=unique(genes1)
  }
  genes1=paste0(genes1,collapse = ",")
  title(paste0("Chromosome ",out1[idx[1],2],": ",genes1,", CpG Island"))
  if (length(unique(pchs))>1)
  {
    legend("topleft",legend=c("not selected by bumphunter","selected by bumphunter"),pch=c(1,4),col=c(1,2))
  }
}
library(qvalue)
qvalues= qvalue(out_20pc[,1])$qvalues
min(qvalues)
#[1] 0.0820771
which(qvalues<=0.1)
#[1] 13573
methylation1[13573,1:4]
# Composite Element REF Chromosome Genomic_Coordinate Gene_Symbol
# cg14947494            cg14947494          6          166796284      BRP44L
out_20pc[13573,1]
#3.539994e-06
save(methylation1,methylationdata,clinicaldata,out_10pc,out_20pc,out_with_residualdisease_10pc,out_with_residualdisease_20pc,
     methylation2,dmrs,file="../result/tmp/methylation_result.RData")
head(dmrs$table[,c(1,2,3,9,11,12)])

hist(out_20pc[,1],xlab="p-value",main="",cex.axis=1.2,cex.lab=1.2)
idx=which(qvalues<=0.1)
qvalues[idx]
#[1] 0.0820771
methylation1$Gene_Symbol[idx]
#[1] "BRP44L"
methylation1$Chromosome[idx]
#[1] "6"
out_20pc[idx,1]
#3.539994e-06 
methylation1$Genomic_Coordinate[idx]
#[1] 166796284

idx=which(qvalues<=0.2 & qvalues>0.1)
methylation1$Gene_Symbol[idx]
#[1] "TMEM55A" "SFT2D1"  "DYNLT1" 
methylation1$Chromosome[idx]
#[1] "8" "6" "6"
out_20pc[idx,1]
#[1] 2.780233e-05 3.142342e-05 1.942106e-05
methylation1$Genomic_Coordinate[idx]
#[1]  92053433 166755991 159066182

idx=which(grepl("BRCA",methylation1$Gene_Symbol))
out_20pc[idx,1]
#[1] 0.7985151 0.7822813 0.2758139 0.7572201 0.6089451 0.0987584 0.1669172 0.9999579 0.8763073 0.4189306 0.7908438

idx=which(grepl("BRCA1",methylation1$Gene_Symbol))
plotzsocres=function(idx,plotgene=1,miny=-5,maxy=3,plotfig=1)
{
  out2 <- matrix(0,0,5)
  colnames(out2)=c("chr","pos","gene","zscore","pvalue")
  for (i in idx)
  {
    yy <- as.numeric(log(methylationdata[i,]/(1-methylationdata[i,]),base=2))
    fit1 <- glm(yy~resistance+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=clinicaldata)
    out2=rbind(out2,c(methylation1$Chromosome[i],methylation1$Genomic_Coordinate[i],methylation1$Gene_Symbol[i],summary(fit1)$coef[2,3],summary(fit1)$coef[2,4]))
  }
  if (plotfig==1)
  {
    plot(out2[,2],out2[,4],xlab="position",ylab="z-score",ylim=c(miny,maxy),type="n",cex.axis=0.75,cex.lab=1.2)
    points(out2[,2],out2[,4],cex=1.5)
    segments(par("usr")[1],1.95,par("usr")[2],1.95,lty=2)
    segments(par("usr")[1],-1.95,par("usr")[2],-1.95,lty=2)
    
    genes=unique(out2[,3])
    genes1=NULL
    for (j in 1:length(genes))
    {
      tmp=unlist(strsplit(genes[j],";"))
      genes1=c(genes1,tmp)
      genes1=unique(genes1)
    }
    genes1=paste0(genes1,collapse = ",")
    if (plotgene==1)
    {
      title(paste0("Chromosome ",out2[1,1],": ",genes1,", CpG Island"))
    }else
    {
      title(paste0("Chromosome ",out2[1,1],": "," CpG Island"))
    }
  }
  
  return(out2)
}
idx=which(grepl("BRCA1",methylation1$Gene_Symbol))
plotzsocres(idx,miny=-3)

idx=which(grepl("BRCA2",methylation1$Gene_Symbol))
plotzsocres(idx,miny=-3)

idx=which(qvalues<0.2 & methylation1$Chromosome=="6")
plotzsocres(idx)

idx1=which(qvalues<0.2 & methylation1$Chromosome=="6")
pos1=min(methylation1$Genomic_Coordinate[idx1])
pos2=max(methylation1$Genomic_Coordinate[idx1])
idx=which(methylation1$Genomic_Coordinate>=pos1 & methylation1$Genomic_Coordinate<=pos2 & methylation1$Chromosome=="6")
plotzsocres(idx,plotgene = 0)

epigenetic_silencing=read.table("../data/TCGA_methylation_epigenetic_silencing168.txt",header=T,sep=" ",stringsAsFactors = F)
sum(epigenetic_silencing$ProbeID %in% methylation1$`Composite Element REF`)
sum(epigenetic_silencing$ProbeID %in% methylation$`Composite Element REF`)
idxmissing=which(epigenetic_silencing$ProbeID %in% methylation$`Composite Element REF` & !epigenetic_silencing$ProbeID %in% methylation1$`Composite Element REF` )
epigenetic_silencing$Gene_Symbol[idxmissing]
#[1] "KLK11"  "CLIP3"  "KLK10"  "GSTM2"  "GIMAP2" "PEG3" 
sum(length(unique(epigenetic_silencing$Gene_Symbol)))
probid=epigenetic_silencing$ProbeID
out2=NULL
for (i in 1:length(probid))
{
  idx=which(methylation1$`Composite Element REF`==probid[i])
  if (length(idx)>0)
  {
    tmp=plotzsocres(idx,plotfig = 0)
    tmp=cbind(probid[i],epigenetic_silencing$Gene_Symbol[i],tmp)
    out2=rbind(out2,tmp)
  }
}  

out3=NULL
for (i in 1:nrow(epigenetic_silencing))
{
  idx=which(grepl(paste0("^",epigenetic_silencing$Gene_Symbol[i],"$"),methylation1$Gene_Symbol) |
                  grepl(paste0("^",epigenetic_silencing$Gene_Symbol[i],";"),methylation1$Gene_Symbol) |
                  grepl(paste0(";",epigenetic_silencing$Gene_Symbol[i],";"),methylation1$Gene_Symbol) |
                  grepl(paste0(";",epigenetic_silencing$Gene_Symbol[i],"$"),methylation1$Gene_Symbol))
  if (length(idx)>0)
  {
    tmp=plotzsocres(idx,plotfig = 0)
    tmp=cbind(probid[i],epigenetic_silencing$Gene_Symbol[i],tmp)
    out3=rbind(out3,tmp)
  }
}
idx2=which(out2[,7]<=0.05)
idx3=which(out3[,7]<=0.05)
sum(out2[idx2,1] %in% out3[idx3,1])
out3[idx3,1][!out3[idx3,1] %in% out2[idx2,1]]
es_res=out3[which(out3[,7]<=0.05),]
es_res=as.data.frame(es_res)
colnames(es_res)[1]="ProbeID"
colnames(es_res)[2]="genename"
colnames(es_res)[5]="probe_target_gene"
length(unique(es_res$genename))
#[1] 12
write.table(es_res,file="../result/tmp/epigenetic_silencing_smallpvalue.txt",sep="\t",quote=F,col.names = T,row.names = F)
out3=as.data.frame(out3)
out3[,7]=as.numeric(as.character(out3[,7]))
plot(-log((1:nrow(out3))/nrow(out3),base=10),-log(out3[order(out3[,7]),7],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="")
abline(0,1)

#read mrna data
methylation1$`Composite Element REF`=as.character(methylation1$`Composite Element REF`)
load("../../Ovarian/TCGAprocesseddata.RData")
gene="FOLR1"

plot_gene_methy_geneexp=function(gene,methylationdata,geneexpdata,resistall)
{
  colnames(geneexpdata)=gsub("A$","",colnames(geneexpdata))
  colnames(geneexpdata)=gsub("B$","",colnames(geneexpdata))
  methylationdata1=methylationdata
  colnames(methylationdata1)=paste0(colnames(methylationdata),"-01")
  commsamples=intersect(colnames(methylationdata1),colnames(geneexpdata))
  print(paste0("number of samples: ",length(commsamples)))
  idx=match(commsamples,colnames(methylationdata1))
  methylationdata1=methylationdata1[,idx]
  idx=match(commsamples,colnames(geneexpdata))
  geneexpdata=geneexpdata[,idx]
  idx=match(commsamples,paste0(resistall$sample,"-01"))
  resist=resistall[idx,]
  #sum(rownames(methylationdata1)==methylation1$`Composite Element REF`)
  idx1=which(methylation1$Gene_Symbol==gene)
  idx2=which(rownames(geneexpdata)==gene)
  colors=rep("red",length(commsamples))
  colors[resist$platinumclass=="Resistant"]="blue"
  if (length(idx1)>0 & length(idx2)==1)
  {
    for (i in idx1[1])
    {
      print(t.test(as.numeric(methylationdata1[i,resist$platinumclass=="Sensitive"]),as.numeric(methylationdata1[i,resist$platinumclass=="Resistant"])))
      x=as.numeric(methylationdata1[i,])
      y=as.numeric(geneexpdata[idx2,])
      fit=lm(y~x)
      summary(fit)
      plot(x,y,col=colors,xlab="Methylation",ylab="Gene expression",cex.lab=1.2,cex.axis=1.2)
      legend("topright",legend=c("Sensitive","Resistant"),col=c("red","blue"),pch=1)
      boxplot(x~factor(resist$platinumclass),ylab="methylation",cex.lab=1.2,cex.axis=1.2)
      boxplot(y~factor(resist$platinumclass),ylab="gene expression",cex.lab=1.2,cex.axis=1.2)
      print(t.test(x[resist$platinumclass=="Sensitive"],x[resist$platinumclass=="Resistant"]))
      print(t.test(y[resist$platinumclass=="Sensitive"],y[resist$platinumclass=="Resistant"],alternative = "greater"))
    }
  }
}
geneexpdata=mrna_gene_hgu133a
plot_gene_methy_geneexp(gene,methylationdata,geneexpdata,resistall)


geneexpdata=mrna_gene_agilent
plot_gene_methy_geneexp(gene,methylationdata,geneexpdata,resistall)
# for (i in 1:ncol(mrna_exon))
# {
#   mrna_exon[,i]=as.numeric(as.character(mrna_exon[,i]))
# }

geneexpdata=mrna_exon
plot_gene_methy_geneexp(gene,methylationdata,geneexpdata,resistall)
