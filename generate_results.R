#!/usr/bin/env Rscript
library(glmnet)
library(gap) #Manhattan plot
library(GenomicRanges) #find overlaps
library(ROCR)
library(pROC)
library(gdata) #read.xls
#include functions
source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
#load data:
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent.RData")
#include data_copynumber and data_mrna.
#data1 for training
dim(data_copynumber$data1)
#the clinical variables
colnames(data_copynumber$data1)[1:13]
#data2 for testing
dim(data_copynumber$data2)

#plot ROC-------------------------------------------------------------------------------------
#for copynumber:
#20 iterations of crossvalidation used
class_copynumber=classify_platinum(data=data_copynumber,platform="copynumber",selgenes=mutation_selgenes,includeclinical=T)
#the result is a list. x,y:training data; x1,y1:testing data; pfit_train:fitted output for training data.
names(class_copynumber)
#[1] "coeff"      "x"          "y"          "x1"         "y1"         "pfit_train" "pfit_test" 
#variables selected:
class_copynumber$coeff

#for mrna
class_mrna=classify_platinum(data=data_mrna,platform="mrna",selgenes=NULL,includeclinical=T)

#consider genes with correlation between CNA and mRNA
cor_mrna_copynumber=compute_corr_mrna_cna()
highcorgenes=names(cor_mrna_copynumber)[order(cor_mrna_copynumber,decreasing = T)]
class_copynumberhighcor=classify_platinum(data=data_copynumber,platform="copynumber",selgenes=highcorgenes[1:1000],includeclinical=T)
class_mrnahighcor=classify_platinum(data=data_mrna,platform="mrna",selgenes=highcorgenes[1:1000],includeclinical=T)

#work on mutation data, it is sparse
data=data_mutation$data1
mutation_sum=colSums(data[,14:ncol(data)])
names(mutation_sum)=colnames(data)[14:ncol(data)]
sum(mutation_sum>=10)
#[1] 22
sum(mutation_sum>2)
#[1] 568
mutation_selgenes=names(mutation_sum)[mutation_sum>2]
hr30genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/HR30genes.txt",sep="\t",header=T)
hr30genes=as.character(hr30genes$Gene)
hr30genes[c(17,23)]
#"MSH2 (+EPCAM)"  "PTEN (+KILLIN)"
hr30genes[17]="MSH2"
hr30genes[23]="PTEN"
hr30genes=c(hr30genes,"EPCAM","KILLIN")
class_mutation=classify_platinum(data=data_mutation,platform="mutation",selgenes=mutation_selgenes,includeclinical=T)
class_mutation=classify_platinum(data=data_mutation,platform="mutation",selgenes=hr30genes,includeclinical=T)
#compute p-values-------------------------------------------------------------------------------
#set runnpermutation=T will run permutations and save the pvalues in 100 permuations in a text file, it takes 2-2.5 hours.
#permutaion results are saved in /fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_permutation1.txt and
# /fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_permutation1.txt

Sys.time()
pvalues_copynumber=compute_pvalues(data=data_copynumber,platform="copynumber",runpermutation=F)
str(pvalues_copynumber)
Sys.time()
pvalues_mrna=compute_pvalues(data=data_mrna,platform="mrna",runpermutation=F)
Sys.time()
#draw histogram of pvalues:
hist(pvalues_copynumber,main="copynumber",xlab="p-value",probability=T)
hist(pvalues_mrna,main="mrna",xlab="p-value",probability=T,ylim=c(0,2))

#compute fdr------------------------------------------------------------------------------------
#read pvalues from permutations
#it takes 10-20 minutes for each computation
#fdr results are also saved in /fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber1.txt and
# /fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_mrna1.txt

Sys.time()
#read the permutation pvalues genreated in the previous step
pvalues_copynumber_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_permutation1.txt",header=F)
fdrs_copynumber=sapply_compute_qvalue_permutation(pvalues=pvalues_copynumber,pvalues_permutation=pvalues_copynumber_permutation,
                                                  numit=100,platform="copynumber")
Sys.time()
#or the result out of 1000 permutations
pvalues_copynumber_1000permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_copynumber_1000permutation.txt",header=T)
fdrs1000p_copynumber=sapply_compute_qvalue_permutation(pvalues=pvalues_copynumber,pvalues_permutation=pvalues_copynumber_1000permutation,
                                                  numit=1000,platform="copynumber")
Sys.time()

sum(fdrs_copynumber<0.05)
Sys.time()

pvalues_mrna_permutation=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/pvalues_mrna_permutation1.txt",header=F)
fdrs_mrna=sapply_compute_qvalue_permutation(pvalues=pvalues_mrna,pvalues_permutation=pvalues_mrna_permutation,
                                                  platform="mrna")
sum(fdrs_mrna<0.05)
Sys.time()

pvalues_mutation=compute_pvalues(data=data_mutation,platform="mutation",runpermutation=F)
#draw Manhattan plot-----------------------------------------------------------------------------
#the fdr threshold to keep genes:
fdrthreshold=0.05
# if the fdr was not generated in previous step, read the previous saved result to save time
if (! exists("fdrs1000p_copynumber"))
{
  fdrs1000p_copynumber=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_1000permutation.txt")
}
draw_manhattan(fdrs=fdrs1000p_copynumber,fdrthreshold=fdrthreshold)

#for mrna
if (! exists("fdrs1000p_mrna"))
{
  fdrs1000p_mrna=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_mrna_1000permutation.txt")
}

draw_manhattan(fdrs=fdrs1000p_mrna,maxy=5,fdrthreshold=fdrthreshold)

#for copynumber_CGH1M
fdrs1000p_copynumber_CGH1M=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_CGH1M_1000permutation.txt")
draw_manhattan(fdrs=fdrs1000p_copynumber_CGH1M,maxy=5,fdrthreshold=fdrthreshold)

#generate tables for genes detected in copynumber with small fdr---------------------------------------------
smallfdr1000p_copynumber=smallfdrtable(fdrs=fdrs1000p_copynumber,data=data_copynumber,threshold=fdrthreshold)
smallfdr1000p_mrna=smallfdrtable(fdrs=fdrs1000p_mrna,data=data_mrna,threshold=fdrthreshold)
#write.table(smallfdr1000p_copynumber,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/copynumber_result_fdrlessthan0.05.txt",row.names = F,col.names = T,sep="\t",quote=F)
