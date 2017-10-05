#! /usr/bin/env Rscript

source(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
pvalues_mrna=read.table(file="../result/pvalues_mrna_platinum.txt",header=T,quote="",sep="\t")
draw_manhattan(pvalues=pvalues_mrna,main="TCGA_mRNA")

load("../data/tcga_data.RData")
fdrs_copynumber=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_copynumber_1000permutation.txt",header=T)
draw_manhattan(fdrs=fdrs_copynumber,main="TCGA_copynumber")
draw_manhattan(pvalues=pvalues_copynumber_allele,main="TCGA_copynumber_DH")
draw_manhattan(pvalues=pvalues_2df,main="TCGA_2df")

load("../data/aocs_data.RData")
draw_manhattan(pvalues=pvalues_aocs_mean_copynumber,main="AOCS_copynumber")
draw_manhattan(pvalues=pvalues_aocs_mean_copynumber_allele,main="AOCS_copynumber_DH")
draw_manhattan(pvalues=pvalues_aocs_mean_copynumber_rmcnv_2df,main="AOCS_copynumber_2df")
