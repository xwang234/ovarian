#load functions
source("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
#draw manhanttan plots
numit=1000
for (i in 1:6)
{
  outputprefix=paste0("TCGA",i)
  fdrs=read.table(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"_",numit,"permutations.txt"))
  postscript(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/manhplot_fdr",i,".ps"),width=12,height=6)
  draw_manhattan(fdrs=fdrs,fdrthreshold=0.05,maxy=6,chrs=NULL,keepres=F,logscale=T)
  dev.off()
}

#save q values
#load the observed p-values
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_platinum_cnloh.RData")
numit=1000
for (i in 1:6)
{
  outputprefix=paste0("TCGA",i)
  fdrs=read.table(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/qvalues_",outputprefix,"_",numit,"permutations.txt"))
  outp=cbind(outp,fdrs[,2])
}
colnames(outp)=c(paste0("p",1:6),paste0("q",1:6))
save(outp,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA_pvalues_qvalues.RData")
