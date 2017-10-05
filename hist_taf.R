#! /usr/bin/env Rscript

folder="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/bcgsc_protectedmutation/bcgsc.ca_OV.Multicenter_mutation_calling_MC3_Cont.Level_2.1.0.0/snvfrqfiles"

files=list.files(folder)
for (i in 1:length(files))
{
  tcgabarcode=unlist(strsplit(files[i],".",fixed=T))
  myfile=paste0(folder,"/",files[i])
  mytable=read.table(myfile,header=T,sep="\t")
  taf=mytable$tumor_alt/(mytable$tumor_ref+mytable$tumor_alt)
  test=which(taf>0.8)
  hist(taf,breaks=20,main=tcgabarcode[3],xlab="Tumor allelic fraction")
}