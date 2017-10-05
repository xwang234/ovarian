
omni_1.0=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni2.5-8v1_C_Gene_Annotation.txt",header=T,sep="\t",stringsAsFactors=F)
omni_1.0=omni_1.0[omni_1.0$Gene.s.!="",]
omni_1.1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/HumanOmni25M-8v1-1_B.annotated.txt",header=T,sep="\t",fill=T,stringsAsFactors=F)
omni_1.1=omni_1.1[omni_1.1$Gene.s.!="",]
