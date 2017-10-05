#!/usr/bin/env Rscript

#salloc -t 0-5 -n 90 mpirun -n 1 R --interactive
library(data.table)
library(ASCAT)

# # helper function to split the genome into parts
# split_genome = function(SNPpos) {
#   # look for gaps of more than 5Mb (arbitrary treshold to account for big centremeres or other gaps) and chromosome borders
#   bigHoles = which(diff(SNPpos[,2])>=5000000)+1
#   chrBorders = which(SNPpos[1:(dim(SNPpos)[1]-1),1]!=SNPpos[2:(dim(SNPpos)[1]),1])+1
#   holes = unique(sort(c(bigHoles,chrBorders)))
#   # find which segments are too small
#   #joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
#   # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
#   #while (1 %in% joincandidates) {
#   # holes=holes[-1]
#   # joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
#   #}
#   #while ((length(holes)+1) %in% joincandidates) {
#   # holes=holes[-length(holes)]
#   # joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
#   #}
#   #while(length(joincandidates)!=0) {
#   # the while loop is because after joining, segments may still be too small..
#   #startseg = c(1,holes)
#   #endseg = c(holes-1,dim(SNPpos)[1])
#   # for each segment that is too short, see if it has the same chromosome as the segments before and after
#   # the next always works because neither the first or the last segment is in joincandidates now
#   #previoussamechr = SNPpos[endseg[joincandidates-1],1]==SNPpos[startseg[joincandidates],1]
#   #nextsamechr = SNPpos[endseg[joincandidates],1]==SNPpos[startseg[joincandidates+1],1]
#   #distanceprevious = SNPpos[startseg[joincandidates],2]-SNPpos[endseg[joincandidates-1],2]
#   #distancenext = SNPpos[startseg[joincandidates+1],2]-SNPpos[endseg[joincandidates],2]
#   # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
#   #joins = ifelse(previoussamechr&nextsamechr,
#   # ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
#   # ifelse(nextsamechr, joincandidates, joincandidates-1))
#   #holes=holes[-joins]
#   #joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
#   #}
#   # if two neighboring segments are selected, this may make bigger segments then absolutely necessary, but I'm sure this is no problem.
#   startseg = c(1,holes)
#   endseg = c(holes-1,dim(SNPpos)[1])
#   chr=list()
#   for (i in 1:length(startseg)) {
#     chr[[i]]=startseg[i]:endseg[i]
#   }
#   return(chr)
# }
# 
# #only load a single sample
# ascat.loadData1=function(Tumor_LogR, Tumor_BAF, Germline_LogR = NULL, 
#           Germline_BAF = NULL, chrs = c(1:22, "X", "Y"), gender = NULL, 
#           sexchromosomes = c("X", "Y")) 
# {
#   # print.noquote("Reading Tumor LogR data...")
#   # #Tumor_LogR <- read.table(Tumor_LogR_file, header = T, row.names = 1, comment.char = "", sep = "\t", check.names = F)
#   # print.noquote("Reading Tumor BAF data...")
#   #Tumor_BAF <- read.table(Tumor_BAF_file, header = T, row.names = 1, comment.char = "", sep = "\t", check.names = F)
#   Tumor_LogR[Tumor_LogR == -Inf] = NA
#   Tumor_LogR[Tumor_LogR == Inf] = NA
#   
#   if (!is.null(Germline_LogR))
#   {
#     Germline_LogR[Germline_LogR == -Inf] = NA
#     Germline_LogR[Germline_LogR == Inf] = NA
#   }
#   # Germline_LogR = NULL
#   # Germline_BAF = NULL
#   # if (!is.null(Germline_LogR_file)) {
#   #   print.noquote("Reading Germline LogR data...")
#   #   Germline_LogR <- read.table(Germline_LogR_file, header = T, 
#   #                               row.names = 1, comment.char = "", sep = "\t", check.names = F)
#   #   print.noquote("Reading Germline BAF data...")
#   #   Germline_BAF <- read.table(Germline_BAF_file, header = T, 
#   #                              row.names = 1, comment.char = "", sep = "\t", check.names = F)
#   #   Germline_LogR[Germline_LogR == -Inf] = NA
#   #   Germline_LogR[Germline_LogR == Inf] = NA
#   # }
#   print.noquote("Registering SNP locations...")
#   SNPpos <- Tumor_LogR[, 1:2]
#   SNPpos = SNPpos[SNPpos[, 1] %in% chrs, ]
#   chrs = intersect(chrs, unique(SNPpos[, 1]))
#   Tumor_LogR = Tumor_LogR[rownames(SNPpos), c(-1, -2), drop = F]
#   Tumor_BAF = Tumor_BAF[rownames(SNPpos), c(-1, -2), drop = F]
#   for (cc in 1:dim(Tumor_LogR)[2]) {
#     Tumor_LogR[, cc] = as.numeric(as.vector(Tumor_LogR[, 
#                                                        cc]))
#     Tumor_BAF[, cc] = as.numeric(as.vector(Tumor_BAF[, cc]))
#   }
#   if (is.null(Germline_LogR)) {
#     Germline_LogR = Germline_LogR[rownames(SNPpos), c(-1, 
#                                                       -2), drop = F]
#     Germline_BAF = Germline_BAF[rownames(SNPpos), c(-1, -2), 
#                                 drop = F]
#     for (cc in 1:dim(Germline_LogR)[2]) {
#       Germline_LogR[, cc] = as.numeric(as.vector(Germline_LogR[, 
#                                                                cc]))
#       Germline_BAF[, cc] = as.numeric(as.vector(Germline_BAF[, 
#                                                              cc]))
#     }
#     
#   }
#   # if (!is.null(Germline_LogR_file)) {
#   #   Germline_LogR = Germline_LogR[rownames(SNPpos), c(-1, 
#   #                                                     -2), drop = F]
#   #   Germline_BAF = Germline_BAF[rownames(SNPpos), c(-1, -2), 
#   #                               drop = F]
#   #   for (cc in 1:dim(Germline_LogR)[2]) {
#   #     Germline_LogR[, cc] = as.numeric(as.vector(Germline_LogR[, 
#   #                                                              cc]))
#   #     Germline_BAF[, cc] = as.numeric(as.vector(Germline_BAF[, 
#   #                                                            cc]))
#   #   }
#   # }
#   last = 0
#   ch = list()
#   SNPorder = vector(length = dim(SNPpos)[1])
#   for (i in 1:length(chrs)) {
#     chrke = SNPpos[SNPpos[, 1] == chrs[i], ]
#     chrpos = chrke[, 2]
#     names(chrpos) = rownames(chrke)
#     chrpos = sort(chrpos)
#     ch[[i]] = (last + 1):(last + length(chrpos))
#     SNPorder[ch[[i]]] = names(chrpos)
#     last = last + length(chrpos)
#   }
#   SNPpos = SNPpos[SNPorder, ]
#   Tumor_LogR = Tumor_LogR[SNPorder, , drop = F]
#   Tumor_BAF = Tumor_BAF[SNPorder, , drop = F]
#   # if (!is.null(Germline_LogR_file)) {
#   #   Germline_LogR = Germline_LogR[SNPorder, , drop = F]
#   #   Germline_BAF = Germline_BAF[SNPorder, , drop = F]
#   # }
#   if (!is.null(Germline_LogR)) {
#     Germline_LogR = Germline_LogR[SNPorder, , drop = F]
#     Germline_BAF = Germline_BAF[SNPorder, , drop = F]
#   }
#   print.noquote("Splitting genome in distinct chunks...")
#   chr = split_genome(SNPpos)
#   if (is.null(gender)) {
#     gender = rep("XX", dim(Tumor_LogR)[2])
#   }
#   return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
#               Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
#               Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
#               SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, samples = colnames(Tumor_LogR), 
#               gender = gender, sexchromosomes = sexchromosomes, failedarrays = NULL))
# }

#load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat.RData")
setwd("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/")
# mpi_ascat=function(jobn)
# {
#   res=FALSE
#   if (!exists("tumorbaf")) load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tcga_ascat.RData")
#   library(ASCAT)
#   if (jobn+2<=ncol(tumorbaf))
#   {
#     Tumor_LogR=tumorlogr[,c(1:2,2+jobn)]
#     #rownames(Tumor_LogR)=tumorlogr[,1]
#     Tumor_BAF = tumorbaf[,c(1:2,2+jobn)]
#     #rownames(Tumor_BAF)=tumorbaf[,1]
#     Germline_LogR = normallogr[,c(1:2,2+jobn)]
#     #rownames(Germline_LogR)=normallogr[,1]
#     Germline_BAF = normalbaf[,c(1:2,2+jobn)]
#     #rm(tumorlogr,tumorbaf,normallogr,normalbaf)
#     #rownames(Germline_BAF)=normalbaf[,1]
#     # mpi.bcast.Robj2slave(Tumor_BAF)
#     # mpi.bcast.Robj2slave(Tumor_LogR)
#     # mpi.bcast.Robj2slave(Germline_BAF)
#     # mpi.bcast.Robj2slave(Germline_LogR)
#     ascat.bc3 = ascat.loadData1(Tumor_LogR, Tumor_BAF, Germline_LogR, Germline_BAF)
#     ascat.bc3 = ascat.GCcorrect(ascat.bc3,"/fh/fast/dai_j/CancerGenomics/Tools/database/other/GC_AffySNP6_102015.txt")
#     ascat.plotRawData(ascat.bc3)
#     Sys.time()
#     ascat.bc3 = ascat.aspcf(ascat.bc3)
#     Sys.time()
#     ascat.plotSegmentedData(ascat.bc3)
#     #Sys.time()
#     ascat.output = ascat.runAscat(ascat.bc3)
#     if (!is.null(ascat.output$segments))
#     {
#       res=TRUE
#       mysample=ascat.bc3$samples
#       write.table(ascat.output$segments,file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.segment"),
#                   row.names = F,col.names = T,sep="\t",quote=F)
#       filecon=file(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.res"),"w")
#       writeLines(paste0("psi ",ascat.output$psi),filecon)
#       writeLines(paste0("aberrantcellfraction ",ascat.output$aberrantcellfraction),filecon)
#       writeLines(paste0("ploidy ",ascat.output$ploidy),filecon)
#       writeLines(paste0("goodnessOfFit ",ascat.output$goodnessOfFit),filecon)
#       close(filecon)
#     }
#   }
#   return(res)
# }

#samples with normals
allsamples=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_samplelist.txt",stringsAsFactors = F)
mpi_ascat1=function(jobn)
{
  res=FALSE
  library(ASCAT)
  mysample=allsamples[jobn,1]
  tumorlogrfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_input/",mysample,"_tumorlogr.txt")
  tumorbaffile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_input/",mysample,"_tumorbaf.txt")
  normallogrfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_input/",mysample,"_normallogr.txt")
  normalbaffile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/ascat_input/",mysample,"_normalbaf.txt")
  ascat.bc3 = ascat.loadData(tumorlogrfile, tumorbaffile, normallogrfile, normalbaffile)
  ascat.bc3 = ascat.GCcorrect(ascat.bc3,"/fh/fast/dai_j/CancerGenomics/Tools/database/other/GC_AffySNP6_102015.txt")
  ascat.plotRawData(ascat.bc3)
  Sys.time()
  ascat.bc3 = ascat.aspcf(ascat.bc3)
  Sys.time()
  ascat.plotSegmentedData(ascat.bc3)
  #Sys.time()
  ascat.output = ascat.runAscat(ascat.bc3)
  if (!is.null(ascat.output$segments))
  {
    res=TRUE
    mysample=ascat.bc3$samples
    write.table(ascat.output$segments,file=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.segment"),
                row.names = F,col.names = T,sep="\t",quote=F)
    filecon=file(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/",mysample,".ascat.res"),"w")
    writeLines(paste0("psi ",ascat.output$psi),filecon)
    writeLines(paste0("aberrantcellfraction ",ascat.output$aberrantcellfraction),filecon)
    writeLines(paste0("ploidy ",ascat.output$ploidy),filecon)
    writeLines(paste0("goodnessOfFit ",ascat.output$goodnessOfFit),filecon)
    close(filecon)
  }
  return(res)
}

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
# mpi.bcast.Robj2slave(split_genome)
# mpi.bcast.Robj2slave(ascat.loadData1)
mpi.bcast.Robj2slave(allsamples)
#res=mpi.parSapply(X=1:(ncol(tumorbaf)-3),FUN=mpi_ascat,tumorlogr=tumorlogr,tumorbaf=tumorbaf,normallogr=normallogr,normalbaf=normalbaf,job.num=njobs)
res=mpi.parSapply(X=1:nrow(allsamples),FUN=mpi_ascat1,job.num=njobs)


## quit program
mpi.close.Rslaves()
mpi.quit()
quit()

#Six failed samples:
failedsamples=c("TCGA-13-1481","TCGA-23-1122","TCGA-25-1323","TCGA-29-2414","TCGA-61-2111","TCGA-23-1119")
failedjobn=which(allsamples[,1] %in% failedsamples)
