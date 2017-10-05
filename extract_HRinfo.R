

#load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/firehose.RData") #rawmutation 461 samples somatic

library(GenomicRanges)
hr30genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/HR30genes.txt",sep="\t",header=T,stringsAsFactors=F)
hr30genes[c(17,23),]
hr30genes$Gene=gsub(" (\\S+)","",hr30genes$Gene,perl=T)
hr30genes[c(17,23),]
hr30genes$chr=rep(NA,nrow(hr30genes))
hr30genes$start=rep(NA,nrow(hr30genes))
hr30genes$end=rep(NA,nrow(hr30genes))
for (i in 1:nrow(hr30genes))
{
  tmp=unlist(strsplit(hr30genes$Regions.targeted.for.capture...UCSC.hg19.[i],":"))
  hr30genes$chr[i]=tmp[1]
  tmp1=unlist(strsplit(tmp[2],"-"))
  hr30genes$start[i]=as.integer(tmp1[1])
  hr30genes$end[i]=as.integer(tmp1[2])
}
hr30genes=hr30genes[,-2]
hr30genes$chr=gsub("chr","",hr30genes$chr)
gr_hr30genes=GRanges(seqnames = hr30genes$chr,ranges=IRanges(start=hr30genes$start,end=hr30genes$end))
#TCGA------------
library(data.table)

filelist=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/GDCdata/copynumber/snp6/estimate_copynumber/filelist.txt",header=T,sep="\t",stringsAsFactors=F)
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/platinum_classificationdata_stringent_filtered.RData")
allsamples=rownames(data_copynumber_filtered)
mutationfiles=list.files(path="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/wustl_protectedmutation/genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.1.0/",patter="genome.*.snv.*.vcf")

getmutations=function(jobn)
{
  chrs=c(1:22,"X")
  samplename=allsamples[jobn]
  res_somatic=res_germline=data.frame(matrix(NA,nrow=1,ncol=nrow(hr30genes)))
  colnames(res_germline)=colnames(res_somatic)=hr30genes$Gene
  errmsg=""
  filenames=mutationfiles[which(grepl(samplename,mutationfiles))]
  if (length(filenames)==0) errmsg="no mutation files found"
  if (errmsg=="")
  {
    filesizes=file.size(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/wustl_protectedmutation/genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.1.0/",filenames))
    filename=filenames[which.max(filesizes)]
  }
  if (errmsg=="")
  {
    if (length(filename)>0)
    {
      res_somatic[1,]=0
      res_germline[1,]=0
      #the mutation calls
      mutations=fread(paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/controlleddata/wustl_protectedmutation/genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_Cont_automated.Level_2.1.1.0/",filename),
                      skip=57,sep="\t",stringsAsFactors = F)
      mutations=as.data.frame(mutations)
      mutations=mutations[mutations$`#CHROM` %in% chrs,]
      # #pick mutation types
      SS=sapply(1:nrow(mutations),function(i){
        res1=NA
        tmp=unlist(strsplit(mutations$FORMAT[i],":",fixed=T))
        idxss=which(tmp=="SS")
        if (length(idxss)>0)
        {
          tmp=unlist(strsplit(mutations[i,10],":",fixed=T))
          res1=tmp[idxss]
        }
        res1
      })
      

      
      # #check status
      FT1=sapply(1:nrow(mutations),function(i){
        res1=FALSE
        tmp=unlist(strsplit(mutations$FORMAT[i],":",fixed=T))
        idxft=which(tmp=="FT")
        if (length(idxft)>0)
        {
          tmp=unlist(strsplit(mutations[i,10],":",fixed=T))
          tmp1=tmp[idxft]
          if (!grepl("fa",tolower(tmp1))) res1=TRUE #no failure or false
        }
        res1
      })
      
      FT=sapply(1:nrow(mutations),function(i){
        res1=FALSE
        tmp=unlist(strsplit(mutations$FORMAT[i],":",fixed=T))
        idxft=which(tmp=="FT")
        if (length(idxft)>0)
        {
          tmp=unlist(strsplit(mutations[i,10],":",fixed=T))
          tmp1=tmp[idxft]
          if (grepl("pass",tolower(tmp1))) res1=TRUE #no failure or false
        }
        res1
      })
      sum(FT)/length(FT)
      
      idxgermline=which(SS %in% c(1,3) & FT)
      germline=data.frame(chr=mutations$`#CHROM`[idxgermline],pos=mutations$POS[idxgermline],ref=mutations$REF[idxgermline],alt=mutations$ALT[idxgermline])
      germline$ref=gsub(",\\w+$","",germline$ref,perl=T)
      germline$alt=gsub(",\\w+$","",germline$alt,perl=T)
      mutationfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tmp_",samplename,".txt")
      write.table(germline,mutationfile,col.names = T,row.names = F,sep="\t",quote=F)
      system(paste0("/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator.sh ",mutationfile))
      germline_anno=read.table(paste0(mutationfile,".oncotator"),skip=2,header=T,sep="\t",fill = T,quote="")
      colnames(germline_anno)[5:6]=c("chr","pos")
      germline_anno=merge(germline_anno,germline,by=c("chr","pos"))
      idx_function=which(!germline_anno$Variant_Classification %in% c("Intron","5'UTR","3'UTR","IGR","5'Flank","Silent"))
      germline_anno=germline_anno[idx_function,]
      gr_germline=GRanges(seqnames =germline_anno$chr, ranges = IRanges(start=germline_anno$pos,width = 1) )
      #gr_germline=GRanges(seqnames =germline$chr, ranges = IRanges(start=germline$pos,width = 1) )
      tmp=distanceToNearest(gr_hr30genes,gr_germline)
      idx=which(mcols(tmp)$distance==0)
      if (length(idx)>0) res_germline[,idx]=1
      idxsomatic=which(SS=="2" & FT)
      somatic=data.frame(chr=mutations$`#CHROM`[idxsomatic],pos=mutations$POS[idxsomatic],ref=mutations$REF[idxsomatic],alt=mutations$ALT[idxsomatic])
      somatic$ref=gsub(",\\w+$","",somatic$ref,perl=T)
      somatic$alt=gsub(",\\w+$","",somatic$alt,perl=T)
      write.table(somatic,mutationfile,col.names = T,row.names = F,sep="\t",quote=F)
      system(paste0("/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator.sh ",mutationfile))
      somatic_anno=read.table(paste0(mutationfile,".oncotator"),skip=2,header=T,sep="\t",fill = T,quote="")
      colnames(somatic_anno)[5:6]=c("chr","pos")
      somatic_anno=merge(somatic_anno,somatic,by=c("chr","pos"))
      if (file.exists(mutationfile)) system(paste0("rm ",mutationfile))
      #if (file.exists(paste0(mutationfile,".oncotator"))) system(paste0("rm ",mutationfile,".oncotator"))
      idx_function=which(!somatic_anno$Variant_Classification %in% c("Intron","5'UTR","3'UTR","IGR","5'Flank","Silent"))
      somatic_anno=somatic_anno[idx_function,]
      gr_somatic=GRanges(seqnames =somatic_anno$chr, ranges = IRanges(start=somatic_anno$pos,width = 1) )
      #gr_somatic=GRanges(seqnames =somatic$chr, ranges = IRanges(start=somatic$pos,width = 1) )
      tmp=distanceToNearest(gr_hr30genes,gr_somatic)
      idx=which(mcols(tmp)$distance==0)
      if (length(idx)>0) res_somatic[,idx]=1
      
    }
  }else
  {
    print(paste0(samplename,": ",errmsg))
  }
  return(list(germline=res_germline,somatic=res_somatic))
}

missingsamples=function()
{
  missingsamples=NULL
  for (jobn in 1:length(allsamples))
  {
    samplename=allsamples[jobn]
    filenames=mutationfiles[which(grepl(samplename,mutationfiles))]
    if (length(filenames)==0) missingsamples=c(missingsamples,samplename)
  }
  
}
#salloc -t 0-5 -n 90 mpirun -n 1 R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(hr30genes)
mpi.bcast.Robj2slave(gr_hr30genes)
mpi.bcast.Robj2slave(mutationfiles)
mpi.bcast.Robj2slave(allsamples)
mpi.bcast.cmd(library(data.table))
mpi.bcast.cmd(library(GenomicRanges))

res=mpi.parSapply(X=1:length(allsamples),FUN=getmutations,job.num=njobs)
save(res,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_HRmutation2.Rdata")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_HRmutation2.Rdata")
res1=matrix(unlist(res),byrow = F,nrow=nrow(hr30genes))
germline1=res1[,seq(1,ncol(res1),2)]
somatic1=res1[,seq(2,ncol(res1),2)]
rownames(germline1)=rownames(somatic1)=hr30genes$Gene
colnames(germline1)=colnames(somatic1)=allsamples
germline1=as.data.frame(germline1)
somatic1=as.data.frame(somatic1)


somatic=germline=data.frame(matrix(NA,nrow=nrow(hr30genes),ncol=length(allsamples)))
rownames(somatic)=rownames(germline)=hr30genes$Gene
colnames(somatic)=colnames(germline)=allsamples
for (i in 1:length(allsamples))
{
  if (i %% 10==0) cat(i,"..")
  mysample=allsamples[i]
  tmp=getmutations(i)
  somatic[,i]=as.numeric(tmp$somatic)
  germline[,i]=as.numeric(tmp$germline)
}

save(germline,somatic,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_HRmutation.Rdata")
#AUS-------------

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/aocs_data.RData")
#find germline vcf files from sequencing data
allsamples=rownames(data_aocs_copynumber)
donortable=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/clinical/donor.tsv",header=T,sep="\t",stringsAsFactors = F)
manifesttable=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/manifest.collaboratory.1480527008951.tsv",header=T,sep="\t",stringsAsFactors = F)
manifestable_london=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/manifest.pcawg-london.1480448530358.xml",header=F,sep="|",stringsAsFactors = F)
idx_analysis=which(grepl("<analysis_id>",manifestable_london[,1])==T)

germlinevcffiles=rep(NA,length(allsamples))
for (i in 1:length(allsamples))
{
  cat(i,"..")
  mysample=allsamples[i]
  mysample1=unlist(strsplit(mysample,"-",fixed=T))
  mysample1=paste0(mysample1[1:2],collapse = "-")
  idx=which(donortable$submitted_donor_id==mysample1)
  mydonor=donortable$icgc_donor_id[idx]
  idx1=which(manifesttable$donor_id.donor_count==mydonor)
  for (j in idx1)
  {
    md5=manifesttable$md5_sum[j]
    idx2=which(grepl(md5,manifestable_london[,1])==T)
    if (length(idx2)>0)
    {
      myfile=manifestable_london[idx2-2,1]
      myfile=gsub(" ","",myfile)
      myfile=gsub("<filename>","",myfile)
      myfile=gsub("</filename>","",myfile)
      if (grepl("germline.snv_mnv.vcf",myfile)==T)
      {
        idx3=which(idx_analysis<idx2)
        idx3=idx_analysis[length(idx3)]
        myfolder=manifestable_london[idx3,1]
        myfolder=gsub(" ","",myfolder)
        myfolder=gsub("<analysis_id>","",myfolder)
        myfolder=gsub("</analysis_id>","",myfolder)
        cnvfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/",myfolder,"/",myfile)
        cnvfile=gsub(".gz$","",cnvfile)
        if (!file.exists(cnvfile))
          system(paste0("gunzip ",paste0(cnvfile),".gz"))
        
        germlinevcffiles[i]=cnvfile
        break
      }
    }
  }
}
names(germlinevcffiles)=allsamples

somaticvcffiles=rep(NA,length(allsamples))
for (i in 1:length(allsamples))
{
  cat(i,"..")
  mysample=allsamples[i]
  mysample1=unlist(strsplit(mysample,"-",fixed=T))
  mysample1=paste0(mysample1[1:2],collapse = "-")
  idx=which(donortable$submitted_donor_id==mysample1)
  mydonor=donortable$icgc_donor_id[idx]
  idx1=which(manifesttable$donor_id.donor_count==mydonor)
  for (j in idx1)
  {
    md5=manifesttable$md5_sum[j]
    idx2=which(grepl(md5,manifestable_london[,1])==T)
    if (length(idx2)>0)
    {
      myfile=manifestable_london[idx2-2,1]
      myfile=gsub(" ","",myfile)
      myfile=gsub("<filename>","",myfile)
      myfile=gsub("</filename>","",myfile)
      if (grepl("germline.snv_mnv.vcf",myfile)==T)
      {
        idx3=which(idx_analysis<idx2)
        idx3=idx_analysis[length(idx3)]
        myfolder=manifestable_london[idx3,1]
        myfolder=gsub(" ","",myfolder)
        myfolder=gsub("<analysis_id>","",myfolder)
        myfolder=gsub("</analysis_id>","",myfolder)
        myfile=gsub("germline","somatic",myfile)
        cnvfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/Australia/ICGC/",myfolder,"/",myfile)
        cnvfile=gsub(".gz$","",cnvfile)
        if (!file.exists(cnvfile) & file.exists(paste0(cnvfile,".gz")))
          system(paste0("gunzip ",paste0(cnvfile),".gz"))
        
        somaticvcffiles[i]=cnvfile
        break
      }
    }
  }
}
names(somaticvcffiles)=allsamples

getausmutations=function(samplename)
{
  res_somatic=res_germline=data.frame(matrix(NA,nrow=1,ncol=nrow(hr30genes)))
  colnames(res_germline)=colnames(res_somatic)=hr30genes$Gene
  idx=which(allsamples==samplename)
  if (!is.na(germlinevcffiles[idx]))
  {
    res_somatic[1,]=0
    res_germline[1,]=0
    germline=read.table(file=germlinevcffiles[idx],skip=42,header=T,sep="\t",comment.char = "",stringsAsFactors = F)
    germline=germline[germline$FILTER=="PASS",]
    
    gr_germline=GRanges(seqnames =germline$X.CHROM, ranges = IRanges(start=germline$POS,width = 1) )
    tmp=distanceToNearest(gr_germline,gr_hr30genes)
    idx1=which(mcols(tmp)$distance==0)
    germline=germline[idx1,]
    germline=germline[,c(1,2,4,5)]
    colnames(germline)=c("chr","pos","ref","alt")
    mutationfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tmp_",samplename,".txt")
    write.table(germline,mutationfile,col.names = T,row.names = F,sep="\t",quote=F)
    system(paste0("/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator.sh ",mutationfile))
    germline_anno=read.table(paste0(mutationfile,".oncotator"),skip=2,header=T,sep="\t",fill = T,quote="")
    colnames(germline_anno)[5:6]=c("chr","pos")
    germline_anno=merge(germline_anno,germline,by=c("chr","pos"))
    idx_function=which(!germline_anno$Variant_Classification %in% c("Intron","5'UTR","3'UTR","IGR","5'Flank","Silent"))
    germline_anno=germline_anno[idx_function,]
    gr_germline=GRanges(seqnames =germline_anno$chr, ranges = IRanges(start=germline_anno$pos,width = 1) )
  
    tmp=distanceToNearest(gr_hr30genes,gr_germline)
    idx1=which(mcols(tmp)$distance==0)
    if (length(idx1)>0) res_germline[,idx1]=1
    
    somatic=read.table(file=somaticvcffiles[idx],skip=42,header=T,sep="\t",comment.char = "",stringsAsFactors = F)
    somatic=somatic[somatic$FILTER=="PASS",]
    gr_somatic=GRanges(seqnames =somatic$X.CHROM, ranges = IRanges(start=somatic$POS,width = 1) )
    tmp=distanceToNearest(gr_somatic,gr_hr30genes)
    idx1=which(mcols(tmp)$distance==0)
    somatic=somatic[idx1,]
    somatic=somatic[,c(1,2,4,5)]
    colnames(somatic)=c("chr","pos","ref","alt")
    mutationfile=paste0("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tmp_",samplename,".txt")
    write.table(somatic,mutationfile,col.names = T,row.names = F,sep="\t",quote=F)
    system(paste0("/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator.sh ",mutationfile))
    somatic_anno=read.table(paste0(mutationfile,".oncotator"),skip=2,header=T,sep="\t",fill = T,quote="")
    colnames(somatic_anno)[5:6]=c("chr","pos")
    somatic_anno=merge(somatic_anno,somatic,by=c("chr","pos"))
    idx_function=which(!somatic_anno$Variant_Classification %in% c("Intron","5'UTR","3'UTR","IGR","5'Flank","Silent"))
    somatic_anno=somatic_anno[idx_function,]
    gr_somatic=GRanges(seqnames =somatic_anno$chr, ranges = IRanges(start=somatic_anno$pos,width = 1) )
    if (file.exists(mutationfile)) system(paste0("rm ",mutationfile))
    if (file.exists(paste0(mutationfile,".oncotator"))) system(paste0("rm ",mutationfile,".oncotator"))
    tmp=distanceToNearest(gr_hr30genes,gr_somatic)
    idx1=which(mcols(tmp)$distance==0)
    if (length(idx1)>0) res_somatic[,idx1]=1
  
  }
  return(list(germline=res_germline,somatic=res_somatic))
}

somatic_aus=germline_aus=data.frame(matrix(NA,nrow=nrow(hr30genes),ncol=length(allsamples)))
rownames(somatic_aus)=rownames(germline_aus)=hr30genes$Gene
colnames(somatic_aus)=colnames(germline_aus)=allsamples
for (i in 1:length(allsamples))
{
  if (i %% 10==0) cat(i,"..")
  mysample=allsamples[i]
  tmp=getausmutations(mysample)
  somatic_aus[,i]=as.numeric(tmp$somatic)
  germline_aus[,i]=as.numeric(tmp$germline)
}

save(germline_aus,somatic_aus,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_HRmutation2.Rdata")

load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_HRmutation.Rdata")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/aus_HRmutation.Rdata")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp/tcga_ascat_loh_cytoband_0.25.RData")
brca1_s=data.frame(sampleID=colnames(somatic),brca1=as.numeric(somatic[6,]))
resist=merge(resist,brca1_s,by="sampleID")
summary(glm(I(resistance=="Sensitive")~brca1+loh_length,family = binomial,data=resist))
brca1_g=data.frame(sampleID=colnames(somatic),brca1_g=as.numeric(somatic[6,]))
resist=merge(resist,brca1_g,by="sampleID")
summary(glm(I(resistance=="Sensitive")~brca1_g+loh_length,family = binomial,data=resist))

#use the suplementary tables
somatictable=read.csv("/fh/fast/dai_j/CancerGenomics/Ovarian/TCGAdata/Mutation/TCGA-OV-mutations/TCGA-OV-Final-Supp-Table-S2.1-13jan2011g.csv",header=T,stringsAsFactors = F)
samples=sapply(1:nrow(somatictable),function(i){
  tmp=unlist(strsplit(somatictable$Tumor_Sample_Barcode[i],"-"))
  res1=paste0(tmp[1:3],collapse = "-")
})
allsamples1=unique(samples)
length(intersect(allsamples,allsamples1))
#[1] 217 overlap

tmp=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_BRCA1_germline.txt",fill=T,stringsAsFactors = F)
brca1_g=NULL
for (i in 1:nrow(tmp))
{
  if (grepl("TCGA",tmp$V1[i])) brca1_g=c(brca1_g,tmp$V1[i])
}
tmp=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_BRCA1_somatic.txt",fill=T,stringsAsFactors = F)
brca1_s=NULL
for (i in 1:nrow(tmp))
{
  if (grepl("TCGA",tmp$V1[i])) brca1_s=c(brca1_s,tmp$V1[i])
}
tmp=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_BRCA2_germline.txt",fill=T,stringsAsFactors = F)
brca2_g=NULL
for (i in 1:nrow(tmp))
{
  if (grepl("TCGA",tmp$V1[i])) brca2_g=c(brca2_g,tmp$V1[i])
}
tmp=read.table("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_BRCA2_somatic.txt",fill=T,stringsAsFactors = F)
brca2_s=NULL
for (i in 1:nrow(tmp))
{
  if (grepl("TCGA",tmp$V1[i])) brca2_s=c(brca2_s,tmp$V1[i])
}
formmutationtable=function(samples=allsamples1,mutationsamples=brca1_g,name="brca1_germline")
{
  res1=data.frame(matrix(0,nrow=length(samples),ncol=1))
  colnames(res1)=name
  rownames(res1)=samples
  res1[samples %in% mutationsamples,1]=1
  return(res1)
}
brcamutation=formmutationtable()
brcamutation=cbind(brcamutation,formmutationtable(mutationsamples = brca1_s,name="brca1_somatic"))
brcamutation=cbind(brcamutation,formmutationtable(mutationsamples = brca2_g,name="brca2_germline"))
brcamutation=cbind(brcamutation,formmutationtable(mutationsamples = brca2_s,name="brca2_somatic"))
brcamutation=cbind(sampleID=allsamples1,brcamutation)
save(brcamutation,file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_BRCA_mutations.RData")
load("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/clinical_371samples.RData")
dat=merge(clinicaltable,brcamutation,by="sampleID")
sum(dat$brca1_germline==1 |dat$brca1_somatic==1|dat$brca2_germline==1|dat$brca2_somatic==1)
#[1] 46
46/nrow(dat)
#[1] 0.2119816
fit=glm(I(resistance=="Sensitive")~I(brca1_germline==1 | brca1_somatic==1 |brca2_germline==1 |brca2_somatic==1),data=dat,family = binomial)
summary(fit)
#                                                                                           Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                                                                  0.5896     0.1596   3.693 0.000221 ***
#  I(brca1_germline == 1 | brca1_somatic == 1 | brca2_germline == 1 | brca2_somatic == 1)TRUE   1.3075     0.4660   2.806 0.005018 **
summary(glm(I(resistance=="Sensitive")~I(brca1_germline==1),data=dat,family = binomial))
summary(glm(I(resistance=="Sensitive")~I(brca1_somatic==1),data=dat,family = binomial))
summary(glm(I(resistance=="Sensitive")~I(brca2_germline==1),data=dat,family = binomial))
#(Intercept)                   0.723      0.150   4.820 1.44e-06 ***
#I(brca2_germline == 1)TRUE    1.916      1.046   1.832    0.067 .  
summary(glm(I(resistance=="Sensitive")~I(brca2_somatic==1),data=dat,family = binomial))
summary(glm(I(resistance=="Sensitive")~I(brca1_germline==1 | brca1_somatic==1),data=dat,family = binomial))
summary(glm(I(resistance=="Sensitive")~I(brca2_germline==1 |brca2_somatic==1),data=dat,family = binomial))
#(Intercept)                                       0.6702     0.1513   4.428  9.5e-06 ***
#  I(brca2_germline == 1 | brca2_somatic == 1)TRUE   2.3744     1.0346   2.295   0.0217 * 
summary(glm(I(resistance=="Sensitive")~I(brca1_germline==1 |brca2_germline==1),data=dat,family = binomial))
#(Intercept)                                        0.6604     0.1564   4.222 2.42e-05 ***
#I(brca1_germline == 1 | brca2_germline == 1)TRUE   1.1314     0.5077   2.228   0.0259 *  
summary(glm(I(resistance=="Sensitive")~I(brca1_somatic==1 |brca2_somatic==1),data=dat,family = binomial))
summary(glm(I(resistance=="Sensitive")~brca1_germline+brca1_somatic+brca2_germline+brca2_somatic,data=dat,family = binomial))
#(Intercept)      0.5877     0.1595   3.685 0.000228 ***
#brca1_germline   0.8089     0.5800   1.395 0.163100    
#brca1_somatic    0.5672     1.1519   0.492 0.622417    
#brca2_germline   1.9810     1.0480   1.890 0.058722 .  
#brca2_somatic   15.9784   906.9427   0.018 0.985944  
fit=glm(I(resistance=="Sensitive")~I(brca1_germline==1 | brca1_somatic==1 |brca2_germline==1 |brca2_somatic==1),data=dat[dat$before2009==1,],family = binomial)
summary(fit)
#(Intercept)                                                                                  0.4678     0.1689   2.769  0.00561 **
#I(brca1_germline == 1 | brca1_somatic == 1 | brca2_germline == 1 | brca2_somatic == 1)TRUE   1.2062     0.4759   2.535  0.01126 *
