#!/usr/bin/env Rscript
library(gage)
library(pathview)
library(gageData)

data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
dcis=dcis[2:length(dcis)]
data(kegg.gs)
data(go.gs)

#kegg test for 1-directional changes
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, 
                        ref = hn, samp = dcis,compare='unpaired')
str(gse16873.kegg.p)
head(gse16873.kegg.p$greater)
head(gse16873.kegg.p$less)
head(gse16873.kegg.p$stats)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

data(egSymb)
sym2eg("TP53")
data=data_mrna
alldata=formwholedataframe(data)

#add indicator of train/test as the last column
#train=rep(FALSE,nrow(alldata))
#train[1:nrow(data$data1)]=TRUE
#alldata=cbind(alldata,train=train)

#remove grade 1(G1) in tumor_grade and Stage1 in clinical_stage
alldata=alldata[is.na(alldata[,"tumor_grade"]) | (! is.na(alldata[,"tumor_grade"]) & alldata[,"tumor_grade"]!="G1"),]
alldata=alldata[is.na(alldata[,"clinical_stage"]) | (! is.na(alldata[,"clinical_stage"]) & alldata[,"clinical_stage"]!="Stage1"),]

#clinicals to remove before analysis, they are not useful for the model. only age and residual_disease are kept
remove_clinicals=c("tumor_grade","clinical_stage","drug_interval_computed","race","initial_pathologic_dx_year","vital_status","death_months_to","treatment_outcome_first_course","progression_free_survival","uselastcontact")
alldata=alldata[,!colnames(alldata) %in% remove_clinicals]
residata=alldata[which(alldata[,"platinumclass"]=="Resistant"),4:ncol(alldata)]
sensdata=alldata[which(alldata[,"platinumclass"]=="Sensitive"),4:ncol(alldata)]
numsens=nrow(sensdata)
alldata=cbind(t(sensdata),t(sensdata))
enid=sapply(1:ncol(residata),function(i){
  res=unname(sym2eg(colnames(residata)[i]))
})
idx=which(!is.na(enid))
alldata=alldata[idx,]
rownames(alldata)=enid[idx]

data.kegg.p <- gage(alldata, gsets = kegg.sets.hs, use.fold=F,
                        ref = 1:numsens, samp = (numsens+1):ncol(alldata),compare='unpaired')
