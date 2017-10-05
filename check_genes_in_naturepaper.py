#!/usr/bin/env python
import re

with open("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/naturepaper.txt","r") as myfile:
	naturedata=myfile.read().replace("\n","")

genes=[]
with open("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/genes_smallqvalues","r") as myfile:
	for line in myfile:
		line=line.replace("\r","")
		line=line.replace("\n","")
		genes.append(line)

count=[]
for i in range(len(genes)):
	res=re.findall("\w*"+genes[i].lower()+"\w*",naturedata.lower())
	count.append(len(res))

for i in range(len(genes)):
	if count[i]>0:
		print i
		print genes[i]

