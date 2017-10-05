#!/usr/bin/env python

import re
from os import listdir
from os.path import isfile, join
#class of vcf
#class vcf:
#	def __init__(self,filename):
#		self.filename=filename
	
def readvcf(filename):
	outputfile=filename+".frq.txt"
	f=open(outputfile,'w')
	s= "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "normal_ref" + "\t" + "tumor_ref" +"\t" + "normal_alt" + "\t" + "tumor_alt" + "\n"
	f.write(s)
	vcffile = open(filename, 'r')
	
	for line in vcffile:
		if not re.search("^#",line):
			lines=line.split("\t")
			chr=lines[0]
			pos=lines[1]
			ref=lines[3]
			alt=lines[4]
			filter=lines[6]
			if filter=="PASS":
				normals=lines[9].split(":")
				tumors=lines[10].split(":")			
				if ref=="A":
					normal_ref=normals[8].split(",")
					normal_ref=normal_ref[0]				
					tumor_ref=tumors[8].split(",")
					tumor_ref=tumor_ref[0]
				elif ref=="C":
					normal_ref=normals[9].split(",")
					normal_ref=normal_ref[0]				
					tumor_ref=tumors[9].split(",")
					tumor_ref=tumor_ref[0]
				elif ref=="G":
					normal_ref=normals[10].split(",")
					normal_ref=normal_ref[0]				
					tumor_ref=tumors[10].split(",")
					tumor_ref=tumor_ref[0]
				else:
					normal_ref=normals[11].split(",")
					normal_ref=normal_ref[0]				
					tumor_ref=tumors[11].split(",")
					tumor_ref=tumor_ref[0]
				if alt=="A":
					normal_alt=normals[8].split(",")
					normal_alt=normal_alt[0]				
					tumor_alt=tumors[8].split(",")
					tumor_alt=tumor_alt[0]
				elif alt=="C":
					normal_alt=normals[9].split(",")
					normal_alt=normal_alt[0]				
					tumor_alt=tumors[9].split(",")
					tumor_alt=tumor_alt[0]
				elif alt=="G":
					normal_alt=normals[10].split(",")
					normal_alt=normal_alt[0]				
					tumor_alt=tumors[10].split(",")
					tumor_alt=tumor_alt[0]
				else:
					normal_alt=normals[11].split(",")
					normal_alt=normal_alt[0]				
					tumor_alt=tumors[11].split(",")
					tumor_alt=tumor_alt[0]
				s= chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + normal_ref + "\t" + tumor_ref +"\t" + normal_alt + "\t" + tumor_alt + "\n"
				f.write(s)
	vcffile.close()
	f.close()

files = [f for f in listdir('../data/controlleddata/bcgsc_protectedmutation/bcgsc.ca_OV.Multicenter_mutation_calling_MC3_Cont.Level_2.1.0.0') if isfile(join('../data/controlleddata/bcgsc_protectedmutation/bcgsc.ca_OV.Multicenter_mutation_calling_MC3_Cont.Level_2.1.0.0', f))]
for myfile in files:
	if re.search("snv.vcf$",myfile):
		print myfile		
		myfile='../data/controlleddata/bcgsc_protectedmutation/bcgsc.ca_OV.Multicenter_mutation_calling_MC3_Cont.Level_2.1.0.0/'+myfile
		readvcf(myfile)


