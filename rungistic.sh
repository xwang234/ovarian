#!/usr/bin/env bash
# run in /fh/fast/dai_j/CancerGenomics/Tools/GISTIC

flag=1
rx=0
#armpeel=0
#firehose use 1
armpeel=1
#brlen=0.98
#firehose use 0.7
brlen=0.7
broad=1
maxseg=3000
#conf=0.95
#firehose use 0.99
conf=0.99
#cnvfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/CNV.hg19.bypos.111213.txt
#firehose use:
cnvfile=/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/examplefiles/SNP6.merged.151117.hg19.CNV.txt

#for TCGA study:
task=2
if [[ $task -eq 1 ]]
then
echo "TCGA"
markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt
if [[ $flag -eq 1 ]]
then
	
	prefix="sensitive"
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_sensitiveseg.txt
	#markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_sensitiveseg_marker.txt

	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi
if [[ $flag -eq 2 ]]
then
	prefix="resistant"
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_resistantseg.txt
	#markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_resistantseg_marker.txt
	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi
if [[ $flag -eq 3 ]]
then
	
	prefix="supersensitive"
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_supersensitiveseg.txt
	#markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_sensitiveseg_marker.txt

	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi
if [[ $flag -eq 4 ]]
then
	prefix="refractory"
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_refractoryseg.txt
	#markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/snp6copynumber_resistantseg_marker.txt
	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi

fi
if [[ $task -eq 2 ]]
then
	echo "AOCS"
	prefix="AOCS"
	markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.markers.txt
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/aocs.cbs.txt
	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi

if [[ $task -eq 3 ]]
then
	echo "TCGA"
	prefix="TCGA_new" #segment generated by cbs
	markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA.markers.txt
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/TCGA.tangent.cbs.txt
	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi

if [[ $task -eq 4 ]]
then
	echo "TCGA"
	prefix="TCGA" #segment generated by TCGA
	markersfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt
	segfile=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/firehose_segment.txt
	./gistic_firehose.sh $rx $maxseg $conf $armpeel $brlen $broad $prefix $segfile $markersfile $cnvfile
fi


