#!/bin/bash
#Rachel Ballantyne
#April 2014

#For one "mixmap input file" showing GWAS SNPs annotated to genomic features,
#    1) Calculate the FDR and Bonferroni corrected p-vals for each SNP based on 
#       the total number of UNIQUE SNPs
#    2) Obtain the counts of the number of FDR-significant and Bonferroni-
#       significant SNPs within each feature
#INPUT: Script is called with two arguments. The first is a simple GWAS name, 
#e.g. GLGC_HDL. The second is a path to the mixmap input file for that GWAS, 
#e.g. /home/raba/stuff/GLGC_HDL_mixmapinput.txt
#The "mixmap input file" is tab-delimited with five columns:
#SNP rs number, gene name, chr, SNP coordinate, and SNP p-value.
#It describes SNPs annotated to genes.
#EXECUTION: You must be in a bsub session on PMACS in order for the cluster to 
#find the Rscript interpreter. You cannot be at consign.

set -e

gwas=$1
mixmapfile=$2

echo "Working on ${gwas}"

#Remove the header and select out rs number and p-val columns
sed '1d' ${mixmapfile} | awk '{print $1,"\t",$5}' | sort | uniq > ${gwas}_SNPpval-uniq.txt

#Create a file called ${gwas}_SNPpval-uniq_corrpvals.txt which contains the following columns:
#SNP rs number, raw p-val, p_bonf, p_fdr (this file does have a header: V1  V2  p_bonf  p_fdr)
./fdrbonf_HELPER_calculatepvals.R ${gwas}_SNPpval-uniq.txt ${gwas}_SNPpval-uniq_corrpvals.txt

#Add the p_bonf and p_fdr to the MixMAP file
python fdrbonf_HELPER_matchpvalswithgenes.py ${gwas}_SNPpval-uniq_corrpvals.txt ${mixmapfile} ${gwas}_7cols_allps.txt

#FDR------------------------
#Sort numerically according to p_fdr, and selected out all SNP-to-feature annotations with p_fdr < 0.05.
touch ${gwas}_7cols_fdrlt005.txt
sort -g -k7,7 ${gwas}_7cols_allps.txt | awk '$7 < 0.05' > ${gwas}_7cols_fdrlt005.txt
echo "${gwas} number of unique SNPs with p_fdr < 0.05 is"
awk '{print $1}' ${gwas}_7cols_fdrlt005.txt | sort | uniq | wc -l

#Take the SNP-to-feature annotations where the p_fdr < 0.05,
#count how many FDR-significant SNPs are in each feature,
#and rank the features according to this count.
awk '{print $2}' ${gwas}_7cols_fdrlt005.txt | sort | uniq -c | sort -k1,1 -g -r | cat -n > ${gwas}_3cols_fdrlt005_perfeat.txt

#BONFERRONI----------------
#Sort numerically according to p_bonf, and select out all SNP-to-feature annotations with p_bonf < 0.05
touch ${gwas}_7cols_bonflt005.txt
sort -g -k6,6 ${gwas}_7cols_allps.txt | awk '$6 < 0.05' > ${gwas}_7cols_bonflt005.txt
echo "${gwas} number of unique SNPs with p_bonf < 0.05 is"
awk '{print $1}' ${gwas}_7cols_bonflt005.txt | sort | uniq | wc -l

#Take the SNP-to-feature annotations where the p_fdr < 0.05,
#count how many Bonferroni-significant SNPs are in each feature,
#and rank the features according to this count.
awk '{print $2}' ${gwas}_7cols_bonflt005.txt | sort | uniq -c | sort -k1,1 -g -r | cat -n > ${gwas}_3cols_bonflt005_perfeat.txt

echo "All done with ${gwas}"
printf "\n"




