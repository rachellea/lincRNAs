#!//bin/bash

#This script is called with one argument: the name of the GWAS file to be
#cleaned, without the .txt extension (e.g. if your file is named
#GLGC_HDL_data.txt, then call this script using
#bash cleanSNPdataset.sh GLGC_HDL_data

#The GWAS file to be cleaned must be tab-delimited, with SNP rs number as the
#first column, and all remaining data related to that SNP contained in the
#second column (use a non-whitespace delimiter between the fields in the
#second column, if you want to include multiple fields in the second column;
#e.g. input file formatted as:
#rs1       pval1*beta1*allele1
#rs2       pval2*beta2*allele2

#This input file does not need to be sorted by anything. If you get a "file not
#in sorted order" error, this could have to do with incompatibilities in your
#versions of 'sort' and 'join'. If you see that error, add
#LANG=en_EN
#before every occurrance of 'sort' and 'join', and that should fix your problem.

#This script makes use of the files RsMergeArch_ready.txt and snp137_ready.txt
#Please put those in the same directory in which you are running this script

#Also, a warning: if you have any files beginning with "DEL" in the same
#directory as this script, they will be automatically deleted at the end
#of this script unless you move them somewhere else temporarily.

#The output will be a file called e.g. GLGC_HDL_data_clean.txt, and will contain
#seven columns:
#column 1: chromosome
#column 2: SNP coordinate in hg19
#column 3: SNP coordinate in hg19
#column 4: 0
#column 5: 0
#column 6: rs number without the rs prefix, in hg19
#column 7: all data corresponding to that SNP that was provided in column 2 of
#the input file (e.g. pval*beta*allele)
#This output file format can be fed directly into the program ANNOVAR.

GWASfile=${1}

touch ${GWASfile}_cleaning_logfile.txt

#FIRST: join GWAS file with RsMergeArch to update the rs numbers================
sort -k1,1 ${GWASfile}.txt > DEL_${GWASfile}_sort.txt
join -1 1 -2 1 -t "	" -o 2.2 1.2 DEL_${GWASfile}_sort.txt RsMergeArch_ready.bcp >DEL_${GWASfile}_1joined.txt
join -a1 -v1 -1 1 -2 1 -t "	" DEL_${GWASfile}_sort.txt RsMergeArch_ready.bcp > DEL_${GWASfile}_1nonjoined.txt
cat DEL_${GWASfile}_1joined.txt DEL_${GWASfile}_1nonjoined.txt > DEL_${GWASfile}_1.txt
echo "Number of lines in ${filename}.txt is " `cat ${GWASfile}.txt | wc -l` > ${GWASfile}_cleaning_logfile.txt
echo "Number of lines remaining after joining with RsMergeArch (to update the rs numbers) should be the same. It is " `cat DEL_${GWASfile}_1.txt | wc -l` >> ${GWASfile}_cleaning_logfile.txt
echo "There were " `cat DEL_${GWASfile}_1joined.txt | wc -l` " SNPs that had their rs numbers updated." >> ${GWASfile}_cleaning_logfile.txt
#note: use cat filename | wc -l instead of wc -l filename because the second way results in the filename being printed as well

#SECOND: give chromosome and position using the processed snp137 file===========
#first have to sort GWASfile again since it's no longer sorted due to cat
sort -k1,1 DEL_${GWASfile}_1.txt > DEL_${GWASfile}_1_sort.txt
join -1 1 -2 4 -t "	" -o 1.2 1.1 2.1 2.3 DEL_${GWASfile}_1_sort.txt snp137_ready.txt > DEL_${GWASfile}_2pre.txt
echo "After joining with snp137_ready to give chr and pos, the  number of lines left is " `cat DEL_${GWASfile}_2pre.txt | wc -l`>> ${GWASfile}_cleaning_logfile.txt
#Delete X, Y, and M chromosomes
grep -v -e X -e Y -e M DEL_${GWASfile}_2pre.txt > DEL_${GWASfile}_2.txt

#THIRD: filtering===============================================================
#first delete all duplicate positions (chr and coord)
sort -k3,4 DEL_${GWASfile}_2.txt | uniq -u -f 2 > DEL_${GWASfile}_3.txt
echo "After deleting all duplicate positions (chr and coord), the number of lines left is " `cat DEL_${GWASfile}_3.txt | wc -l` >> ${GWASfile}_cleaning_logfile.txt
#then rearrange the order so rs number is at the end, sort by rs, and delete all duplicate rs numbers
awk '{print $1, "\t", $3, "\t", $4, "\t", $2}' DEL_${GWASfile}_3.txt | sort -k4,4 | uniq -u -f 3 > DEL_${GWASfile}_4.txt
echo "After deleting all the duplicate rs numbers, the number of lines left is " `cat DEL_${GWASfile}_4.txt | wc -l` >> ${GWASfile}_cleaning_logfile.txt
#current order of columns is pval, chr, coordinate, rsnumber

#FOURTH: put into annovar format; delete files==================================
awk '{print $2,"\t",$3,"\t",$3,"\t",0,"\t",0,"\t",$4,"\t",$1}' DEL_${GWASfile}_4.txt | tr -d "chr" > DEL_${GWASfile}_clean.txt

#Renaming to protect final file from deletion:
mv DEL_${GWASfile}_clean.txt ${GWASfile}_clean.txt

rm DEL_*
