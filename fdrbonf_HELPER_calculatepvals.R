#! /usr/bin/Rscript
#This script is called by doing ./<scriptname>.R <argument1> <argument2>
#<argument1> is a path to a tab-delimited file with no header containing two columns:
#SNP rs number and SNP p-value.
#<argument2> is a path indicating the location and file name for this script's output file.

args<-commandArgs(trailingOnly=TRUE)
SNPpvaluniqfile<-args[1]
outputfilename<-args[2]

thedata<-read.table(SNPpvaluniqfile,header=FALSE) #SNP rs number col will be named V1; p-val col will be V2
pbonf<-p.adjust(thedata$V2,"bonferroni")
cat("The number of pbonfs<0.05 is",sum(pbonf<0.05),"\n")
pfdr<-p.adjust(thedata$V2,"fdr")
cat("The number of pfdrs<0.05 is",sum(pfdr<0.05),"\n")
thedata["p_bonf"]<-pbonf #add column to table
thedata["p_fdr"]<-pfdr
write.table(thedata,file=outputfilename,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
cat("Done with R script\n")
rm(thedata)
rm(pbonf)
rm(pfdr)
