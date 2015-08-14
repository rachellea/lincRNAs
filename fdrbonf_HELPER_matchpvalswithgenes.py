#Rachel Ballantyne
#April 2014

"""Match corrected p-values with the appropriate SNPs annotated to genomic
features.
INPUT: The script is called with three arguments.
-->>The first argument is a path to a tab-delimted file formatted with four
    columns and a header line:
    <SNP rs number>    <raw p value>    <p_bonf>    <p_fdr> 
-->>The second argument is a path to a tab-delimited file formatted with five
    columns and a header line:
    SNP rs number, gene name, chromosome, SNP coordinate, SNp p-value.
    E.g. rs10747505      AGL*NM_000642   1       100355560       0.581799
-->>The third argument is a path and filename where this script's output
    will be saved.

OUTPUT:
The output of this module is the merging of the contents of the files specified
by the first two arguments. There will be seven columns:
SNP rs number, gene name, chromosome, coordinate, raw p-val, p_bonf, p_fdr"""

import sys

def matchpvals(corrfilename,MMfilename,outputfilename):
    pvaldict={}
    corrpvals=open(corrfilename,'r')
    for line in corrpvals:
        lineaslist=line.rsplit()
        if lineaslist[0]=="V1":
            pass #ignore first line
        else:
            rsnumber=lineaslist[0]
            p_raw=lineaslist[1]
            p_bonf=lineaslist[2]
            p_fdr=lineaslist[3]
            pvaldict[rsnumber]=[p_raw,p_bonf,p_fdr]
    outputfile=open(outputfilename,'w')
    mixmapfile=open(MMfilename,'r')
    for line in mixmapfile:
        lineaslist=line.rsplit()
        if lineaslist[0]=="snp": #if reading the header line
            outputfile.write(line.replace("\n","")+"\tp_bonf\tp_fdr\n")
            #write a new header line in outputfile
        else:
            pvallist=pvaldict[lineaslist[0]] #look up rsnumber stored in
            #lineaslist[0] in the pvaldict created from corrpvals file
            assert float(pvallist[0])==float(lineaslist[4]), ("Error:"
                    +" pval for "+lineaslist[0]+" from mixmap file was "
                    +lineaslist[4]+" but pval from rightcorrpvals file was "
                    +pvallist[0])
            #assert that the raw p-val stored with the rs number in the dict
            #matches the raw p-val stored with the same rs number in the file
            xp_bonf=pvallist[1]
            xp_fdr=pvallist[2]
            outputfile.write(line.replace("\n","")+"\t"+xp_bonf+"\t"
                             +xp_fdr+"\n")

if __name__ == '__main__':
    corrfyle=sys.argv[1]
    MMfyle=sys.argv[2]
    outputfyle=sys.argv[3]
    matchpvals(corrfyle,MMfyle,outputfyle)
    print "Done with script 'fdrbonf_HELPER_matchpvalswithgenes.py'"
