#Module categorizeSNPsII.py
#Rachel Ballantyne

#Finish prepping this for publication - check all docstrings/comments especially

"""Create files of GWAS p-values based on the genomic features containing SNPs.

----------Module Input----------------------------------------------------------
This is an executable python module, called as follows:
    python categorizeSNPs.py <path to configuration file>
    e.g. python categorizeSNPs.py /home/raba/GWASready/categorizeSNPs_config.txt

Example config file:
#======================================#
# Specify the paths to the GWA Studies #
#======================================#
#list the names and paths to the GWAS's that will be analyzed
#Start each path on a new line with the format: GWAS: <GWASpath>

GWAS: /home/raba/GWAS_cleaned/Cardiogram_All/Cardiogram_All_clean.txt
GWAS: /home/raba/GWAS_cleaned/Cardiogram_MI/Cardiogram_MI_clean.txt

#=======================================================#
# Specify paths to files needed for making chrom arrays #
#=======================================================#
#the path to the file containing a pickled Python list of SmallFeature
#objects (e.g. representing lincRNAs):
lincRNApath: /home/raba/definedlincs.bin

#the path to the RefGene.txt file:
mRNApath: /home/raba/refGene.txt

#the path to the Gencode GTF file:
gencodepath: /home/raba/gencode.v18.annotation.gtf

#the path to the Hangauer S1 dataset file:
hang1path: /home/raba/Hangauer_S1_ExtendedProteinCodingGeneBoundaryFilter_HG19.txt

#the path to the Yale pseudogenes txt file:
pseudpath: /home/raba/Human_Pseudogene.txt

#the path to the nonintergenic regions pickled Python list:
nonintpath: /home/raba/categorizeSNPs_nonintergenic_regions_11-26_ext0.bin

#==============================================================#
# Specify additional parameters needed for making chrom arrays #
#==============================================================#
#an integer 0 <= bases <= 20000
bases: 0

#True or False, to determine whether actual chromosome sizes are used
userealsizes: True

#=================================================#
# Specify where the module's output will be saved #
#=================================================#
#the path to the dir in which the results dir will be saved (include final /)
outputlocation: /home/raba/Cardiogram_Redo/enrichment/
#END CONFIG FILE

Notes about config file format: the config file must include
-->>One or more lines beginning with "GWAS: " followed by the path to a
    tab-delimited GWAS data file with the following seven columns:
    chromosome as an int, coordinate as an int, coordinate as an int, 0, 0,
    rs number without the rs prefix, and GWAS p-value as a float.
    Each GWAS provided will be analyzed separately.
-->>Paths to various datasets.
    The module will initially assume that the relevant chromosome arrays have
    been created and that they are located in the same dir as the config file.
    If there are no chromosome arrays in the same dir as the config file, they
    will be created, making use of the following files which should be specified
    in the config file:
    -->>one line beginning with "lincRNApath: " followed by the path to a
        pickled Python list of SmallFeature objects that represent lincRNAs.
    -->>one line beginning with "mRNApath: " followed by the path to the file
        refGene.txt. 
    -->>one line beginning with "gencodepath: " followed by the path to
        Gencode GTF file
    -->>one line beginning with "hang1path: " followed by the path to the
        Hangauer S1 dataset
    -->>one line beginning with "pseudpath: " followed by the path to the file
        Human_Pseudogene.txt
    -->>one line beginning with "nonintpath: " followed by the path to the
        pickled Python list of "nonintergenic" regions
        (categorizeSNPs_NOTintergenic_locslist.bin) which was created using
        the lincRNA, mRNA, gencode, hang1file, and pseudogenes you specified.
    -->>one line beginning with "bases: " followed by the number of bases which
        mRNA and lincRNA locations will be expanded by.
    -->>one line beginning with "userealsizes: " followed by True or False.
        This line will be needed if the chromosome arrays have not been created
        before calling outputpvals(). It will be used to determine whether
        chromosome arrays should be created at their "real" size (i.e. using the
        actual sizes of human chromosomes) or using a fake, 'testing' size of
        100 bases. (Because in testing, I don't want to be dealing with
        chromosomes that are way bigger than they need to be.)
-->>The final line will determine where the module's output will be saved.
Any lines in the configfile beginning with # are ignored.

----------Module Output----------
The output of the module will be five files, which are described in the
docstring for outputpvals(). These files are saved in a dir called
"categorizeSNPsII_outputpvals_<date>_<GWAS>" in the location specified by the
"outputlocation" line of the config file. Each of these files contains
p-values (separated by newlines) for GWAS SNPs that fall within the
specified location (lincRNA, mRNA, both lincRNA and mRNA, nonintergenic, or
intergenic.) Using the files created by this module, the comparisons of
p-value distributions can be made in R, or using permutation testing with
the module enrichment.py"""

import cPickle
import sys
import copy
import time
import os
import os.path

import GTFparser_general
import FeatureClass_Small
import RefGene_parserII
import PseudogeneParser
import chromoarray

def do_snipification(theGWASpath):
    #Tested. I examined the original file of AcuteVSChronic_clean.txt and the
    #SNPs list generated from it at the beginning, middle, and end. The SNPs
    #list was created perfectly.
    """Return a Python list of SNPs (as lists) based on <GWASpath>.
    ----------Input---------
    <GWASpath> is the full filename and path to the annovar input file, e.g.
    /home/raba/GWAS_cleaned/Cardiogram_All/CardiogramGWAS_clean.txt
    This annovar input file can be created by the script cleanSNPdataset.sh
    The annovar input file is tab-delimited with seven columns:
          CHR     POS     POS     0   0   rsnumber    pvalue
    where CHR is the chromosome number given just as a number, and rsnumber is
    the rs number given just as a number (there are no "chr" or "rs" prefixes).
    
    ----------Output----------
    Each line of the GWAS file is turned into a list that represents a SNP:
        [rsnumber,chromosome,position,pvalue,category]
    so accessing [0] of a SNP will give the rsnumber, accessing [1] will give
    the chromosome, etc. Sorry that using five-member lists for SNPs isn't very
    readable; it's just faster than using a custom SNP object."""
    anvfile=open(theGWASpath,'r')
    snpslist=[]
    linecount=0
    for line in anvfile:
        linecount+=1
        lineaslist=line.rsplit("\t")
        chrm=int(lineaslist[0])
        pos=int(lineaslist[1])
        rsnum=int(lineaslist[5])
        pval=float(lineaslist[6])
        snpslist.append([rsnum,chrm,pos,pval])
    anvfile.close()
    assert linecount==len(snpslist), ("Error: "
            +"number of lines in annovarfile is "
            +`linecount`+" but number of SNPs in PyList is "+`len(snpslist)`)
    return snpslist

def _collectlocations(rnalist):
    """Return a list of RNA locations.
    ----------Function Input----------
    <rnalist> is a list of SmallFeature objects that represent RNAs
    
    ----------Function Output----------
    This function puts all of the RNA locations into a list and returns that
    list. Format of location list:
        [[chrm1],[chrm2],...,[chrm23,[chrm24]]
    Where a given chromosome consists of [int,int],[int,int],...,[int,int],
    that is, a bunch of integer intervals representing the locations of the
    RNAs on that chromosome.
    Note: chromosome X = chromosome 23 {index 22 of location list};
    chromosome Y = chromosome 24 {index 23 of location list}.
    Only autosome will be used later so X and Y are actually irrelevant."""
    for chromo in rnalist:
        GTFparser_general.giveranger(chromo)
    locations=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
        [],[],[]]
    #Careful: index 0 corresponds to chromosome 1
    chrcounter=0
    while chrcounter<22:
        for RNA in rnalist[chrcounter]:
            locations[chrcounter].append(RNA.ranger)
        chrcounter+=1
    return locations


def outputpvals(GWASpath,lincpath,mrnapath,nonintpath,bases,configpath,
                whether_to_use_real_sizes,outputloc):
    #See file "11-25-13_testing_the_module_categorizeSNPsII" for description
    #of how this function was tested
    """Create five files containing p-values.
    ----------Function Input----------
    The input is read out of the specified configuration file when the module
    is called.
    
    ----------Function Implementation----------
    This function makes use of arrays representing chromosomes to classify SNPs
    into five categories. If the arrays have not been created already,
    they will be created. One array is created per autosome. Also, a chromosome
    can have multiple arrays, if this function is called with different
    <bases> parameters. (see NOTE_8213)
    
    ----------Function Output----------
    This function creates five text files of p-values for the GWAS specified:
    prefix=<GWAS>_<date>_ext<bases>
    <prefix>_linconly.txt       : p-values for SNPs within a lincRNA only
    <prefix>_mRNAonly.txt       : p-values for SNPs  within an mRNA only
    <prefix>_both_linc+mRNA.txt : p-values for SNPs within both a lincRNA
                                           and an mRNA. 
    <prefix>_noninter.txt       : p-values for SNPs not in a lincRNA or mRNA,
                                           but ARE inside of "nonintergenic"
                                           regions defined by the file
                                           created by the function
                                           notintergenic() (see below).
    <prefix>_intergenic.txt     : p-values for all remaining SNPs (not in mRNA
                                           and/or lincRNA, not in
                                           "nonintergenic" regions)
    The ext<bases> in the output file names indicates by how many bases the
    size of the RNA was extended in either direction. A custom array must be
    created for each bases parameter."""
    #Extract simple GWAS name from path, and open the files that will have
    #pvals stored in them
    pathaslist=GWASpath.rsplit("/") #should the separator ever be "\\"
    GWASfilename=pathaslist[-1] 
    GWAS=GWASfilename[:-10] #WARNING: THIS ASSUMES THE FORMAT OF THE FILENAME
    #IS <GWAS>_clean.txt
    print "sorting GWAS SNPs by chrom"
    SNPs=do_snipification(GWASpath)
    SNPsbychrom=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
        [],[],[],[]]
    for SNP in SNPs:
        SNPsbychrom[(SNP[1])-1].append(SNP)
    SNPsbychrom=SNPsbychrom[:22] #Remove X and Y chromosomes so that list
    #index isn't out of range for nonintplaces in while loop below
    #Assign every SNP to a category and write to results files, which are all
    #saved into subdir of the GWAS dir
    prettydate=(`time.localtime().tm_mon`+"-"+`time.localtime().tm_mday`
                +"-"+`time.localtime().tm_year`)
    #creates new dir for results files
    _ensure_dir(outputloc+"/categorizeSNPsII_outputpvals_"+prettydate+"_"+GWAS)
    prefix=(outputloc+"/categorizeSNPsII_outputpvals_"+prettydate+"_"
            +GWAS+"/"+GWAS+"_"+prettydate+"_ext"+`bases`)
    fmrnas=open(prefix+"_mRNAonly.txt",'w')
    fintergenic=open(prefix+"_intergenic.txt",'w')
    flincs=open(prefix+"_linconly.txt",'w')
    fnonint=open(prefix+"_noninter.txt",'w')
    fboth=open(prefix+"_both_linc+mRNA.txt",'w')
    print "assigning SNPs to categories and writing to files"
    overallcount=0; lccount=0; mRcount=0
    bothcount=0; intcount=0; nonintcount=0
    for chrmdex in range(22): #ignore X and Y which have been removed anyway
        print "chromosome "+`chrmdex+1`+" is being worked on"
        chrarray=_loadchrarray(chrmdex,bases,lincpath,mrnapath,nonintpath,
                               configpath,whether_to_use_real_sizes)
        snpstoclassify=SNPsbychrom[chrmdex]
        for snip in snpstoclassify:
            overallcount+=1
            category=chrarray[snip[2]-1]#snip[2] is the snip's coordinate;
            #look up classification of the relevant base
            #have to subtract one because the snip's coordinate is one-based but
            #Python arrays have a "zeroth" index
            if category=="i":#Write to file based on category
                fintergenic.write(`snip[3]`+"\n")#snip[3] is pvalue
                intcount+=1
            elif category=="m":
                fmrnas.write(`snip[3]`+"\n")
                mRcount+=1
            elif category=="l":
                flincs.write(`snip[3]`+"\n")
                lccount+=1
            elif category=="n":
                fnonint.write(`snip[3]`+"\n")
                nonintcount+=1
            elif category=="b":
                fboth.write(`snip[3]`+"\n")
                bothcount+=1
            else:
                assert False, ("this should never happen: '"+category
                               +"' is not a valid category")
    fmrnas.close()
    fintergenic.close()
    flincs.close()
    fnonint.close()
    fboth.close()
    logg=open(prefix+"_LOGFILE_for_categorizeSNPs.txt",'w')#logfile has counts
    logg.write("overallcount: "+`overallcount`+
               "\nmRNA only: "+`mRcount`+
               "\nlincRNA only: "+`lccount`+
               "\nboth lincRNA and mRNA: "+`bothcount`+
               "\nintergenic: "+`intcount`+
               "\nnonintergenic: "+`nonintcount`)
    subcountsum=mRcount+lccount+bothcount+intcount+nonintcount
    assert overallcount==subcountsum, ("Error: overallcount was "
            +`overallcount`+" but the sum of the subcounts was "+`subcountsum`)

def _loadchrarray(chrmdex1,bases1,lincpath1,mrnapath1,nonintpath1,configpath1,
                  whether_to_use_real_sizes1):
    """Return the chrarray called "chromosome<chrmdex1+1>_ext<bases1>_array.bin"
    (for example, "chromosome1_ext0_array.bin").
    The chrarray should be located in the same directory as the configuration file.
    If the chrarray file already exists, it is loaded and returned.
    If the chrarray file does not exist, it is created using
    chromoarray.create_22chrom_arrays() for the
    correct chrmdex, and saved in the same directory as the configuration file.
    The creation of the chrarray file requires the nonintergenic regions to be
    defined. If they have been defined, they are loaded and used. If they have
    not been defined, then they are defined before being used."""
    confFilepathaslist=configpath1.rsplit("/") #should separator ever be "\\"?)
    confFiledirpath="/".join(confFilepathaslist[:-1])
    try:
        chrarray=cPickle.load(open(confFiledirpath+"/chromosome"
                +`chrmdex1+1`+"_ext"+`bases1`+"_array.bin",'rb'))
                #load the pre-created array for the appropriate chromosome
                #and bases1 parameter
    except IOError as myerr:
        if "No such file or directory" in myerr.strerror: #if you haven't
            #created the chromosome arrays yet, make them now
            lincRNAslist=cPickle.load(open(lincpath1,'rb'))
            linkLocations=_collectlocations(lincRNAslist)
            for chromzome in linkLocations:
                _clean_up_locations(chromzome) #If you don't clean up the
                #locations, then you have overlapping intervals, and it's just
                #cleaner to NOT have overlapping intervals
            mRNAslist=GTFparser_general.sortSmallFeatsbychrom(
                RefGene_parserII.returnRefGenelist(mrnapath1,["NM"]))
            mrnaLocations=_collectlocations(mRNAslist)
            for chromzome2 in mrnaLocations:
                _clean_up_locations(chromzome2)
            try: #similar setup for the nonintLocations file: create it if it
                #doesn't exist
               nonintLocations=cPickle.load(open(nonintpath1,'rb'))
            except IOError as myerr2:
                if "No such file or directory" in myerr2.strerror:
                    notintergenic(configpath1,nonintpath1)#name the file
                    #whatever is specified
                    #e.g. categorizeSNPs_NOTintergenic_locslist_<suffix>.bin
                    nonintLocations=cPickle.load(open(nonintpath1,'rb'))
                else:
                    raise IOError
            chromoarray.create_22chrom_arrays(mrnaLocations,linkLocations,
                                nonintLocations,bases1,[chrmdex1+1],
                                whether_to_use_real_sizes1,
                                dirtosavein=confFiledirpath)
            chrarray=cPickle.load(open(confFiledirpath+"/chromosome"
                            +`chrmdex1+1`+"_ext"+`bases1`+"_array.bin",'rb'))
        else:
            raise IOError #any IOError BESIDES not having created the chromosome
            #arrays will get reported
    return chrarray


def _ensure_dir(f):
    """Create directory <f> if it does not exist already.
    <f> is the complete path to a dir"""
    if not os.path.exists(f):
        os.makedirs(f)

def notintergenic(pathtoconfig,filename):
    #Tested in "SECOND TEST" as described in "10-15-13 testing categorizeSNPs
    #with fake files.docx"
    """Create a pickled Python list that describes 'notintergenic regions'
    
    ----------Function Output----------
    Creates a nested list of locationsin the same format as the output of
    _collectlocations() for regions that are going to be excluded from
    "intergenic" classification. The nested list is saved as a pickled Python
    list in the file "filename" in the same directory as the configuration file
    (specified by <pathtoconfig>) The choice of what to include in
    "notintergenic" is based on Hangauer paper. Within each chromosome the
    locations have been collapsed such that the minimum number of intervals
    is reported (to help speed up the process of determining whether a SNP falls
    in a 'notintergenic' region.)
    Regions included as "notintergenic":
        RefSeq NR genes (RNA) and XR genes (RNA predicted model)
        Gencode: everything in version 18
        Hangauer's extended protein coding gene structures derived from RNA-seq
        data (downloaded as S1)
        Pseudogenes from Yale
    Note: This is based on what Hangauer et al did but with a few differences,
    including: The "mappability" of the genome is not considered here
    (Hangauer only considered mappable regions); H-Invitational data is not
    used; alternative and extended 5' and 3' UTRs from UTRdb are not used;
    everything in Gencode was used instead of just "genes"
    Note: "the default human gene set in the Ensembl browser is therefore also
    the current version of GENCODE." Which is why it's weird to me that
    Hangauer mentions Gencode and Ensembl separately...?"""
    confFile=open(pathtoconfig,'r')
    confFilepathaslist=pathtoconfig.rsplit("/") #should separator ever be "\\"?
    confFiledirpath="/".join(confFilepathaslist[:-1])
    mRNApath=""; gencodepath=""; hang1path=""; pseudpath=""
    for line in confFile:
        if line[0]!="#":
            lineaslist=line.rsplit(" ")
            if lineaslist[0]=="mRNApath:":
                mRNApath=lineaslist[1]
                mRNApath=stripnewlineEnd(mRNApath)
            elif lineaslist[0]=="gencodepath:":
                gencodepath=lineaslist[1]
                gencodepath=stripnewlineEnd(gencodepath)
            elif lineaslist[0]=="hang1path:":
                hang1path=lineaslist[1]
                hang1path=stripnewlineEnd(hang1path)
            elif lineaslist[0]=="pseudpath:":
                pseudpath=lineaslist[1]
                pseudpath=stripnewlineEnd(pseudpath)
            else:
                pass
    if mRNApath=="" or gencodepath=="" or hang1path=="" or pseudpath=="":
        assert False, ("There was an error"
            +" parsing the configuration file: \n\tmRNA path was "
            +`mRNApath`+"\n\tgencodepath was "+`gencodepath`
            +"\n\thang1path was "+`hang1path`+"\n\tpseudpath was "+`pseudpath`)
    #Obtain all RefSeq NR and XR genes:
    refseqz=RefGene_parserII.returnRefGenelist(mRNApath,["NR","XR"])
    #Take everything from Gencode
    gencodez=GTFparser_general.makeSmallfeatures_fromGTF(gencodepath,True,"ALL")
    GTFparser_general.giveranger(gencodez)
    #Use Hangauer's extended protein coding gene structures
    assert "19" in hang1path, ("Error: file name did not contain 19--"
                    +"use the Hangauer S1 file that's been converted to hg19!")
    h1file=open(hang1path,'r')
    h1list=[]
    h1fails=open(confFiledirpath
                 +"/h1fails_from_notintergenic_in_categorizeSNPs.txt",'w')
    for line in h1file:
        try:
            lineaslist=line.rsplit("\t")
            chrom=int(lineaslist[0][3:]) #chop off the chr prefix
            startpoz=int(lineaslist[1])
            stoppos=int(lineaslist[2])
            h1list.append(FeatureClass_Small.SmallFeature("",".",chrom,
                                [startpoz],[stoppos],[1],[startpoz,stoppos]))
        except:
            h1fails.write(line) #for example, if the chromosome was 16_random,
            #casting to an int would fail, so the line would not contribute to
            #the "nonintergenic" regions and would get written to h1fails file
    h1fails.close()
    #Get Yale pseudogenes
    pseuds=PseudogeneParser.return_pseudogenes(pseudpath)
    #Extract the locations from all, and make unified location list by chrom
    allfeats=refseqz+gencodez+h1list+pseuds #TO DO add in Ensembl
    locilist=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
        [],[],[]]
    for feat in allfeats:
        locilist[(feat.chromosome)-1].append(feat.ranger)
    #Now make simplest possible list of positions for each chromosome that
    #covers all the necessary bases
    chrcounter=1
    for chromo in locilist:
        print "cleaning up the locations for chromosome "+`chrcounter`
        _clean_up_locations(chromo)
        chrcounter+=1
    print "pickling the cleaned locations list"
    cPickle.dump(locilist,open(filename,'wb'),2)
    print "All Done :D"

def stripnewlineEnd(word):
    """Strip away annoying newlines to help in parsing the configuration file"""
    #NOTE: better way of doing this is with string.replace("\n","")
    if word.endswith("\n"):
        word=word[:-1]
    return word

def _clean_up_locations(localist):
    #A nearly-identical version of this function was tested in
    #definelincs_test.test9()
    """Merge locations that overlap by one or more bases, to transform a
    list of overlapping locations into a list of distinct non-overlapping
    locations that cover exactly the same regions of the genome.
    localist is a two-dimensional list e.g. [[12,34],[45,62],30,70]] where each
    sub-list has two members, a start and a stop position.
    A function called _clean_up_locus_locations() in definelincs.py has exactly
    the same logic as this function and passed testing by
    definelincs_test.test9()"""
    pos1=0
    pos2=1
    while pos1<len(localist)-1 and pos2<len(localist):
        locusA=localist[pos1]
        locusB=localist[pos2]
        startA=locusA[0]; stopA=locusA[1]
        startB=locusB[0]; stopB=locusB[1]
        if startA <= startB <= stopA:
            if stopB >= stopA: #loci offset from each other
                locusB[0]=startA
                localist.remove(locusA)
                pos2=pos1+1
            else: #locusB is "engulfed by" locusA
                assert stopB <= stopA
                locusB[0]=startA
                locusB[1]=stopA
                localist.remove(locusA)
                pos2=pos1+1
        elif startB <= startA <= stopB:
            if stopB >= stopA: #locusA is "engulfed by" locusB
                localist.remove(locusA)
                pos2=pos1+1
            else: #loci offset from each other
                assert stopB <=stopA
                locusB[1]=locusA[1]
                localist.remove(locusA)
                pos2=pos1+1
        else:#They don't overlap and were not combined
            #ChangePositionAfterNoCombo
            length=len(localist)
            if pos2 < length-1 and pos1 < length-2:
                pos2=pos2+1
            elif pos2==length-1 and pos1<length-2:
                pos1=pos1+1
                pos2=pos1+1
            else:
                assert pos2==length-1 and pos1==length-2, ("ERROR: "
                        +"pos2 is actually "+`pos2`+" and pos1 is "+`pos1`)
                pos1=length #stop the while loop
    localist.sort(key=lambda x: x[0]) #sort the loci by start position
    #Note: this function does not return the localist; it modifies the
    #localist in place

if __name__=='__main__':
    confpath=sys.argv[1]
    conffile=open(confpath,'r')
    myGWASstudies=[]; mylincpath=""; mymrnapath=""; mybases=-1
    mynonintpath=""; myuserealsizes=""; myoutputlocation=""
    for line in conffile:
        if line[0]!="#":
            lineaslist=line.rsplit(" ")
            if lineaslist[0]=="GWAS:":
                myGWAZ=stripnewlineEnd(lineaslist[1])
                myGWASstudies.append(myGWAZ)
            elif lineaslist[0]=="lincRNApath:":
                mylincpath=stripnewlineEnd(lineaslist[1])
            elif lineaslist[0]=="mRNApath:":
                mymrnapath=stripnewlineEnd(lineaslist[1])
            elif lineaslist[0]=="nonintpath:":
                mynonintpath=stripnewlineEnd(lineaslist[1])
            elif lineaslist[0]=="bases:":
                mybases=int(stripnewlineEnd(lineaslist[1]))
            elif lineaslist[0]=="userealsizes:":
                if stripnewlineEnd(lineaslist[1])=="True":
                    myuserealsizes=True
                elif stripnewlineEnd(lineaslist[1])=="False":
                    myuserealsizes=False
                else:
                    raise stupiderror(Exception)
            elif lineaslist[0]=="outputlocation:":
                myoutputlocation=stripnewlineEnd(lineaslist[1])
            else:
                pass
    if (myGWASstudies==[] or mylincpath=="" or mymrnapath=="" or mybases==-1 or
        mynonintpath=="" or myuserealsizes=="" or myoutputlocation==""):
        assert False, ("There was an error"
                       +" in reading the config file. \nlincpath was "
                       +`mylincpath`+
                       "; \nmrnapath was "+`mymrnapath`+"; \nnonintpath was "
                       +`mynonintpath`+
                       "; \nbases was "+`mybases`+"; \nuserealsizes was "
                       +`myuserealsizes`+
                       "; \noutputlocation was "+`myoutputlocation`)
    assert type(mybases) is int and 0<=mybases<=20000, ("Error:"
            +" mybases must be an int, and 0 <= mybases <= 20000. "
            +"You chose mybases="+`mybases`)
    for myGWASpath in myGWASstudies:
        print "running outputpvals for "+myGWASpath+" with extbases "+`mybases`
        outputpvals(myGWASpath,mylincpath,mymrnapath,mynonintpath,mybases,
                    confpath,myuserealsizes,myoutputlocation)
        print "done running outputpvals"

