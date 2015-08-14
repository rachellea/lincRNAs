#Rachel Ballantyne 2014
#analyzeLD.py

"""Report the linkage disequilibrium between top GWA SNPs in sets of features.

#===============================================================================
           PARAMETERS, CONFIGURATION FILE, AND DESCRIPTION OF ANALYSES:
#===============================================================================
This module is called with one argument: the path to a configuration file.
e.g.
python analyzeLD.py /path/to/config.txt

Here is an example configuration file:

   -----------------------------------------------------
   | #########################                         |
   | # analyzeLD config file #                         |
   | #########################                         |
   | RunFocusPeriphAnalysis: YES                       |
   | RunTopRegionSNPAnalysis: YES                      |
   | RunSpecificPairingAnalysis: YES                   |
   |                                                   |
   | FocusFeatures: /path/to/file1                     |
   | PeripheralFeatures: /path/to/file2                |
   | SpecificPairing: /path/to/file3                   |
   |                                                   |
   | GWAS: GLGC_HDL: /path/to/GLGC_HDL_clean.txt       |
   | GWAS: GLGC_TC: /path/to/GLGC_TC_clean.txt         |
   | GWAS: DIAGRAM_T2D: /path/to/DIAGRAM_T2D_clean.txt |
   | LDfiles: /path/to/LDchr#.txt                      |
   | Keyword: mykeyword                                |
   | OutputDir: /where/output/will/go/                 |
   | ExtBases: 500000                                  |
   -----------------------------------------------------

Note that there must be a space after every colon in the config file.
All fields are required to have a value. If a field is not applicable to you
(see below), write NONE as its value.
You must have at least one GWAS field; there is no upper limit on the number of
GWAS fields.
Lines beginning with # are ignored.

The first three fields, RunFocusPeriphAnalysis, RunTopRegionSNPAnalysis, and
RunSpecificPairingAnalysis, will determine which analyses the module will run.
Each of these fields can be YES or NO. If all these three fields are
NO, then the module will do nothing.
Note: 
    > If RunFocusPeriphAnalysis is YES, then the fields FocusFeatures and
      PeripheralFeatures must have valid paths. SpecificPairing can be NONE.
    > If RunTopRegionSNPAnalysis is YES, then the field FocusFeatures must have
      a valid path. PeripheralFeatures and SpecificPairing can both be NONE.
    > If RunSpecificPairingAnalysis is YES, then the field SpecificPairing must
      have a valid path. FocusFeatures and PeripheralFeatures can both be NONE.

<FocusFeatures> indicates the path to a file of features which will be the focus
of the FocusPeriphAnalysis and RunTopRegionSNPAnalysis. The file must be
tab-delimited with five columns: feature name, GWAS name, chromosome as an int,
start position, and stop position. The start/stop positions should be in a
one-based inclusive coordinate system.
You may include the same name more than once in the file with a different GWAS
but the same position info throughout. Do not include the same name with
different positions, or you may get inconsistent results. If you want to look at
isoforms of a gene, please give each isoform a distinct name.

<PeripheralFeatures> is the path to a file with the same tab-delimited
five-column format as the FocusFeatures file. It describes the
PeripheralFeatures.

<SpecificPairing> is the path to a tab-delimited file with the following
ten columns: the name, GWAS, chr, start, and stop of feature 1, followed by the
name, GWAS, chr, start, and stop of feature 2. It can have as many rows
as you would like, and they can all refer to different pairs of features.
If some features listed in this file are also FocusFeatures described in a
FocusFeatures file, then describe those features in the first five columns
(rather than the last five columns) of each row.

If RunFocusPeriphAnalysis is YES, then the following analysis ("FocusPeriphA")
will be run:
    The region around each FocusFeature will be searched to
    determine whether any of the PeripheralFeatures are present. If any
    PeripheralFeatures are in the region, the LD between the top GWAS SNP in the
    FocusFeature and the top GWAS SNP in the PeripheralFeature will be
    calculated. The "top GWAS SNP" in a feature means the SNP with the lowest
    GWAS p-value in that feature. The GWAS p-values used will be for the GWAS
    of the FocusFeature (regardless of what the GWAS of the PeripheralFeature
    is. The GWAS of the PeripheralFeature will still be noted in the output
    file, though.)
    For example, if FocusFeatures are lincRNAs and PeripheralFeatures are GWAS
    candidate genes, then this module will search for candidate genes in the
    region around lincRNAs, and calculate the relevant LD. PeripheralFeatures
    could also be genome-wide Bonferroni-significant SNPs, in which case this
    module will search for genome-wide significant SNPs in the region of each
    FocusFeature and calculate LD between that SNP and the FocusFeature.
    Note that if the PeripheralFeature is a SNP, then its start position will
    equal its end position.

If RunTopRegionSNPAnalysis is YES, then the following analysis ("TopRegSNPA")
will be run:
    The module will search the region around each FocusFeature to determine what
    the top GWAS SNP in that region is (based on the p-values from the GWAS of
    the FocusFeature.) If the top GWAS SNP is inside of the FocusFeature, the
    FocusFeature will be assigned a "status" of 1. If the top GWAS SNP is
    outside of the FocusFeature, then LD will be calculated between that
    external top SNP and whatever the SNP with the lowest p-value within the
    FocusFeature is.
    If this LD is <=0.3, the FocusFeature will receive a status of 2; if this
    LD is 0.3 < R2 < 0.8, status 3; if LD is >=0.8, status 4. These status codes
    will be reported in the analysis output.

If RunSpecificPairingAnalysis is YES, then the following analysis ("SpecPairA")
will be run:
    The pairs of features described in the ten-column file at the path indicated
    by <SpecificPairing> will be analyzed. Specifically, the LD will be
    calculated between the top GWAS SNP in feature 1 versus the top GWAS SNP in
    feature 2, for every row of the file.
    If you want to calculate feature 1 vs feature 2, and feature 1 vs feature 3,
    then include two lines in the file, one with feature 1 and feature 2 info,
    and the other with feature 1 and feature 3 info. This file can include as
    many lines and feature combinations as you would like. Use this analysis
    instead of the FocusPeriphAnalysis if you care about LD only between
    very specific pairings of features, rather than LD between sets of features
    that are in close proximity. An example application of the
    RunSpecificPairingAnalysis is to calculate LD between the top SNP in a gene
    and the top SNP in a neighboring gene.

The following configuration file fields are always required to have valid
values, regardless of which set of analyses you are running:

The <GWAS> field indicates what GWA studies are being considered. You must have
one GWAS field for every GWAS that is mentioned in your FocusFeatures and
PeripheralFeatures files. Immediately after the <GWAS:> key, record the
GWAS name you want to use (it cannot include spaces). Then put a colon and a
space followed by the path to a file describing the results of that GWAS.
This results file should be tab-delimited with 7 columns: chromosome as an int,
SNP coordinate, SNP coordinate, 0, 0, SNP rs number or other identifier,
and SNP p-value (this is also the input format for the software ANNOVAR). 
Make sure the GWAS names you use in the configuration file match the GWAS names
you use when describing the GWAS of the FocusFeatures or PeripheralFeatures.
I.e. if you referred to the GLGC total cholesterol GWAS as "GLGC_TC" in your
FocusFeatures and PeripheralFeatures files, then you should use "GLGC_TC" as
the name for it in the configuration file-- GWAS: GLGC_TC: path/to/file. If you
used "GLGC_TotalChol" instead, that would generate an error. 

<LDfiles> indicates where the PLINK output files storing LD information can
be found. You must have a separate file for every chromosome that you want
to consider. Name all the files after the same pattern, e.g. chr2_Eur.ld,
chr10_Eur.ld, chr17_Eur.ld
Then in the <LDfiles> field record the path to one of the LD files and replace
the chromosome number with the # symbol, e.g. /path/to/LD/chr#_Eur.ld
The only requirements for the file naming are that the pattern of file naming
be the same for all chromosomes, and the chromosome number be included
in the file name.
The LD files should be in the PLINK output format: tab-delimited with the
following six columns: CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B, R2. The CHR
columns indicate chromosome as an integer; the BP columns indicate SNP
coordinate as an integer; the SNP columns include a SNP identifier (e.g. an
rs number or Thousand Genomes chr##:#### identifier); the R2 column indicates
the LD between SNP_A and SNP_B, as a float 0 <= R2 <= 1. The file should be
redundant in that if rs12345 appears as SNP_A against rs98760 as SNP_B once,
then this R2 should also be recorded elsewhere with rs12345 as SNP_B and
rs98760 as SNP_A (this redundancy is the default behavior of PLINK when
calculating LD files like this.)
The example LD files provided with this module was generated using plink 
and 1000 Genomes European data. Only SNP pairs with Rsquared > 0.3 are included
in the file (otherwise the file size would be un-usably huge). Here is an
example of the plink command that was used:
    plink --noweb --file /project/mingyaolab/Ped_Map/EUR22 --r2 --ld-snp-list
    /home/raba/European_LD/chr22rs.id --ld-window-kb 1000 --ld-window 99999
    --ld-window-r2 0.3 --out chr22_Jan15

<Keyword> is a phrase without spaces that will be included in the file names
of module output files.

<OutputDir> is the path to a directory where output files will be stored.

<ExtBases> is a number, in bases, that dictates the region size around every
feature in FocusFeatures. Only the PeripheralFeatures that fall inside of
the region around a FocusFeature will be used. ExtBases will be used to
extend the FocusFeature by the specified number of bases in both directions
(i.e. if ExtBases is 500 then for a FocusFeature located at coordinates
1200-1800, the region for this FocusFeature will include 700-2300.)
Warning: if ExtBases is chosen to be too large, this module may take a very
long time to run, because there will be many PeripheralFeatures in the region
of each FocusFeature.

#===============================================================================
                             OUTPUT FILE FORMATS:
#===============================================================================
The FocusPeriphAnalysis (FocusPeriphA) and SpecificPairingAnalysis (SpecPairA)
each have two output files: "compact" (C) and "non-compact" (NC).
They are named:
    > analyzeLD_C_FocusPeriphA_<keyword>_<Month-Day>.txt
    > analyzeLD_NC_FocusPeriphA_<keyword>_<Month-Day>.txt
    > analyzeLD_C_SpecPairA_<keyword>_<Month-Day>.txt
    > analyzeLD_NC_SpecPairA_<keyword>_<Month-Day>.txt
The TopRegionSNPAnalysis (TopRegSNPA) only has a "compact" output file:
    > analyzeLD_C_TopRegSNPA_<keyword>_<Month-Day>.txt
since there is guaranteed to be only one 'other feature' in each Region and
there is thus no difference between the C and NC output file formats.

The "compact" and "non-compact" output files contain the same information
except that the "compact" output file includes GWAScategory and
LDcategory columns and the "non-compact" file does not.
The "compact" file describes one Region per line, and the "non-compact"
file has multiple lines for each Region (one line for every 'candidate'
in that region.)

Compact format: there is one line for every FocusFeature. Columns:
    > focusfeatname
    > focusfeatgwas
    > focusfeatposition: chr#:start-stop format
    > focusfeatTopSNPname
    > focusfeatTopSNPpval
    > OtherFeatsCount: the total count of other features ('Candidates',
      see implementation) within the FocusFeature's region
    > TotalOtherFeatNames: the names of all PeripheralFeatures within the
      FocusFeature's region, separated by asterisks: NAME1*NAME2*NAME3
  The following columns all have the same format: the name of an 'other
  feature' (PeripheralFeature, regional top SNP, whatever) attached to
  some additional info about that feature by an equals sign:
    > OtherFeatGWASs: e.g. GENE1=GWAS1*GENE2=GWAS2*GENE3=GWAS3
    > OtherFeatPositions:
      e.g. GENE1=chr5:1234-6789*GENE2=chr6:100-200*GENE3=chr19:80300-80600
    > OtherFeatTopSNPname:
      e.g. GENE1=rs495867*GENE2=rs10293854*GENE3=rs99488456 when the SNP has
      an rs number. If the SNPs have 1000 Genomes identifiers:
      GENE1=chr15:9847556*GENE2=chr21:9298438*GENE3=chr8:2948384
    > OtherFeatTopSNPpval: e.g. GENE1=0.002*GENE2=0.98*GENE3=0.000002354
    > Distance: this is the distance between the focusfeatTopSNP and the
      OtherFeatTopSNP. Format: GENE1=1002*GENE2=9385*GENE3=9823
    > LD: this is the LD (R squared) between the focusfeatTopSNP and the
      OtherFeatTopSNP. Format: GENE1=0.47*GENE2=0.33*GENE3=0.55
    > GWAScategory: This is either 1 or 2. A GWAScategory of 1 indicates
      that the FocusFeature contains the SNP with the lowest p-value out of
      all the top SNPs in each candidate in the region, for the candidates
      relevant to that particular analysis. A GWAScategory of 2 indicates that
      the FocusFeature does NOT contain this SNP.
      For the TopRegSNP analysis, a GWAScategory of 1 means that the
      FocusFeature contains the top SNP in the entire Region. For the
      FocusPeriph analysis and the SpecificPairing analysis, a GWAScategory of 1
      only means that the FocusFeature's top SNP has a lower p-value than the
      top SNPs of the other candidates for that analysis.
    > LDcategory: This is either 1, 2, 3, or 9.
        - LDcategory 1: considering all the candidates in the Region for a
          particular analysis, the highest LD between the top SNP in a
          candidate and the top SNP in a FocusFeature is R2 <= 0.3
        - LDcategory 2: 0.3 < R2 < 0.8
        - LDcategory 3: R2 >= 0.8
        - LDcategory 9: no R2 was found between the top SNP in any candidate and
          the top SNP in the FocusFeature.
        - Note: this categorization could be misleading if your LD input files
          are incomplete. For example, consider a FocusFeature with two
          candidates in the region. If this FocusFeature has high LD with
          Candidate X, but this LD data is not included in the LD input file,
          and the FocusFeature has low LD with Candidate Y, which IS included in
          the input file, then the FocusFeature will be categorized as having
          low LD with candidates in the region, EVEN THOUGH it really has high
          LD with Candidate X! So if you know you have a lot of missing data
          in your LD input files, be careful about relying on LD category for
          filtering your results.

Non-compact format: here is one line for every FocusFeature/other feature
pairing. Columns:
    > focusfeatname
    > focusfeatgwas
    > focusfeatposition
    > focusfeatTopSNPname
    > focusfeatTopSNPpval
    > otherfeatName
    > otherfeatGWAS
    > otherfeatposition
    > otherfeatTopSNPname
    > otherfeatTopSNPpval
    > distance (between top SNP in FocusFeature and top SNP in other feat)
    > LD (between top SNP in FocusFeature and top SNP in other feat)"""

import copy
import cPickle
import os.path
import shelve
import sys
import time

#===============================================================================
#----------FUNCTIONS------------------------------------------------------------
#===============================================================================
def returnCandidates(candidatespath):
    """Return the PeripheralFeatures ('Candidates') in file at <candidatespath>.
    PARAMETERS:
    <candidatespath> is the path to the file for PeripheralFeatures. See
        description of PeripheralFeatures file in module docstring.
    OUTPUT: 
    This function returns a list of lists of dictionaries, where each of the
    22 sub-lists represents an autosome, and each dictionary describes a
    candidate gene on that autosome. See Region Class for more info about
    the dictionary format. Do not include the same candidate more than once."""
    print "Initializing PeripheralFeatures from file at "+candidatespath
    candidatesfile=open(candidatespath,'r')
    allcandidates=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[]]
    for line in candidatesfile:
        lineaslist=line.rsplit()
        allcandidates[int(lineaslist[2])-1].append({"name":lineaslist[0],
                    "GWAS":lineaslist[1],"chrom":int(lineaslist[2]),
                    "start":int(lineaslist[3]),"stop":int(lineaslist[4])})
    print "Done initializing PeripheralFeatures"
    return allcandidates

def initializeRegionobjects(whichanalyses,gwaspathsdict):
    """Initialize the Region objects for the files and analysis specified.
    PARAMETERS:
    <whichanalysie> is a dictionary indicating the file(s) that should be used
        for creating the Region objects. The keys are analysis types and
        the values are paths to files. See module execution for details
    <GWASpathsdict> is a dictionary where the keys are simple GWAS names and the
        values are paths to GWAS files. See description of GWAS files in the
        module docstring.
    OUTPUT:
    This function returns a list of Region objects based on the file
    specified, where the FocusFeature of each Region object is initialized.
    Note that the Region objects will also include candidate information, but
    only for the Regions relevant to a 'SpecificPairingAnalysis' -- otherwise,
    addCandidatestoRegionobjects() will have to be used on this function's
    output to add the candidate information."""
    print "Initializing FocusFeatures"
    regobjects=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[]]
    for analysistype in whichanalyses: #analysistypes=keys in {whichanalyses}
        featsfile=open(whichanalyses[analysistype],'r')
        for line in featsfile:
            #Read what FocusFeature is being described
            lineaslist=line.rsplit()
            featname=lineaslist[0]
            featgwas=lineaslist[1]
            featchrom=int(lineaslist[2])
            featstart=int(lineaslist[3])
            featstop=int(lineaslist[4])
            featuresnps=getTopSNP(featchrom,featstart,featstop,
                                  gwaspathsdict[featgwas])
            newobject=Region({"name":featname,"chrom":featchrom,
                        "start":featstart,"stop":featstop,"SNPs":featuresnps,
                        "GWAS":featgwas},[],{},{},{})
            #Check whether a Region with the same FocusFeature already exists
            #so that if needed you can change what the variable 'newobject'
            #is pointing to
            alreadyexists=False
            for existingob in regobjects[featchrom-1]:
                if newobject.focusfeat==existingob.focusfeat:
                    newobject=existingob
                    alreadyexists=True
                    break
            #Update the trackanalysis property of the Region object
            if "self" not in newobject.trackanalyses:
                newobject.trackanalyses["self"]=[analysistype]
            else:
                if analysistype not in newobject.trackanalyses["self"]:
                    newobject.trackanalyses["self"].append(analysistype)
            #Add it if it doesn't already exist
            if (not(alreadyexists)):
                regobjects[featchrom-1].append(newobject)
            #If applicable, now parse the additional lines of SpecPairA files:
            if analysistype=="SpecPairA":
                newcnd={"name":lineaslist[5],"chrom":int(lineaslist[7]),
                        "start":int(lineaslist[8]),"stop":int(lineaslist[9]),
                        "GWAS":lineaslist[6]}
                if newcnd not in newobject.candidates:
                    newobject.candidates.append(newcnd)
                else:
                    for c in newobject.candidates:
                        if c==newcnd:
                            newcnd=c #make the variable 'newcnd' refer to the
                            #already-existing object
                            break
                nc=newcnd["name"]
                if nc not in newobject.candidateGWASpvals:
                    newobject.candidateGWASpvals[nc]=getTopSNP(newcnd["chrom"],
                        newcnd["start"],newcnd["stop"],gwaspathsdict[featgwas])
                if nc not in newobject.trackanalyses:
                    newobject.trackanalyses[nc]=(["SpecPairA"])
                else: 
                    if ("SpecPairA" not in newobject.trackanalyses[nc]):
                        newobject.trackanalyses[nc].append("SpecPairA")
        featsfile.close()
    return regobjects

def addCandidatestoRegionobjects(extbases,regobjects,gwaspathsdict,
                                 whichanalysis,candidatespath=""):
    """Assign 'Candidates' to each FocusFeature for the FocusPeriphA or
    TopRegSNPA analyses.
    PARAMETERS:
    <extbases> determines how far each FocusFeature will be extended in
        either direction to create the region. See description of
        ExtBases in module docstring.
    <regobjects> is the output of initializeRegionobjects()
    <gwaspathsdict> is a dictionary where the keys are simple GWAS names and the
        values are paths to GWAS files. See description of GWAS files in the
        module docstring.
    <whichanalysis> is either 'FocusPeriphA' or 'TopRegSNPA' and determines
        this function's behavior. (For 'SpecPairA' the candidates were
        already initialized in initializeRegionObjects().)
    <candidatespath> is the path to the file for PeripheralFeatures. This
        argument is only required if <whichanalysis> is 'FocusPeriphA'. See
        description of PeripheralFeatures file in module docstring.
    OUTPUT: This function modifies <regobjects> in places. It does not have an
    explcit return value.
    If <whichanalysis> is 'FocusPeriphA':
        For each Region object, this function looks through all the candidates
        (PeripheralFeatures) and adds them to the 'candidates' property
        of the Region object if they are within the region. Then it fills in the
        'candidateGWASpvals' on the basis of which candidates are in the region.
        The candidate must be completely within the region (both start and stop)
        in order to be included.
    If <whichanalysis> is 'TopRegSNPA':
        For each Region object, this function looks for all GWAS SNPs in each
        Region. If the top SNP in the Region is outside of the FocusFeature,
        then this SNP is added to the 'candidates' property."""
    if whichanalysis=="FocusPeriphA":
        allcandidates_bychrom=returnCandidates(candidatespath)
        #for each tblobject, determine if there are any candidates in its region
        print "Finding the PeripheralFeatures in the region of each FocusFeat"
        for index in xrange(22):
            for ob in regobjects[index]:
                regionstart=ob.focusfeat["start"]-extbases
                regionstop=ob.focusfeat["stop"]+extbases
                for can in allcandidates_bychrom[index]:
                    if ((regionstart <= can["start"] <= regionstop)
                        and (regionstart <= can["stop"] <= regionstop)):
                        if can not in ob.candidates: #only add it if it's
                            #not already there
                            ob.candidates.append(can)
                        #Now update trackanalysis property of the Region
                        cn=can["name"]
                        if cn not in ob.trackanalyses:
                            ob.trackanalyses[cn]=["FocusPeriphA"]
                        else:
                            if "FocusPeriphA" not in ob.trackanalyses[cn]:
                                ob.trackanalyses[cn].append("FocusPeriphA")
        #Now update the candidateGWASpvals property
        #you CANNOT update the 'trackanalyses' property here because you're
        #iterating over ALL candidates currently existing, which may include
        #some SpecPairA candidates that you do NOT want to also classify as
        #FocusPeriphA
        print "Finding the GWAS p-values for each PeripheralFeature"
        for index2 in xrange(22):
            for ob in regobjects[index2]:
                for candi in ob.candidates:
                    cn=candi["name"]
                    if cn not in ob.candidateGWASpvals:
                        ob.candidateGWASpvals[cn]=getTopSNP(
                            candi["chrom"],candi["start"],candi["stop"],
                            gwaspathsdict[ob.focusfeat["GWAS"]])
    if whichanalysis=="TopRegSNPA":
        print "Finding the top SNP in each Region"
        for index in xrange(22):
            for ob in regobjects[index]:
                topregSNP=getTopSNP(ob.focusfeat["chrom"],
                                    ob.focusfeat["start"]-extbases,
                                    ob.focusfeat["stop"]+extbases,
                                    gwaspathsdict[ob.focusfeat["GWAS"]])[0]
                if ob.focusfeat["SNPs"][0]!=topregSNP and topregSNP[0]!="rsNA":
                    newcnd={"name":topregSNP[0],"chrom":index+1,
                            "start":topregSNP[2],"stop":topregSNP[2],
                            "GWAS":ob.focusfeat["GWAS"]}
                    if newcnd not in ob.candidates:
                        ob.candidates.append(newcnd)
                    else: #newcnd is already there
                        for c in ob.candidates:
                            if c==newcnd:
                                newcnd=c #make the variable 'newcnd'
                                #refer to the already-existing object
                                break
                    cn=newcnd["name"]
                    if cn not in ob.candidateGWASpvals:
                        ob.candidateGWASpvals[cn]=[topregSNP]
                    if cn not in ob.trackanalyses:
                        ob.trackanalyses[cn]=["TopRegSNPA"]
                    else:
                        if "TopRegSNPA" not in ob.trackanalyses[cn]:
                            ob.trackanalyses[cn].append("TopRegSNPA")

def shelve_LD_files(myregions,LDlocation,chrom):
    """Create a shelve file for <chrom> containing LD info from <LDlocation>.
    PARAMETERS:
    <myregions> is a list of region objects that will require LD lookup.
    <LDlocation> is the path to an LD file, with the specific chromosome
        number replaced by #. See LDfiles description in module docstring.
    <chromnum> is an int for chromosome (chr1 = 1)
    NOTES: Shelve keys are SNP_A rs number and values are lists of strings
    where each string has the format "SNP_B_name R2" for every SNP_B that
    SNP_A is in LD with. The usefulness of this function's output relies
    entirely on having consistent SNP identifiers between the GWAS files and
    the PLINK LD output files, because the output will be searched for
    GWAS SNPs on the basis of rs number.
    WARNING: the LD files at <LDlocation> must have their entries grouped by
    SNP_A name (not necessarily alphanumerically sorted, but every entry
    related to SNP_A should appear consecutively). Producing output grouped
    by SNP_A name is the default behavior of PLINK, but if different software
    was used to acquire LD information, you must specifically make sure it is
    grouped by SNP_A name."""
    print "Creating LDshelve for LD lookup, chr"+`chrom+1`
    oldLDfile=open(LDlocation.replace("#",`chrom+1`),'r')
    LDshelve=shelve.open(LDlocation.replace("#",`chrom+1`)+"shelved.db",
                  protocol=2,writeback=False)
    previousrsnum=""
    previousval=[]
    for line in oldLDfile:
        lineaslist=line.rsplit()
        SNPAname=lineaslist[2]
        info=lineaslist[5]+" "+lineaslist[6] #SNP_B name, R2
        if SNPAname==previousrsnum:
            previousval.append(info)
        else: #new rs num
            LDshelve[previousrsnum]=previousval
            previousrsnum=SNPAname
            previousval=[info]
    #once done looping over file, need to write out the info for the very last
    #SNP_A recorded in the file: (otherwise it will never be written out)
    LDshelve[previousrsnum]=previousval
    oldLDfile.close()
    LDshelve.close()

def editRegions_re_LD(myregions,LDlocation):
    """Add LD information to Region objects.
    PARAMETERS:
    <myregions> is the list of Region objects that will have LD info added.
    <LDlocation> is the path to an LD file, with the specific chromosome
        number replaced by #. See LDfiles description in module docstring.
    NOTES: This function looks for SNPs in the LD files on the basis of
    location, NOT SNP identifier. It is thus extra important to make sure that
    all of your input files have a consistent coordinate system, especially
    the LD files and GWAS files, otherwise it will be impossible to properly
    match SNPs across files in the manner necessary to complete these analyses.
    OUTPUT: This function modifies <myregions> in place, to add LD information.
    There is no explicit return value."""
    #Figure out which chromosomes are relevant (don't want to go thru work
    #of making a shelve file for a chromosome unless there is at least one
    #Region object on that chromosome)
    relevant=[]
    for index in xrange(22):
        if len(myregions[index])>0:
            relevant.append(index)
    pairnotfound=0
    LDnotfound=0
    for index in relevant:
        chrm=myregions[index]
        if not os.path.exists(LDlocation.replace("#",`index+1`)+"shelved.db"):
            shelve_LD_files(myregions,LDlocation,index)
        LDshelve=shelve.open(LDlocation.replace("#",`index+1`)+"shelved.db",
                  protocol=2,writeback=False)
        print "Looking up LD for chr"+`index+1`
        for ob in chrm:
            assert len(ob.focusfeat["SNPs"])==1 #Note that the only SNP stored
            #in focusfeat["SNPs"] will be the SNP in the focusfeat with the
            #smallest GWAS p-value. See getTopSNP() which is the function used
            #to populate the "SNPs" value.
            focusfeatSNPname=ob.focusfeat["SNPs"][0][0]
            try:
                excerptfromfile=LDshelve[focusfeatSNPname]
                #excerptfromfile is a list of strings describing all the rows in
                #the LD file where SNP_A is the top SNP in the FocusFeat.
                #Format of strings: "SNP_B_name R2"
            except KeyError:
                #it's possible that PLINK didn't find any LD between the top
                #SNP in the linc and any other SNPs, in which case the
                #LDshelve will not contain a key of SNP_A coordinate
                excerptfromfile=[]
            #Now, look for the instances where SNP_B is the top SNP in one
            #of the candidates
            for candidate in ob.candidates:
                candidateSNPname=ob.candidateGWASpvals[candidate["name"]][0][0]
                for row in excerptfromfile:
                    rowaslist=row.rsplit()
                    if candidateSNPname==rowaslist[0]:
                        ob.candidateLD[candidate["name"]]=[(focusfeatSNPname,
                                candidateSNPname,float(rowaslist[1]))]
                    #tuple format (SNP_A rs, SNP_B rs, R-squared)
                if candidate["name"] not in ob.candidateLD:
                    ob.candidateLD[candidate["name"]]=[("rsNA","rsNA",2)]
                    if len(excerptfromfile)!=0:
                        #LD was found between focusfeatSNP and some other SNPs,
                        #but not with this candidate's top SNP in particular
                        pairnotfound+=1
                    else: #LD was not found for focusfeatSNP with any SNPs
                        LDnotfound+=1
            assert (len(set([c["name"] for c in ob.candidates]))
                    ==len(set(ob.candidateLD.keys()))), ("Error:"
                    +" names of ob.candidates were "
                    +`set([c["name"] for c in ob.candidates])`
                    +" but keys in ob.candidateLD were "
                    +`ob.candidateLD.keys()`)
                #use 'set' rather than comparing lengths directly because you
                #could have two candidates with the same name but different
                #internal properties, e.g. GENE1 in GWASA and GENE1 in GWASB,
                #in which case you would have two entries in ob.candidates
                #but only one entry in ob.candidateLD
    print ("Count: LD for top FocusFeature SNP was not found with "
            +"any SNPs: "+`LDnotfound`)
    print ("Count: LD for top FocusFeature SNP was found with some SNPs, "
            +"but not with a particular candidate gene's top SNP: "
            +`pairnotfound`)

def summarizeFocusPeriphAnalysis(regobjects):
    """Print summary stats about the number of candidates for Region objects."""
    cand0=0;cand1=0; cand2=0; cand3=0; cand4=0; candgt5=0
    for chrom in regobjects:
        for ob in chrom:
            if len(ob.candidates)==0:
                cand0+=1
            elif len(ob.candidates)==1:
                cand1+=1
            elif len(ob.candidates)==2:
                cand2+=1
            elif len(ob.candidates)==3:
                cand3+=1
            elif len(ob.candidates)==4:
                cand4+=1
            elif len(ob.candidates)>=5:
                candgt5+=1
    summarystring=("Summarizing table entries:\n\tNumber of features Y"
            +" having X candidate genes in region, printed as X: Y"
            +"\n\t0: "+`cand0`+"\n\t1: "+`cand1`+"\n\t2: "+`cand2`+"\n\t3: "
            +`cand3`+"\n\t4: "+`cand4`+"\n\tgt5: "+`candgt5`
            +"\n\tTotal number of table entries: "
            +`len([x for chromo in regobjects for x in chromo])`)
    return summarystring

def writeCompactOutputFile(regobjects,keyword,dirforoutput,analysisname):
    """Write a compact output file containing the results of an analysis.
    See description of compact output file in module docstring."""
    datestring="-".join(time.ctime().rsplit()[1:3])
    print "Writing compact output file for "+analysisname
    if analysisname=="FocusPeriphA":
        outputfile=open(os.path.join(dirforoutput,
                "analyzeLD_C_FocusPeriphA_"+keyword+"_"+datestring+".txt"),'w')
    elif analysisname=="TopRegSNPA":
        outputfile=open(os.path.join(dirforoutput,
                "analyzeLD_C_TopRegSNPA_"+keyword+"_"+datestring+".txt"),'w')
    elif analysisname=="SpecPairA":
        outputfile=open(os.path.join(dirforoutput,
                "analyzeLD_C_SpecPairA_"+keyword+"_"+datestring+".txt"),'w')
    else:
        assert StandardError, ("Error:"
            +" analysisname must be either focusperiph, topregionsnp,"
            +" or specificpairing, not "+analysisname)
    outputfile.write("focusfeatname\tfocusfeatgwas\tfocusfeatposition\t"
                     +"focusfeatTopSNPname\tfocusfeatTopSNPpval\t"
                     +"OtherFeatsCount\tOtherFeatNames\t"
                     +"OtherFeatGWASs\tOtherFeatPositions\t"
                     +"OtherFeatTopSNPname\tOtherFeatTopSNPpval\t"
                     +"Distance\tLD\tGWAScategory\tLDcategory\n")
    for chrm in regobjects:
        for ob in chrm:
            ct=0; nms=[]; gwss=[]; pozs=[]; tsnp=[]; tpval=[]; dist=[]; ld=[]
            keepcands=[]
            for candi in ob.candidates:
                if analysisname in ob.trackanalyses[candi["name"]]:
                    keepcands.append(candi["name"])
                    ct+=1
                    nms.append(candi["name"])
                    gwss.append(candi["name"]+"="+candi["GWAS"])
                    pozs.append(candi["name"]+"="+"chr"
                                    +`candi["chrom"]`+":"+`candi["start"]`
                                    +"-"+`candi["stop"]`)
            for name in keepcands:
                if name in ob.candidateGWASpvals:
                    tsnp.append(name+"="
                                            +ob.candidateGWASpvals[name][0][0])
                    tpval.append(name+"="
                                        +`ob.candidateGWASpvals[name][0][1]`)
                    dist.append(name+"="
                        +`1+abs(ob.candidateGWASpvals[name][0][2]
                              -ob.focusfeat["SNPs"][0][2])`)
                        #add 1 to distance because these are one-based
                        #inclusive coordinates
                if name in ob.candidateLD:
                    ld.append(name+"="+`ob.candidateLD[name][0][2]`)
            #add dummy values
            if nms==[]: nms=["none"]
            if gwss==[]: gwss=["none=none"]
            if pozs==[]: pozs=["none=chr0:0-0"]
            if tsnp==[]: tsnp=["none=rsNA"]
            if tpval==[]: tpval=["none=2.0"]
            if dist==[]: dist=["none=0"]
            if ld==[]: ld=["none=2"]
            #sort alphabetically (applicable when non-dummy values)
            nms.sort();gwss.sort();pozs.sort();tsnp.sort();tpval.sort()
            dist.sort();ld.sort()
            #write out only if it has at least one relevant candidate or no
            #candidates at all and yet the focusfeature is supposed to be
            #included in the relevant analysis(this is to avoid writing out a
            #feature that only has irrelevant candidates that are part of a
            #different analysis)
            if (len(keepcands)>0
                or ((len(keepcands)==0)
                    and analysisname in ob.trackanalyses["self"])):
                outputfile.write(ob.focusfeat["name"]+"\t"
                    +ob.focusfeat["GWAS"]+"\t"
                    +"chr"+`ob.focusfeat["chrom"]`+":"+`ob.focusfeat["start"]`
                        +"-"+`ob.focusfeat["stop"]`+"\t"
                    +ob.focusfeat["SNPs"][0][0]+"\t"
                    +`ob.focusfeat["SNPs"][0][1]`+"\t"
                    +`ct`+"\t")
                for lyst in ([nms,gwss,pozs,
                              tsnp,tpval,dist]):
                    outputfile.write("*".join(lyst)+"\t")
                outputfile.write("*".join(ld)+"\t")
                #calculate GWAScategory
                focusfeatSNPpval=ob.focusfeat["SNPs"][0][1]
                minp=2.0
                keepcands_withGWASinfo=[]
                for candname in keepcands:
                    if candname in ob.candidateGWASpvals.keys():
                        keepcands_withGWASinfo.append(candname)
                        #realistically, every candidate should have
                        #SOMETHING in ob.candidateGWASpvals, but in one of my
                        #tests I didn't want to write up fake GWAS or LD info
                        #so I'm doing this extra step here instead of just
                        #iterating directly over keepcands :)
                for nombre in keepcands_withGWASinfo: #do not directly use
                    #the names in ob.candidateGWASpvals because you want to
                    #calculate GWAScategory only based on the candidates for
                    #that particular analysis, NOT based on all candidates
                    #for any analysis
                    p=ob.candidateGWASpvals[nombre][0][1]
                    if p < minp:
                        minp=p
                if focusfeatSNPpval < minp:
                    GWAScategory=1
                else:
                    GWAScategory=2
                outputfile.write(`GWAScategory`+"\t")
                #calculate LDcategory
                maxLD=-1
                keepcands_withLDinfo=[]
                for candname2 in keepcands:
                    if candname2 in ob.candidateLD.keys():
                        keepcands_withLDinfo.append(candname2)
                for nombre2 in keepcands_withLDinfo:
                    ld=ob.candidateLD[nombre2][0][2]
                    if maxLD < ld < 2:
                        maxLD=ld
                if 0 <= maxLD <= 0.3:
                    LDcategory=1
                elif 0.3 < maxLD < 0.8:
                    LDcategory=2
                elif maxLD >= 0.8:
                    LDcategory=3
                else: #if only LD present is dummy LD of 2.0 or if no LD was
                    #assigned at all and the maxLD is still set to -1
                    LDcategory=9
                outputfile.write(`LDcategory`+"\n")
    outputfile.close()

def writeNonCompactOutputFile(regobjects,keyword,dirforoutput,analysisname):
    """Write a non-compact output file containing the results of an analysis.
    See description of non-compact output file in module docstring."""
    datestring="-".join(time.ctime().rsplit()[1:3])
    print "Writing non-compact output file for "+analysisname
    if analysisname=="FocusPeriphA":
        outputfile=open(os.path.join(dirforoutput,
                        "analyzeLD_NC_FocusPeriphA_"+keyword+"_"+datestring
                        +".txt"),'w')
    elif analysisname=="SpecPairA":
        outputfile=open(os.path.join(dirforoutput,
                        "analyzeLD_NC_SpecPairA_"+keyword+"_"+datestring
                        +".txt"),'w')
    else:
        assert StandardError, ("Error:"
            +" analysisname must be either focusperiph or specificpairing,"
            +"not "+analysisname)
    outputfile.write("focusfeatname\tfocusfeatgwas\tfocusfeatposition\t"
                     +"focusfeatTopSNPname\tfocusfeatTopSNPpval\t"
                     +"otherfeatName\totherfeatGWAS\totherfeatposition\t"
                     +"otherfeatTopSNPname\totherfeatTopSNPpval\tdistance\t"
                     +"LD\n")
    #then use the special property of the candidates to determine which output
    #file they will go in to
    for chrm in regobjects:
        for ob in chrm:
            if analysisname in ob.trackanalyses["self"]:
                onecandwritten=False
                redundantwrite=(ob.focusfeat["name"]+"\t"
                            +ob.focusfeat["GWAS"]+"\t" 
                            +"chr"+`ob.focusfeat["chrom"]`+":"
                                +`ob.focusfeat["start"]`
                                +"-"+`ob.focusfeat["stop"]`+"\t"
                            +ob.focusfeat["SNPs"][0][0]+"\t"
                            +`ob.focusfeat["SNPs"][0][1]`+"\t")
                towrite=[]
                for cand in ob.candidates:
                    if analysisname in ob.trackanalyses[cand["name"]]:
                        x=(redundantwrite+cand["name"]+"\t"
                            +cand["GWAS"]+"\t"+"chr"+`cand["chrom"]`+":"
                            +`cand["start"]`+"-"+`cand["stop"]`+"\t")
                        written2=False
                        for c2 in ob.candidateGWASpvals:
                            if c2==cand["name"]:
                                x=x+(ob.candidateGWASpvals[c2][0][0]
                                    +"\t"+`ob.candidateGWASpvals[c2][0][1]`+"\t"
                                    +`1+abs(ob.focusfeat["SNPs"][0][2]
                                        -ob.candidateGWASpvals[c2][0][2])`+"\t")
                                        #distance has +1 because these are
                                        #one-based inclusive coordinates
                                written2=True
                                break
                        if not written2:
                            raise StandardError,("Error: was writing NonCompact"
                                +" output file and could not find the "
                                +"information for the top SNP in candidate "
                                +cand)
                        written3=False
                        LDwrite="2.0"
                        for c3 in ob.candidateLD:
                            if c3==cand["name"]:
                                LDwrite=`ob.candidateLD[c3][0][2]`
                        x=x+(LDwrite+"\n")
                        onecandwritten=True
                        towrite.append(x)
                towrite.sort()
                for l in towrite:
                    outputfile.write(l)
                if (not onecandwritten
                    and analysisname in ob.trackanalyses["self"]):
                    #if there were no candidates, but this FocusFeature was part
                    #of the specified analysis, then write dummy values
                    outputfile.write(redundantwrite+"none\tnone\tchr0:0-0\t"
                                     +"rsNA\t2.0\t0\t2.0\n")
    outputfile.close()

def getTopSNP(chrom,start,stop,gwaspath):
    """Return the SNP with the lowest GWAS p-val within the specified interval.
    
    PARAMETERS:
    <chrom> is an int indicating which chromosome to search on.
    <start> is an int indicating the starting point of the interval to look
        for SNPs in.
    <stop> is an int indicating the endpoint of the interval.
    <gwaspath> is the path to the GWAS file that will be searched for SNPs.
    
    NOTES: First this function pulls out all of the SNPs within the specified
    interval, including SNPs on the exact endpoints of the interval.
    Then it selects whichever has the lowest p-value, and returns it
    as a tuple inside a list: [("SNPidentifier",p-value,coordinate)]"""
    myfile=open(gwaspath,'r')
    outputlist=[]
    counter=0
    for line in myfile:
        lineaslist=line.rsplit()
        if int(lineaslist[0])==chrom:
            if start <= int(lineaslist[1]) <= stop:
                outputlist.append(("rs"+lineaslist[5],float(lineaslist[6]),
                                   int(lineaslist[1])))
                counter+=1
    assert len(outputlist)==counter, "Error: THIS WAS A STOOPID PHAIL"
    if len(outputlist)==0:
        return [("rsNA",2.0,0)] #dummy default
    else:
        outputlist.sort(key=lambda x: x[1])
        return [outputlist[0]] #return SNP with smallest p-val

#===============================================================================
#----------Helper functions for __eq___() in Region CLASS--------------------
#===============================================================================
def featuredictsequal(dict1,dict2):
    """Return True if dict1 is equal to dict2, else return False.
    Each dict must have the following keys: name, chrom, start, stop, GWAS.
    name and GWAS are STRINGS; chrom, start, and stop are each INTS"""
    if (dict1["name"]!=dict2["name"]
        or dict1["chrom"]!=dict2["chrom"]
        or dict1["start"]!=dict2["start"]
        or dict1["stop"]!=dict2["stop"]
        or dict1["GWAS"]!=dict2["GWAS"]):
        return False
    else:
        return True

def SNPsEqual(SNPlist1,SNPlist2):
    """Return True if SNPlist1 is equal to SNPlist2, else return False.
    Each SNPlist must be a list of tuples with this format:
    (SNP rs number, pval, coordinate)"""
    if (set([snip[0] for snip in SNPlist1])!= set([snip[0]
            for snip in SNPlist2])):
        #if each focusfeat does NOT have the same set of SNP rs numbers
        return False
    elif abs(sum(snip[1] for snip in SNPlist1)- sum(snip[1]
            for snip in SNPlist2))> 0.0001:
        #sum up the p-values for every SNP in SELF and in OTHER;
        #if sum(self)-sum(other) > 0.0001,
        #assume that the p-value sets are different and therefore the objects
        #are different. Can't use direct equality because of weirdness
        #with floats
        return False
    elif (set([snip[2] for snip in SNPlist1])!=set(snip[2]
            for snip in SNPlist2)):
        #if each focusfeat does NOT have the same set of SNP coordinates
        return False
    else:
        return True

def LDlistsequal(LDlist1,LDlist2):
    """Return True if LDlist1 is equal to LDlist2, else return False.
    Each LDlist must be a list of tuples with this format:
    (SNP_A rs number, SNP_B rs number, R2)"""
    if (set([snip[0] for snip in LDlist1])!= set([snip[0]
            for snip in LDlist2])):
        return False
    elif (set([snip[1] for snip in LDlist1])!= set([snip[1]
            for snip in LDlist2])):
        return False
    elif abs(sum(snip[2] for snip in LDlist1)- sum(snip[2]
            for snip in LDlist2))> 0.0001:
        return False
    else:
        return True

#===============================================================================
#----------Region CLASS------------------------------------------------------
#===============================================================================
class Region(object):
    """A Region is centered on a focus feature, and contains information about
    other genomic features nearby (called 'candidates'), including their
    location, GWAS p-values, and LD with the focus feature.
    A Region has five properties: focusfeat, candidates, candidateGWASpvals,
    candidateLD, and trackanalyses.
    
    focusfeat:
        {"name":"EMPTY","chrom":0,"start":0,"stop":0,"SNPs":[],"GWAS":""}
        format of the value stored in "SNPs" key: list of tuples; each SNP
        is a tuple (SNP rs number,GWAS p value, coordinate)
        "GWAS" value is a string indicating the GWAS of the focusfeat
    
    candidates:
        a list of dictionaries, where each dictionary has the
        same format as the dictionary for focusfeat, except there is no "SNPs"
        key. There is one dictionary for every candidate gene in the
        focusfeat's region.
    
    candidateGWASpvals:
        candidateGWASpvals is a dictionary.
        It contains keys that are candidate names and values that are lists
        of SNPs in the tuple format:
        [(rs#, GWAS_p-val, coord),(rs#2, GWAS_p-val2, coord2)]
        for SNPs in that candidate gene,
        with p-values taken from THE GWAS OF THE FOCUSFEAT which is not
        necessarily the same as the GWAS of the candidate.
    
    candidateLD:
        candidateLD is a dictionary. It contains keys that are candidate names
        and values that record LD between that candidate and the focusfeat.
        The values are lists of tuples.
        candidateLD={"candidatename":[(~,~,~),(~,~,~),(~,~,~)]}
        Each tuple is: (SNP_A rs number, SNP_B rs number, R-squared)
        where SNP_A comes from focusfeat and SNP_B comes from the candidate
    
    trackanalyses:
        trackanalyses is a dictionary.
        It contains keys that are candidate names or "self"
        and values that are lists of analysis names, indicating which
        analysis/analyses each candidate gene or the FocusFeat ("self" key)
        is relevant to for this specific Region.
        I cannot denote this in the candidates themselves because
        the actual candidate dictionary objects are shared between
        Regions to save memory, and the "candidates" property is just a list
        of references to the actual objects. This is the same reason that
        there is no "SNPs" property in a candidate dictionary; the same
        candidate might be connected with FocusFeatures from different GWASs,
        in which case the GWAS SNPs for that candidate are different depending
        on which FocusFeature it's paired with."""
    
    _focusfeat=None
    _candidates=None
    _candidateGWASpvals=None
    _candidateLD=None
    _trackanalyses=None
    
    #Note that the way this module is currently written, there is no need for
    #any of the "SNPs" keys to have lists as values, since I only record one
    #SNP for each FocusFeat or Candidate: the top GWAS SNP. However, I am going
    #to keep these as lists, because it means later if I want to expand
    #the number of SNPs I record as associated with each feature, I can do that
    #easily without having to change syntax everywhere.
    
    @property
    def focusfeat(self):
        return self._focusfeat
    @focusfeat.setter
    def focusfeat(self,value):
        if value is None:
            self._focusfeat={}
            #initialize to empty dictionary instead of field default.
            #You cannot have {} as the field default otherwise the same empty
            #dictionary object will be used for EVERY Region object.
            #Likewise for the other properties.
        else:
            self._focusfeat=value
    
    @property
    def candidates(self):
        return self._candidates
    @candidates.setter
    def candidates(self,value):
        if value is None:
            self._candidates=[]
        else:
            assert type(value) is list, ("Error:"
                                    +" cannot set candidates to "+`value`
                                    +". Candidates must be a list.")
            self._candidates=value
    
    @property
    def candidateGWASpvals(self):
        return self._candidateGWASpvals
    @candidateGWASpvals.setter
    def candidateGWASpvals(self,value):
        if value is None:
            self._candidateGWASpvals={}
        else:
            assert type(value) is dict, ("Error:"
                            +" cannot set candidateGWASpvals to "+`value`)
            self._candidateGWASpvals=value
    
    @property
    def candidateLD(self):
        return self._candidateLD
    @candidateLD.setter
    def candidateLD(self,value):
        if value is None:
            self._candidateLD={}
        else:
            assert type(value) is dict, ("Error: cannot set candidateLD to "
                                    +`value`+". CandidateLD must be a dict.")
            self._candidateLD=value
    
    @property
    def trackanalyses(self):
        return self._trackanalyses
    @trackanalyses.setter
    def trackanalyses(self,value):
        if value is None:
            self._trackanalyses={}
        else:
            assert type(value) is dict, ("Error: cannot set trackanalyses to "
                                    +`value`+". Trackanalyses must be a dict.")
            self._trackanalyses=value
    
    def __init__(self,focusfeat,candidates,candidateGWASpvals,candidateLD,
                 trackanalyses):
        self.focusfeat=focusfeat
        self.candidates=candidates
        self.candidateGWASpvals=candidateGWASpvals
        self.candidateLD=candidateLD
        self.trackanalyses=trackanalyses
    
    def __str__(self):
        """Return a string representation of the Region object"""
        def returnchromstring(dictionary):
            """for a dictionary with a chrom, start, and stop key, this
            function returns the values stored in those keys in the format
            chr##:start-stop"""
            return ("chr"+`dictionary["chrom"]`
                    +":"+`dictionary["start"]`+"-"+`dictionary["stop"]`)
        s=""
        s=s+("\nFOCUSFEAT: "+self.focusfeat["name"]+", "
              +returnchromstring(self.focusfeat)+"\n"
              +"\tGWAS: "+self.focusfeat["GWAS"]+"\n"
              +"\tSNPs: "+str(self.focusfeat["SNPs"])+"\nCANDIDATES: ")
        for item in self.candidates:
            s=s+"\n"+str(item)
        s=s+"\nCANDIDATE GWAS P-VALUES: "
        for key in self.candidateGWASpvals:
            s=s+"\n"+key+": "+str(self.candidateGWASpvals[key])
        s=s+"\nCANDIDATE LD: "
        for key in self.candidateLD:
            s=s+"\n"+key+": "+str(self.candidateLD[key])
        s=s+"\nTRACKANALYSES: "
        for key in self.trackanalyses:
            s=s+"\n"+key+": "+str(self.trackanalyses[key])
        return s+"\n--------------------------------------"
    
    def __eq__(self,other):
        #compare focusfeat property
        if not featuredictsequal(self.focusfeat,other.focusfeat):
            return False
        elif not SNPsEqual(self.focusfeat["SNPs"],other.focusfeat["SNPs"]):
            return False
        #compare candidates property
        elif ((set([candidate["name"] for candidate in self.candidates])
             !=set([candidate["name"] for candidate in other.candidates])) or
            (set([candidate["chrom"] for candidate in self.candidates])
              !=set([candidate["chrom"] for candidate in other.candidates])) or
            (set([candidate["start"] for candidate in self.candidates])
              !=set([candidate["start"] for candidate in other.candidates])) or
            (set([candidate["stop"] for candidate in self.candidates])
              !=set([candidate["stop"] for candidate in other.candidates]))):
            return False
        elif (set([gwasname for candidatedict in self.candidates for gwasname
                   in candidatedict["GWAS"]])
        != set([gwasname for candidatedict in other.candidates for gwasname
                in candidatedict["GWAS"]])):
            return False
        #compare candidateGWASpvals property
        elif not SNPsEqual([snip for key in self.candidateGWASpvals for snip
                            in self.candidateGWASpvals[key]],
            [snip for key in other.candidateGWASpvals
             for snip in other.candidateGWASpvals[key]]):
            #Format: [snip for key in dictionary for snip in dictionary[key]]
            #aggregate all the SNPs within each object's candidateGWASpvals
            #property, and compare those aggregate lists to see if they
            #have the same members
            return False
        #compare candidateLD property
        elif not LDlistsequal([LDtuple for key in self.candidateLD
                               for LDtuple in self.candidateLD[key]],
            [LDtuple for key in other.candidateLD
             for LDtuple in other.candidateLD[key]]):
            return False
        #compare trackanalyses property
        elif (set([x for x in self.trackanalyses])
                       !=set([y for y in other.trackanalyses])):
            return False
        elif (set(a for key in self.trackanalyses
                  for a in self.trackanalyses[key])!=
            set(b for kiy in other.trackanalyses
                for b in other.trackanalyses[kiy])):
            return False
        #if all properties were the same, the objects are the same; return True
        else:
            return True

#===============================================================================
#----------EXECUTION------------------------------------------------------------
#===============================================================================
if __name__ == '__main__':
    #initialize variables the config file will fill
    runfocusperiphanalysis=False
    runtopregionsnpanalysis=False
    runspecificpairinganalysis=False
    focusfeaturespath=""
    peripheralfeaturespath=""
    specificpairingpath=""
    GWASes={}
    LDfilepattern=""
    keyword=""
    outputdir=""
    extbas=0
    grepavailable=False
    #read the config file
    configfile=open(sys.argv[1],'r')
    for line in configfile:
        if line[0]=="#" or len(line)<2:
            pass
        else:
            lineaslist=line.rsplit()
            if lineaslist[0]=="RunFocusPeriphAnalysis:":
                if lineaslist[1].upper()=="YES":
                    runfocusperiphanalysis=True
                elif lineaslist[1].upper()=="NO":
                    runfocusperiphanalysis=False
                else:
                    raise StandardError, ("Error:"
                        +" RunFocusPeriphAnalysis must be YES or NO, not "
                        +lineaslist[1])
            elif lineaslist[0]=="RunTopRegionSNPAnalysis:":
                if lineaslist[1].upper()=="YES":
                    runtopregionsnpanalysis=True
                elif lineaslist[1].upper()=="NO":
                    runtopregionsnpanalysis=False
                else:
                    raise StandardError, ("Error:"
                        +" RunTopRegionSNPAnalysis must be YES or NO, not "
                        +lineaslist[1])
            elif lineaslist[0]=="RunSpecificPairingAnalysis:":
                if lineaslist[1].upper()=="YES":
                    runspecificpairinganalysis=True
                elif lineaslist[1].upper()=="NO":
                    runspecificpairinganalysis=False
                else:
                    raise StandardError, ("Error:"
                        +" RunSpecificPairingAnalysis must be YES or NO, not "
                        +lineaslist[1])
            elif lineaslist[0]=="FocusFeatures:":
                focusfeaturespath=lineaslist[1]
            elif lineaslist[0]=="PeripheralFeatures:":
                peripheralfeaturespath=lineaslist[1]
            elif lineaslist[0]=="SpecificPairing:":
                specificpairingpath=lineaslist[1]
            elif lineaslist[0]=="GWAS:":
                if lineaslist[1].replace(":","") in GWASes:
                    raise StandardError, ("Error:"
                            +" Do not include the same GWAS name more than once"
                            +" in your config file.")
                else:
                    GWASes[lineaslist[1].replace(":","")]=lineaslist[2]
            elif lineaslist[0]=="LDfiles:":
                LDfilepattern=lineaslist[1]
            elif lineaslist[0]=="Keyword:":
                keyword=lineaslist[1]
            elif lineaslist[0]=="OutputDir:":
                outputdir=lineaslist[1]
            elif lineaslist[0]=="ExtBases:":
                extbas=int(lineaslist[1])
            else:
                pass
    #Check that the config file had values for all fields
    if (not(runfocusperiphanalysis)
        and not(runtopregionsnpanalysis)
        and not(runspecificpairinganalysis)):
        raise StandardError, "Error: Choose at least one analysis to run"
    if (focusfeaturespath=="" or peripheralfeaturespath==""
        or specificpairingpath=="" or GWASes=={} or LDfilepattern==""
        or keyword=="" or outputdir=="" or extbas==0):
        raise StandardError, "Error: you forgot a needed value in the config"
    #Check that the config file included the necessary information
    if runfocusperiphanalysis:
        assert os.path.exists(focusfeaturespath)
        assert os.path.exists(peripheralfeaturespath)
    if runtopregionsnpanalysis:
        assert os.path.exists(focusfeaturespath)
    if runspecificpairinganalysis:
        assert os.path.exists(specificpairingpath)
    #Run analyses:
    start_time=time.time()
    whattorun={}
    if runfocusperiphanalysis:
        whattorun["FocusPeriphA"]=focusfeaturespath
    if runtopregionsnpanalysis:
        whattorun["TopRegSNPA"]=focusfeaturespath
    if runspecificpairinganalysis:
        whattorun["SpecPairA"]=specificpairingpath
    myregionobjs=initializeRegionobjects(whattorun,GWASes)
    if runfocusperiphanalysis:
        print "Starting the FocusPeriphAnalysis"
        addCandidatestoRegionobjects(extbas,myregionobjs,GWASes,"FocusPeriphA",
                        peripheralfeaturespath)
        print summarizeFocusPeriphAnalysis(myregionobjs)
    if runtopregionsnpanalysis:
        print "Starting the TopRegionSNPAnalysis"
        addCandidatestoRegionobjects(extbas,myregionobjs,GWASes,"TopRegSNPA")
    editRegions_re_LD(myregionobjs,LDfilepattern)
    #Write to output files:
    if runfocusperiphanalysis:
        writeCompactOutputFile(myregionobjs,keyword,outputdir,"FocusPeriphA")
        writeNonCompactOutputFile(myregionobjs,keyword,outputdir,"FocusPeriphA")
    if runtopregionsnpanalysis:
        writeCompactOutputFile(myregionobjs,keyword,outputdir,"TopRegSNPA")
    if runspecificpairinganalysis:
        writeCompactOutputFile(myregionobjs,keyword,outputdir,"SpecPairA")
        writeNonCompactOutputFile(myregionobjs,keyword,outputdir,"SpecPairA")
    print "Pickling the Region objects"
    cPickle.dump(myregionobjs,
        open(os.path.join(outputdir,"analyzeLD_pickledRegions_"+keyword+"_"
             +"-".join(time.ctime().rsplit()[1:3])+".bin"),'wb'),2)
    print ("All done. Runtime "+`(time.time()-start_time)/60.0`+" minutes")

