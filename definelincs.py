#Rachel Ballantyne
#definelincs.py

"""Filter and merge input lincRNA annotation datasets to create a comprehensive,
unified lincRNA dataset.

definelincs.py requires one argument: the full path to the configuration file.

Example configuration file:
#########################
# definelincs.py config #
#########################
broad8000: path/lincRNAs_transcripts.gtf
gencode: path/gencode.v18.annotation.gtf
sigova0: path/SigovaData_S1_hESClincRNAonly_minRPKM0_final_hg19-REPAIRED.bed
hangauerS3: path/Hangauer_DatasetS3_FPKM1_hg19.bed
pseudogenes: path/Human_Pseudogene.txt
refgene: path/refGene.txt
output: path/to/outputdir/
use: broad8000 gencode sigova0 hangauerS3
#END CONFIG FILE

All fields shown above are required in the config file. There must be a space
after every colon. The "use:" field can be followed by any combination of the
four input dataset names shown in the example, separated by spaces. The
pipeline will be run on only the datasets indicated in the "use:" field."""

import time
import cPickle
import copy
import sys
import os.path

import GTFparser_general
import RefGene_parserII
import BEDparser
import FeatureClass_Small

#===============================================================================
#----------FUNCTIONS------------------------------------------------------------
#===============================================================================
def define_lincRNAs(datasets,outputpath):
    """Use the specified <datasets> to create a unified dataset of lincRNAs.
    
    PARAMETERS:
    <datasets> is a dictionary with:
        > some combination of the keys "broad8000","gencode","sigova0", and
          "hangauerS3", OR this key alone: "testing"
        > Plus these keys: "refgene", "pseudogenes"
        The values are paths to the corresponding file (e.g. the value at
        "gencode" is a path to the Gencode GTF file; see module docstring.)
    <outputpath> is the path to the directory where output files will be stored.
    
    NOTES: This function runs the entire lincRNA processing pipeline by
    calling upon other functions in this module. 
    
    OUTPUT: The following files are created:
        > definedlincs_<datasetnames>_<Month>-<Day>.bed: a BED file describing
          the final set of lincRNAs created.
        > definedlincs_<datasetnames>_<Month>-<Day>.bin: a pickled Python list
          describing the final set of lincRNAs created. See function
          pickle_lincs() for details.
        > definelincs_logfile_<Month>-<Day>.txt: a log file describing the
          lincRNA processing pipeline and providing summary statistics from
          each step in the pipeline.
        > quantified_linc_overlap_<datasetnames>_PERCENTS_<Month>-<Day>.txt
          and quantified_linc_overlap_<datasetnames>_DETAILS_<Month>-<Day>.txt
          These files describe the exon overlap between lincRNAs that were
          merged. See function _quantify_linc_exons_overlap() for details.
        > definedlincs_<datasetnames>_FullNamesAndSources_<Month>-<Day>.txt: a txt file
          describing the abbreviated and full names of the lincRNAs as well
          as their sources. See function pickle_lincs() for details."""
    if "nowrite" in datasets: #So that log file from testing can be deleted;
        #see definelincs_test.py
        logfile=open(os.path.join(outputpath,"delthis.txt"),'w')
    else:
        logfile=open(os.path.join(outputpath,
                    "definelincs_logfile_"+"-".join(time.ctime().rsplit()[1:3])
                    +".txt"),'w')
    logfile.write("\n\n"+"------------------------------------------------"
                  +time.ctime()+"\n")
    logfile.close()
    print "RUNNING return_initial_linclist()"
    mylinks=return_initial_linclist(datasets,outputpath)
    print "RUNNING removeproteinoverlap()"
    removeproteinoverlap(mylinks,outputpath,datasets)
    print "RUNNING collapseifexonoverlap()"
    collapseifexonoverlap(mylinks,datasets,outputpath)
    print "RUNNING pickle_lincs()"
    pickle_lincs(mylinks,datasets,outputpath)
    print "RUNNING calculate_genomecoverage()"
    calculate_genomecoverage(mylinks,datasets,outputpath)
    print "RUNNING makeBED_fromlincsbychrom()"
    makeBED_fromchromlist(mylinks,datasets,outputpath)
    print "ALL DONE"

#===========STEP ONE: return_initial_linclist===================================
def return_initial_linclist(datasetdict,outputpath):
    """Return a list of SmallFeature objects describing lincRNAs.
    
    PARAMETERS: 
    <datasetdict> is a dictionary with:
        > some combination of the keys "broad8000","gencode","sigova0", and
          "hangauerS3", OR this key alone: "testing"
        > Plus these keys: "refgene", "pseudogenes"
        The values are paths to the corresponding file (e.g. the value at
        "gencode" is a path to the Gencode GTF file; see module docstring.)
    <outputpath> is the path to the directory where output files will be stored.
    
    NOTES: LincRNAs will be initialized from all of the lincRNA
    datasets present in <datasetdict>. Note that many of these lincRNAs may
    overlap each other and/or mRNAs. All lincRNAs on sex chromosomes are
    removed. LincRNAs less than 200 nt in length are removed as well. 
    
    OUTPUT: The 2-D list of lincRNAs returned is sorted into sub-lists by
    chromosome--lincRNAs on the same chromosome will be together in the same
    sub-list. LincRNAs are represented by SmallFeature objects.
    The following chars are used in the source property of each SmallFeature
    object to indicate where the lincRNA is derived from:
        H=hangauerS3, B=broad8000, G=gencodeV18, S=sigova0"""
    if "nowrite" in datasetdict:
        logfile=open(os.path.join(outputpath,"delthis.txt"),'w')
    else:
        logfile=open(os.path.join(outputpath,
                    "definelincs_logfile_"+"-".join(time.ctime().rsplit()[1:3])
                    +".txt"),'a')
    datasetnames="-".join(x for x in sorted(datasetdict.keys()) if x in
                          ["hangauerS3","broad8000","gencode","sigova0"])
    logfile.write("Datasets used: "+datasetnames)
    totallincs=[]
    if "broad8000" in datasetdict:
        totallincs=add_dataset("B",totallincs,
                    GTFparser_general.makeSmallfeatures_fromGTF,
                    [datasetdict["broad8000"],True,"ALL"])
    if "gencode" in datasetdict:
        totallincs=add_dataset("G",totallincs,
                    GTFparser_general.makeSmallfeatures_fromGTF,
                    [datasetdict["gencode"],True,"lincRNA"])
    if "sigova0" in datasetdict:
        totallincs=add_dataset("S",totallincs,
                    makeSmallfeatures_fromSIGOVA,[datasetdict["sigova0"]])
    if "hangauerS3" in datasetdict:
        totallincs=add_dataset("H",totallincs,
                    BEDparser.makeSmallfeatures_fromBED,
                    [datasetdict["hangauerS3"]])
    GTFparser_general.giveranger(totallincs)
    totallincs_bychrom=GTFparser_general.sortSmallFeatsbychrom(totallincs)
    logfile.write("\nThere are "+
                  `GTFparser_general.countchromifiedlist(totallincs_bychrom)`
                  +" lincRNAs on chr 1-22, X, and Y in total.")
    xcount=0; ycount=0
    for linc in totallincs_bychrom[22]:
        xcount+=1
    for linc in totallincs_bychrom[23]:
        ycount+=1
    totallincs_bychrom[22]=[] #delete any lincs on "chromosome 23" (X)
    totallincs_bychrom[23]=[] #delete any lincs on "chromosome 24" (Y)
    toosmallcount=0 #delete any lincs that are <200 base pairs long
    for chromo in totallincs_bychrom:
        posish=0
        while posish < len(chromo):
            linc=chromo[posish]
            lincsize=(linc.ranger[1]-linc.ranger[0])+1
            if lincsize < 200:
                chromo.remove(linc)
                toosmallcount+=1
            else:
                posish+=1
    logfile.write("\nAll "+`xcount`+" lincs on X and all "+`ycount`+
                  " lincs on Y have been deleted.")
    logfile.write("\nAll "+`toosmallcount`+
                  " lincs that were < 200 base pairs have been deleted.")
    logfile.write("\nThere are "+
                  `GTFparser_general.countchromifiedlist(totallincs_bychrom)`
                  +" lincRNAs on chr 1- 22.")
    logfile.close()
    return totallincs_bychrom

def add_dataset(sourcechar,TotalLincs,openerfunction,openerparams):
        """Add to <totallincs> the lincRNAs described in file opened by
        <openerfunction> called with <openerparams>.
        Each of these added lincRNAs will have <sourcechar> added to their
        source property."""
        mydata=openerfunction(*openerparams)
        for linky in mydata:
                linky.source=sourcechar
        TotalLincs=TotalLincs+mydata
        return TotalLincs

#=====================STEP TWO: remove protein overlap==========================
def removeproteinoverlap(lincsbychrom,outputpath,datasetdict):
    """Remove lincRNAs from <lincsbychrom> that are near protein-coding genes.
    
    PARAMETERS:
    <lincsbychrom> is a two-dimensional Python list of SmallFeature objects
        describing lincRNAs. See return_initial_linclist() OUTPUT description.
    <outputpath> is the path to the directory where output files will be stored.
    <datasetdict> is described in return_initial_linclist() parameters. Note
        that this dictionary contains the  keys "refgene" and "pseudogenes."
        The corresponding values indicate where the RefGene.txt and
        Human_Pseudogene.txt files are stored. NM transcript locations are
        obtained from the file RefGene.txt. Pseudogenes are obtained from the
        file Human_Pseudogene.txt
    
    NOTES:
        STRANDED OR NONSTRANDED TRANSCRIPTS are removed if:
            They overlap a pseudogene on either strand.
        STRANDED TRANSCRIPTS are removed if:
            They are within 1 kb of RefSeq NM genes on the same strand
            They overlap an NM gene that is on the opposite strand by at least
            one base
        TRANSCRIPTS WITHOUT STRAND INFO are removed if:
            They are within 1 kb of RefSeq NM genes on either strand
    
    OUTPUT: This function adds to the pipeline log file and modifies the
    <lincsbychrom> list in place. There is no explicit return value."""
    if "nowrite" in datasetdict:
        logfile=open(os.path.join(outputpath,"delthis.txt"),'w')
    else:
        logfile=open(os.path.join(outputpath,"definelincs_logfile_"+
                              "-".join(time.ctime().rsplit()[1:3])+".txt"),'a')
    kntr=0
    for chromos in lincsbychrom:
        for lincRNA in chromos:
            lincRNA.misc=True #initialize all misc properties to True, which
                              #here means "do not delete"
            kntr+=1 #count all lincRNAs too
    logfile.write("\nFUNCTION removeproteinoverlap():")
    logfile.write("\n\tlincRNAs before running this function: "+`kntr`)
    #FIRST: mark for removal any lincs that overlap a pseudogene
    pseudogenes=_return_pseudogenes(datasetdict["pseudogenes"])
    GTFparser_general.giveranger(pseudogenes)
    pseudogenesbychrom=GTFparser_general.sortSmallFeatsbychrom(pseudogenes)
    removed_bcz_pseudoverlap=0
    for chromosome in range(0,22): #ignores sex chroms; those are empty anyway
        for lincrna in lincsbychrom[chromosome]:
            for psgene in pseudogenesbychrom[chromosome]:
                if _is_overlap(lincrna.ranger[0],lincrna.ranger[1],
                               psgene.ranger[0],psgene.ranger[1]):
                #if the lincrna and psgene overlap by one or more bases, then:
                    if lincrna.misc:
                        lincrna.misc=False
                        removed_bcz_pseudoverlap=removed_bcz_pseudoverlap+1
    #SECOND: mark for removal any lincs that are near/overlap an NM gene
    mrnaslist=RefGene_parserII.returnRefGenelist(datasetdict["refgene"],["NM"])
    GTFparser_general.giveranger(mrnaslist)
    mrnasbychrom=GTFparser_general.sortSmallFeatsbychrom(mrnaslist)
    removed_bcz_mrnaoverlap=0
    for chromosome in range(0,22):
        print "chromosome "+`chromosome+1`+" is being worked on."
        for lincrna in lincsbychrom[chromosome]:
            for mrna in mrnasbychrom[chromosome]:
                if (("+" in mrna.strand and "+" in lincrna.strand) or
                    ("-" in mrna.strand and "-" in lincrna.strand) or
                    ("+" not in lincrna.strand and "-" not in lincrna.strand)):
                    #if mrna/lincrna on same strand, or linc has no strand info
                    if _is_overlap(lincrna.ranger[0],lincrna.ranger[1],
                                   mrna.ranger[0],mrna.ranger[1],500):
                    #Note: make it 500 instead of 1000 because each feature
                    #is expanded by <distance> bases; and 500+500=1000
                    #if the lincRNA and mRNA overlap by at least one base when
                    #the lincRNA is extended by 1000 bp in either direction,then
                        if lincrna.misc:
                            lincrna.misc=False #if it hasn't been set to False
                                #by pseudogene overlap,set lincrna.misc to False
                            removed_bcz_mrnaoverlap+=1
                            break #this is VERY important;
                            #you don't want to check one linc over and over
                            #if you already know it overlaps one mrna
                elif (("+" in mrna.strand and "-" in lincrna.strand)
                    or ("-" in mrna.strand and "+" in lincrna.strand)):
                    #mrna/lincrna are on opposite strands
                    if _is_overlap(lincrna.ranger[0],lincrna.ranger[1],
                                   mrna.ranger[0],mrna.ranger[1]):
                    #if the lincRNA and mRNA overlap by at least one base
                        if lincrna.misc:
                            lincrna.misc=False
                            removed_bcz_mrnaoverlap+=1
                            break
                else: #you should never go in this else
                    assert False, "ERROR in removeproteinoverlap"
    #THIRD: actually remove lincRNAs that need to be removed
    delkntr=0
    for chromosome in range(0,22):
        poz=0
        while poz < len(lincsbychrom[chromosome]):
            lincrna=lincsbychrom[chromosome][poz]
            if not(lincrna.misc): #if lincrna.misc has been changed to False
                lincsbychrom[chromosome].remove(lincrna)
                delkntr+=1
            else:
                poz+=1 #only increment the pos if you didn't just delete a linc
    logfile.write("\n\tlincRNAs removed due to overlap with pseudogene: "
                  +`removed_bcz_pseudoverlap`)
    logfile.write("\n\tlincRNAs removed due to overlapping"
                  +"or being within 1000 kb of RefSeq NM gene: "
                  +`removed_bcz_mrnaoverlap`)
    logfile.write("\n\tTotal lincRNAs removed: "+`delkntr`)
    lincsleft=GTFparser_general.countchromifiedlist(lincsbychrom)
    logfile.write("\n\tTotal lincRNAs remaining: "+`lincsleft`)
    logfile.close()
    for chromosome in lincsbychrom:
        for lincRNA in chromosome:
            assert lincRNA.misc, ("Error: You forgot to delete a lincRNA that"
                                  +"has a False misc property!")
    assert delkntr==removed_bcz_pseudoverlap+removed_bcz_mrnaoverlap, ("Error:"
                            +"the removal counts do not add up properly!")
    assert lincsleft==kntr-delkntr, "Error: Counts are not summing properly"
    #Note: you don't need to return anything. This function modifies the list.

def _return_pseudogenes(pseudogenepath):
    """Return a list of SmallFeature objects read out of Human_Pseudogene.txt.
    These SmallFeatures contain chromosome, start, and stop information."""
    pseudfile=open(pseudogenepath,'r')
    pseudlist=[]
    invalidlinecount=0
    validlinecount=0
    totallinecount=0
    for line in pseudfile:
        totallinecount+=1
        lineaslist=line.rsplit("\t")
        try:
            chrm=int(lineaslist[1]) #ignore X, Y, M, and any other non-autosome
            strt=int(lineaslist[2])
            stp=int(lineaslist[3])
            newpseud=FeatureClass_Small.SmallFeature("XZP",".",chrm,[strt],
                                                     [stp],[1],[strt,stp],[])
            pseudlist.append(newpseud)
            validlinecount+=1
        except:
            invalidlinecount+=1
    assert totallinecount==invalidlinecount+validlinecount, ("Error in parsing"
                +" Human_Pseudogene.txt: not all lines were parsed")
    return pseudlist

#=========STEP THREE: collapse overlapping transcripts using exon info==========
def collapseifexonoverlap(lincsbychrom,datasetdict,outputpath):
    """Collapse lincRNAs together based on exon overlap.
    
    PARAMETERS:
    <lincsbychrom> is a two-dimensional Python list of SmallFeature objects
        describing lincRNAs. See return_initial_linclist() OUTPUT description.
        This list should have been processed by removeproteinoverlap() before
        collapseifexonoverlap() is applied.
    <datasetdict> is described in return_initial_linclist() parameters.
    <outputpath> is the path to the directory where output files will be stored
    
    NOTES: linc1 is merged with linc2 if 50% of an exon of either linc overlaps
    an exon of the other linc.
        collapseifexonoverlap() uses several helper functions.
        First, collapseifexonoverlap() calls upon _is_overlap() using
            general lincRNA start and stop positions to determine if the general
            region occupied by linc1 overlaps the general region occupied by
            linc2. If the general regions do overlap, then linc1 and linc2 will
            be compared using _linc_exons_overlap() which determines whether
            linc1 linc2 have any 50% overlapping exons (in other words, it
            determines whether linc1 and linc2 should be merged.)
            _linc_exons_overlap() will merge linc1 and linc2 if they should be
            merged.
        The reason for having the first check be for general overlap
            rather than using_linc_exons_overlap() immediately is that most
            lincs aren't even in the same vicinity as each other, and
            _linc_exons_overlap() takes a MUCH greater amount of time than
            looking for general overlap. It is only worth it to call
            _linc_exons_overlap() if you know that the lincs are generally
            overlapping.
        Strandedness: If a linc has strand information, that information will be
            used to make sure it is never combined with a linc that is
            specifically on the opposite strand, and it will additionally never
            be combined with a nonstranded linc that has been combined with a
            linc on the opposite strand. See function
            _is_legal_strand_combo()
    
    OUTPUT: This function adds to the pipeline log file and modifies the
        <lincsbychrom> list in place. One of the helper functions creates
        the quantified_linc_overlap output txt files. There is no explicit
        return value."""
    if "nowrite" in datasetdict:
        logfile=open(os.path.join(outputpath,"delthis.txt"),'w')
    else:
        logfile=open(os.path.join(outputpath,"definelincs_logfile_"+
                              "-".join(time.ctime().rsplit()[1:3])+".txt"),'a')
    logfile.write("\nFUNCTION collapseifexonoverlap():")
    logfile.write("\n\tlincRNAs before running this function: "+
                  `GTFparser_general.countchromifiedlist(lincsbychrom)`)
    genoverlap_count=0
    genoverlap_yes_exon_ovlp_cnt=0
    genoverlap_no_exon_ovlp_cnt=0
    no_overlap_count=0
    #Open more logfiles, which will be written to
    #by _quantify_linc_exons_overlap()
    datasetnames="-".join(x for x in sorted(datasetdict.keys()) if x in
                          ["hangauerS3","broad8000","gencode","sigova0"])
    quantfile_details=open(os.path.join(outputpath,"quantified_linc_overlap_"
                +datasetnames
                +"_DETAILS_"+"-".join(time.ctime().rsplit()[1:3])+".txt"),'w')
    quantfile_percents=open(os.path.join(outputpath,"quantified_linc_overlap_"
                +datasetnames
                +"_PERCENTS_"+"-".join(time.ctime().rsplit()[1:3])+".txt"),'w')
    quantfile_percents.write("GREATER\tLESSER\tMERGED\n")
    logfylelist=[quantfile_details,quantfile_percents]
    for chrom in lincsbychrom:
        for linck in chrom:
            _clean_up_exon_locations(linck)
            linck.misc=linck.featurename
    for chromosome in range(0,22): #ignores sex chrs, which are empty anyway
        position1=0
        position2=1
        while (position1 < len(lincsbychrom[chromosome])-1 and
               position2 < len(lincsbychrom[chromosome])):
            linc1=lincsbychrom[chromosome][position1]
            linc2=lincsbychrom[chromosome][position2]
            if (_is_overlap(linc1.ranger[0],linc1.ranger[1],linc2.ranger[0],
                            linc2.ranger[1]) #general overlap test
                and _is_legal_strand_combo(linc1.strand,linc2.strand)):
                genoverlap_count+=1
                if _linc_exons_overlap(linc1,linc2,lincsbychrom,chromosome,
                                       logfylelist):
                    #if _linc_exons_overlap() was True, linc1 and linc2 were
                    #just combined, linc1 was deleted, and 
                    #_quantify_linc_exons_overlap() was called
                    #within _linc_exons_overlap()
                    genoverlap_yes_exon_ovlp_cnt+=1
                    position2=position1+1
                else: #linc1 and linc2 were not combined;
                    #increment positions accordingly
                    #ChangePositionAfterNoCombo:
                    genoverlap_no_exon_ovlp_cnt+=1 #I expect this to be small
                    length=len(lincsbychrom[chromosome])
                    if position2 < length-1 and position1 < length-2:
                        position2=position2+1
                    elif position2==length-1 and position1<length-2:
                        position1=position1+1
                        position2=position1+1
                    else:
                        assert (position2==length-1 and
                                position1==length-2),("ERROR: "+
                                    "pos2 is actually "+`position2`+
                                    " and pos1 is "+`position1`)
                        position1=length #stop the while loop
            else: #they do not overlap at all, or are on different strands,
                #and thus should not be combined
                no_overlap_count+=1
                #ChangePositionAfterNoCombo:
                length=len(lincsbychrom[chromosome])
                if position2 < length-1 and position1 < length-2:
                    position2=position2+1
                elif position2==length-1 and position1<length-2:
                    position1=position1+1
                    position2=position1+1
                else:
                    assert position2==length-1 and position1==length-2, "ERROR"
                    position1=length #stop the while loop
    logfile.write("\n\t\toccurrences of general overlap: "+`genoverlap_count`
                  +"\n\t\toccurrences of general overlap with exon overlap: "
                  +`genoverlap_yes_exon_ovlp_cnt`
                  +"\n\t\toccurrences of general overlap with NO exon overlap: "
                  +`genoverlap_no_exon_ovlp_cnt`
                  +"\n\t\toccurrences of no overlap at all: "+`no_overlap_count`
                  +"\n\tlincRNAs after running this function: "
                  +`GTFparser_general.countchromifiedlist(lincsbychrom)`)
    logfile.close();quantfile_details.close();quantfile_percents.close()

def pickle_lincs(lincsbychrom,datasetdict,outputpath):
    """Save info in <lincsbychrom> to a bin file and a txt file.
    
    PARAMETERS:
    <lincsbychrom> is a two-dimensional Python list of SmallFeature objects
        describing lincRNAs. See return_initial_linclist() OUTPUT description.
        Before applying pickle_lincs() to <lincsbychrom>, the <lincsbychrom>
        should be processed by removeproteinoverlap() and
        collapseifexonoverlap().
    <datasetdict> is described in return_initial_linclist() parameters.
    <outputpath> is the path to the directory where output files will be stored.
    
    OUTPUT:
        > definedlincs_<datasetnames>_<Month>-<Day>.bin: this is <lincsbychrom>
          pickled and saved. The lincRNAs are stored as SmallFeature objects.
          Additional information present in this file that will not be
          included in the BED file output by this module is the original source
          dataset(s) of the lincRNA, and the full lincRNA name which is a
          concatenation of all of the contributing lincRNA names. The full
          lincRNA name is stored in the misc property of each SmallFeature
          object. In the BED file, the shortest name of all the names of
          contributing lincRNAs is chosen to refer to the final
          composite lincRNA.
        > definedlincs_<datasetnames>_FullNamesAndSources_<Month>-<Day>.txt:
          the first column is the short lincRNA name used to refer to the
          lincRNA in the BED file.
          The second column is the corresponding long name, which is a
          concatenation of the full names of all the original lincRNAs that were
          merged to create the final lincRNA. For the many lincRNAs that come
          from just one source dataset, the first and second columns will be
          identical. The third column is the "source" property of the
          SmallFeature object, which contains characters indicating the
          data files that the lincRNA was drawn from. See OUTPUT description of
          return_initial_linclist() for letter code info.
          This txt file is created so that all information output by this module
          is available in BED or TXT format and nobody has to deal with the
          pickled BIN file if they don't want to."""
    datasetnames="-".join(x for x in sorted(datasetdict.keys()) if x in
                          ["hangauerS3","broad8000","gencode","sigova0"])
    fyl=open(os.path.join(outputpath,"definedlincs_"+datasetnames+"_"
                +"-".join(time.ctime().rsplit()[1:3])+".bin"),'wb') #Must be wb
    print "Beginning to pickle the final list"
    cPickle.dump(lincsbychrom,fyl,2)
    print "Done pickling the final list"
    fyl.close()
    fyl2=open(os.path.join(outputpath,"definedlincs_"+datasetnames
                           +"_FullNamesAndSources_"
                           +"-".join(time.ctime().rsplit()[1:3])
                           +".txt"),'w')
    for chrom in lincsbychrom:
        for linc in chrom:
            fyl2.write(linc.featurename+"\t"+linc.misc+"\t"+linc.source+"\n")
    fyl2.close()

def calculate_genomecoverage(lincsbychrom,datasetdict,outputpath):
    """Return the total number of bases covered by lincRNAs in <lincsbychrom>.
    
    PARAMETERS:
    <lincsbychrom> is a two-dimensional Python list of SmallFeature objects
        describing lincRNAs. See return_initial_linclist() OUTPUT description.
        Before applying calculate_genomecoverage() to <lincsbychrom>, the
        <lincsbychrom> should be processed by removeproteinoverlap() and
        collapseifexonoverlap().
    <datasetdict> is described in return_initial_linclist() parameters.
    <outputpath> is the path to the directory where output files will be stored.
    
    NOTES: calculate_genomecoverage() does not consider exon
    information--it merely uses the ranger (the overall start and overall
    end position) of each lincRNA to see how many bases of the genome are
    covered. calculate_genomecoverage() makes a deep copy of <lincsbychrom>
    so that it can collapse everything according to general overlap without
    affecting the actual <lincsbychrom> list. The reason to collapse by general
    overlap is so in the case of lincRNAs that are separate (due to no exon
    overlap) but actually lie in a similar location (general overlap) the same
    region of genome will not be counted twice towards the total coverage.
    A lot of this function was copied right out of the body of
    collapseifexonoverlap(). Since the coordinate system I use is closed
    (the endpoints are included in the interval) I need to add one when I
    calculate length from endpoint subtraction.
    
    OUTPUT: This function writes the genome coverage result to the pipeline
    log file. There is no explicit return value."""
    print "Calculating genome coverage"
    modylist=copy.deepcopy(lincsbychrom)
    #First collapse everything by general overlap
    genoverlapcnt=0
    for chromosome in range(0,22): #ignores sex chrs, which are empty anyway
        position1=0
        position2=1
        while (position1 < len(modylist[chromosome])-1
               and position2 < len(modylist[chromosome])):
            linc1=modylist[chromosome][position1]
            linc2=modylist[chromosome][position2]
            if _is_overlap(linc1.ranger[0],linc1.ranger[1],linc2.ranger[0],
                           linc2.ranger[1]):
                genoverlapcnt+=1
                _mergelincs(linc1,linc2)
                modylist[chromosome].remove(linc1)
                position2=position1+1
            else: #linc1 and linc2 should not be combined;
                #increment positions accordingly
                #ChangePositionAfterNoCombo:
                length=len(modylist[chromosome])
                if position2 < length-1 and position1 < length-2:
                    position2=position2+1
                elif position2==length-1 and position1<length-2:
                    position1=position1+1
                    position2=position1+1
                else:
                    assert (position2==length-1
                            and position1==length-2), ("ERROR:"
                                +" pos2 is actually "+`position2`
                                +" and pos1 is "+`position1`)
                    position1=length #stop the while loop
    #Then calculate the coverage of the genome
    coverage=0
    for chromo in modylist:
        for linc in chromo:
            coverage=coverage+(linc.ranger[1]-linc.ranger[0])+1
            #add one since my coordinate system has the endpoints included
    if "nowrite" in datasetdict:
        logfile=open(os.path.join(outputpath,"delthis.txt"),'w')
    else:
        logfile=open(os.path.join(outputpath,"definelincs_logfile_"
                              +"-".join(time.ctime().rsplit()[1:3])+".txt"),'a')
    logfile.write("\nFUNCTION returngenomecoverage(): \n\tThe total number of"
                  +" bases in the genome that are covered by a lincRNA is "
                  +`coverage`)
    print ("The total number of bases in the genome that are covered by a"
            +" lincRNA is "+`coverage`)
    logfile.close()
    return coverage

def makeBED_fromchromlist(lincsbychrom,datasetdict,outputpath):
    """Save a BED file from lincRNAs in <lincsbychrom> at <outputpath>.
    
    PARAMETERS:
    <lincsbychrom> is a two-dimensional Python list of SmallFeature objects
        describing lincRNAs. See return_initial_linclist() OUTPUT description.
        Before applying makeBED_fromchromlist() to <lincsbychrom>, the
        <lincsbychrom> should be processed by removeproteinoverlap() and
        collapseifexonoverlap().
    <datasetdict> is described in return_initial_linclist() parameters.
    <outputpath> is the path to the directory where output files will be stored.
    
    NOTES:
    The SmallFeature objects have their positions in the one-based
    closed interval coordinate system. The positions will be manipulated
    appropriately in this function to put them into the half-open zero-based
    BED file coordinate system.
    
    OUTPUT:
    A BED file describing lincRNAs in <lincsbychrom>. There is no header in the
    BED file created by this function.
    See <http://genome.ucsc.edu/FAQ/FAQformat.html> for information about the
    BED file format."""
    datasetnames="-".join(x for x in sorted(datasetdict.keys()) if x in
                          ["hangauerS3","broad8000","gencode","sigova0"])
    BEDfile=open(os.path.join(outputpath,"definedlincs_"+datasetnames+"_"
                              +"-".join(time.ctime().rsplit()[1:3])+".bed"),'w')
    for chrom in lincsbychrom:
        for feat in chrom:
            BEDfile.write("chr"+`feat.chromosome`+"\t") #chrom
            BEDfile.write(`feat.ranger[0]-1`+"\t")#chromStart,0 based,start incl
            BEDfile.write(`feat.ranger[1]`+"\t")  #chromEnd, 0 based, end excl
            BEDfile.write(feat.featurename+"\t")        #name
            BEDfile.write("500\t")                      #score
            if feat.strand=="%" or feat.strand=="*":
                BEDfile.write(".\t")  #strand. Loss of partial strand info
                #since BED format doesn't use my symbols % and *
            else:
                BEDfile.write(feat.strand+"\t")
            BEDfile.write(`feat.ranger[0]`+"\t")        #thickStart=chromStart
            BEDfile.write(`feat.ranger[1]`+"\t")        #thickEnd=chromEnd
            BEDfile.write("255,0,0\t")              #itemRgb is set to 255,0,0
            BEDfile.write(`len(feat.exons)`+"\t")       #blockCount
            #calculate the start positions
            starts_zerobased=""
            for startpos in feat.start:
                starts_zerobased=(starts_zerobased
                                  +`startpos-1-(feat.ranger[0]-1)`+",")
            #calculate the sizes
            sizes=""
            for index in range(len(feat.exons)):
                exonsize=feat.stop[index]-(feat.start[index]-1)
                sizes=sizes+`exonsize`+","
            #write the sizes and starts to the file
            BEDfile.write(sizes+"\t")                   #blockSizes
            BEDfile.write(starts_zerobased+"\t\n")      #blockStarts
    BEDfile.close()

#===============================================================================
#----------HELPER FUNCTIONS-----------------------------------------------------
#===============================================================================
def _is_legal_strand_combo(symbol1,symbol2):
    """Return True if strand symbol1 can be combined with strand symbol2,
    else return False.
    Possible Strand Symbols:
        + : on the positive strand
        - : on the positive strand
        . : no strand information provided
        * : result of a + and . concatenation
        % : result of a - and . concatenation
    Combination Rules:
        Legal combinations: every symbol in the group can be combined with
        any other symbol in that group
        Group 1: + and . and *
        Group 2: - and . and %"""
    legallist1=["+",".","*"]
    legallist2=["-",".","%"]
    if symbol1 in legallist1 and symbol2 in legallist1:
        return True
    elif symbol1 in legallist2 and symbol2 in legallist2:
        return True
    else:
        return False

def _is_overlap(startA,stopA,startB,stopB,distance=0):
    """Return the number of bases by which unit A and unit B overlap
    (a non-zero integer, True). If unit A and unitB do not overlap, return
    0 (False).
    
    PARAMETERS:
    All parameters are ints, with startA < stopA and startB < stopB. 
    The provided start and stop positions could represent lincRNA
    start and stop positions (for testing general lincRNA overlap) or they
    could represent the start/stop for exons (for testing exon overlap.)
    If <distance>, an integer, is provided, it will decrease each start position
    by <distance> and increase each stop position by <distance>, thus extending
    each unit by <distance> in either direction. The default for <distance>
    is zero."""
    #The reason for the +1 is because of your coordinate system of choice
    startA=startA-distance
    stopA=stopA+distance
    startB=startB-distance
    stopB=stopB+distance
    if startA <= startB <= stopA:
        if stopB >= stopA:
            return (stopA-startB)+1
        if stopB <= stopA:
            return (stopB-startB)+1
    elif startB <= startA <= stopB:
        if stopB >= stopA:
            return (stopA-startA)+1
        if stopB <= stopA:
            return (stopB-startA)+1
    else:
        return 0 #They don't overlap

def _linc_exons_overlap(linc1,linc2,lincsbychrom,chromosome,logfilelist):
    """Determine whether two lincRNAs should be merged, and if yes, merge them.
    
    PARAMETERS:
    <linc1> and <linc2> are SmallFeature objects.
    <lincsbychrom> is a two-dimensional Python list of SmallFeature objects
        describing lincRNAs. See return_initial_linclist() OUTPUT description.
    <chromosome> is an int for which chromosome the lincRNA is on. Note that
        this is a Python list index, so what is chromosome 5 on the genome
        will be at index 4 (since indices start at 0 in Python.)
    <logfilelist> is a two-member list of open file objects to which
        _quantify_linc_exons_overlap() will write.
    
    NOTES: This function uses _is_overlap() to determine whether linc1 and
    linc2 have any exons that overlap by >=50%. If they do, then it uses
    _mergelincs() which is called by _quantify_linc_exons_overlap() to merge
    linc1 into linc2, and returns True. If linc1 and linc2 do NOT have any exons
    that overlap by >=50%, then they will not be modified at all, and this
    function will return False."""
    for linc1dex in range(len(linc1.exons)):
        start1=linc1.start[linc1dex]
        stop1=linc1.stop[linc1dex]
        for linc2dex in range(len(linc2.exons)):
            start2=linc2.start[linc2dex]
            stop2=linc2.stop[linc2dex]
            amount_of_overlap=_is_overlap(start1,stop1,start2,stop2)
            if (amount_of_overlap>((stop1-start1+1)/2)
                or amount_of_overlap>((stop2-start2+1)/2)): #looking for whether
                #either exon is overlapped 50% or more
                _quantify_linc_exons_overlap(linc1,linc2,lincsbychrom,
                                             logfilelist)
                lincsbychrom[chromosome].remove(linc1) #now you can remove linc1
                #(it's already part of linc2)
                return True
    return False

def _quantify_linc_exons_overlap(linc1,linc2,lincsbychrom,logfilelist):
    """Write info about merging to output files, and clean up the merged linc.
    
    PARAMETERS: see parameters of _linc_exons_overlap()
    
    NOTES:
    This function is called before _linc_exons_overlap() returns True--in
    other words, when either linc1 or linc2 has an exon that is 50% overlapped
    by an exon of the other linc.
    
    OUTPUT:
        > quantified_linc_overlap_<datasets>_PERCENTS.txt: A three-column
          file with a GREATER, LESSER, and MERGED column. Each row
          contains info about a merge event between two lincRNAs. The percent
          by which linc1 is overlapped by linc2 is calculated, as is the percent
          by which linc2 is overlapped by linc1. Whichever of these percents is
          greater gets written to the "GREATER" column, and the other percent
          is written to the "LESSER" column. The "MERGED" column contains
          the amount of exonic overlap between linc1 and linc2 divided by
          the total exonic length of the merged linc.
        > quantified_linc_overlap_<datasets>_DETAILS.txt: A four-column
          file to which two lines are added for each merge event. Column 1
          is lincRNA name, column 2 is total length of the lincRNA's exons,
          column 3 is the total bases of the lincRNA's exons that are overlapped
          by the other lincRNA's exons, and column 4 is the percent of the
          lincRNA's exons that are overlapped by the other lincRNA's exons.
          A separate line is written for linc1 and linc2 of the merge event."""
    detailslog=logfilelist[0]
    percentslog=logfilelist[1]
    linc1_total_exon_size=0
    linc2_total_exon_size=0
    exonic_overlap=0
    for linc1dex in range(len(linc1.exons)):
        start1=linc1.start[linc1dex]
        stop1=linc1.stop[linc1dex]
        linc1_total_exon_size+=(stop1-start1)+1
        for linc2dex in range(len(linc2.exons)):
            start2=linc2.start[linc2dex]
            stop2=linc2.stop[linc2dex]
            exonic_overlap+=_is_overlap(start1,stop1,start2,stop2)
    for linc2dex in range(len(linc2.exons)):
        start2=linc2.start[linc2dex]
        stop2=linc2.stop[linc2dex]
        linc2_total_exon_size+=(stop2-start2)+1 #cannot have this inside the
        #previous for loop or it gets counted over too many times
    percentlinc1=float(exonic_overlap)/float(linc1_total_exon_size)
    percentlinc2=float(exonic_overlap)/float(linc2_total_exon_size)
    if percentlinc1>1.0:
        percentlinc1=1.0
    if percentlinc2>1.0:
        percentlinc2=1.0 #I don't want to be writing percent overlaps that are
        #greater than 1 since that looks weird
    detailslog.write(linc1.featurename+"\t"+`linc1_total_exon_size`+"\t"
                     +`exonic_overlap`+"\t"+`percentlinc1`
                     +"\n"+linc2.featurename+"\t"+`linc2_total_exon_size`
                     +"\t"+`exonic_overlap`+"\t"+`percentlinc2`+"\n")
    percentslog.write(`max(percentlinc1,percentlinc2)`+"\t"
                      +`min(percentlinc1,percentlinc2)`+"\t")
    _mergelincs(linc1,linc2)
    _clean_up_exon_locations(linc2)
    linc2_merged_exon_size=0
    for newlinc2dex in range(len(linc2.exons)):
        newstart2=linc2.start[newlinc2dex]
        newstop2=linc2.stop[newlinc2dex]
        linc2_merged_exon_size+=(newstop2-newstart2)+1
    assert exonic_overlap<=linc2_merged_exon_size,("Error"
        +" with _clean_up_exon_locations: merged exon size is "
        +`linc2_merged_exon_size`+" which is greater than exonic overlap of "
        +`exonic_overlap`)
    percentslog.write(`float(exonic_overlap)/float(linc2_merged_exon_size)`
                      +"\n")

def _mergelincs(linc1,linc2):
    """Merge linc1 into linc2 by appending the contents of key properties of
    linc1 into the corresponding properties of linc2.
    
    PARAMETERS: <linc1> and <linc2> are separate SmallFeature objects.
    
    NOTE: the function _clean_up_exon_locations() is used on the output of
    _mergelincs() to clean up the exon structure of the lincRNA resulting
    from the merge.
    
    OUTPUT: This function modifies linc2 in place. There is no explicit
    return value."""
    assert linc1.chromosome==linc2.chromosome,("Error:"
                        +" linc1 and linc2 are on different chromosomes.")
    strand1=linc1.strand; strand2=linc2.strand
    if strand1!=strand2:
        legallistplus=["+",".","*"]
        legallistminus=["-",".","%"]
        if strand1 in legallistplus and strand2 in legallistplus:
            linc2.strand="*"
        elif strand1 in legallistminus and strand2 in legallistminus:
            linc2.strand="%"
        #Note that if linc1.strand==linc2.strand, then you'll leave
        #linc2.strand alone, thus preserving the strand info
    if len(linc2.featurename)>len(linc1.featurename): #only change the
        #featurename if you can shorten it by making it the other featurename
        linc2.featurename=linc1.featurename
    linc2.start=linc2.start+linc1.start
    linc2.stop=linc2.stop+linc1.stop
    linc2.exons=linc2.exons+linc1.exons
    if linc1.ranger[0] < linc2.ranger[0]:
        linc2.ranger[0]=linc1.ranger[0]
    if linc1.ranger[1] > linc2.ranger[1]:
        linc2.ranger[1]=linc1.ranger[1]
    linc2.misc=linc2.misc+","+linc1.misc
    linc2.source=linc2.source+linc1.source

def _clean_up_exon_locations(linc):
    """Clean the exon locations for a lincRNA that has resulted from a merge.
    
    PARAMETER: <linc> is a SmallFeature object that has resulted from a merge
    event and needs to be 'cleaned up.'
    
    NOTES:
    _clean_up_exon_locations() is run on a linc that has just had another linc
    appended to it. The function examines all the exons in the newly-created
    lincRNA, some of which are guaranteed to overlap (because you just merged
    two lincs BECAUSE of exon overlap). This function merges together exons that
    overlap by one or more bases, and modifies the start and stop properties of
    the lincRNA accordingly so that the minimal number of exons are described,
    and none of the final set of exons described overlap each other.
    This function has EXTREMELY SIMILAR logic to the function
    collapseifexonoverlap() except that it is collapsing exons of a single linc
    together, instead of collapsing different lincs together.
    
    OUTPUT: This function modifies <linc> in place. There is no explicit
    return value."""
    #Make the list of all exons, as Exon objects
    listofExonObjects=[]
    for position in range(len(linc.exons)):
        exonbeginning=linc.start[position]
        exonend=linc.stop[position]
        listofExonObjects.append(Exon(exonbeginning,exonend))
    #Collapse the exon list
    pos1=0
    pos2=1
    while pos1<len(listofExonObjects)-1 and pos2<len(listofExonObjects):
        exonA=listofExonObjects[pos1]
        exonB=listofExonObjects[pos2]
        startA=exonA.exonSTART; stopA=exonA.exonSTOP
        startB=exonB.exonSTART; stopB=exonB.exonSTOP
        if startA <= startB <= stopA:
            if stopB >= stopA: #exons offset from each other
                exonB.exonSTART=startA
                listofExonObjects.remove(exonA)
                pos2=pos1+1
            else: #exonB is "engulfed by" exonA
                assert stopB <= stopA
                exonB.exonSTART=startA
                exonB.exonSTOP=stopA
                listofExonObjects.remove(exonA)
                pos2=pos1+1
        elif startB <= startA <= stopB:
            if stopB >= stopA: #exonA is "engulfed by" exonB
                listofExonObjects.remove(exonA)
                pos2=pos1+1
            else: #exons offset from each other
                assert stopB <=stopA
                exonB.exonSTOP=exonA.exonSTOP
                listofExonObjects.remove(exonA)
                pos2=pos1+1
        else:#They don't overlap and were not combined
            #ChangePositionAfterNoCombo
            length=len(listofExonObjects)
            if pos2 < length-1 and pos1 < length-2:
                pos2=pos2+1
            elif pos2==length-1 and pos1<length-2:
                pos1=pos1+1
                pos2=pos1+1
            else:
                assert pos2==length-1 and pos1==length-2, ("ERROR:"
                        +" pos2 is actually "+`pos2`+" and pos1 is "+`pos1`)
                pos1=length #stop the while loop
    #Update the linc's start, stop, and exon properties to reflect
    #the collapsed exon list
    newlincstart=[]
    newlincstop=[]
    newlincexons=[]
    exoncount=1
    listofExonObjects.sort(key=lambda x: x.exonSTART) #sort the exons by start
    #if you don't sort them by start position, then the Genome Browser can't
    #open your BED file because the Genome Browser assumes your
    #exons are provided in order :P
    for ekzon in listofExonObjects:
        newlincstart.append(ekzon.exonSTART)
        newlincstop.append(ekzon.exonSTOP)
        newlincexons.append(exoncount)
        exoncount+=1
    linc.start=newlincstart
    linc.stop=newlincstop
    linc.exons=newlincexons

def makeSmallfeatures_fromSIGOVA(sigovafilename):
    """Return a list of SmallFeature objects based on the Sigova data file.
    
    NOTES:
    The Sigova data file is in a format similar to BED but not identical. 
    NOTE about the coordinate system used in Sigova data:
    They do not specify in the paper what coordinate system they are using,
    but they do say this in thier supporting information: "We first mapped the
    genomic origin of each transcript within a locus using BLAT"
    I went to look up BLAT, and BLAT uses the same coordinate system as UCSC
    (the same format used in BED files.) Therefore I did the same conversions
    that I do on BED files in order to put the Sigova data into the one-based
    inclusive system I use internally for the rest of the files that go through
    definelincs.py.
    SEE: http://genome.ucsc.edu/goldenPath/help/blatSpec.html
    BLAT output format is psl, and psl format uses coordinate system of UCSC."""
    sigovafile=open(sigovafilename,'r')
    featslistsigova=[]
    for line in sigovafile:
        lineaslist=line.rsplit(); chromstr=lineaslist[0]
        if "X" in chromstr:
                chromset=23 #for X
        elif "Y" in chromstr:
            chromset=24 #for Y
        else: #chr is an autosome
            try:
                chromset=int(chromstr[3:])
            except:
                try: #this is for the weird formatting of "chr#_stuff"
                    underscoreindex=chromstr.find("_")
                    chromset=int(chromstr[3:underscoreindex])
                except:
                    chromset=0 #This includes all chrUn
        ranger=[int(lineaslist[1])+1,int(lineaslist[2])] #plus one on the start
        #to make it one-based
        overallstart=int(lineaslist[1]) #this is a zero-based start
        overallstop=int(lineaslist[2]) #needed for sigova file only
        featurename=lineaslist[3]
        strand=lineaslist[5]
        if "." in strand:
            strand="."
        if "+" in strand:
            strand="+"
        if "-" in strand:
            strand="-"
        startslist=[overallstart+1] #add one to make it one-based included
        stopslist=[overallstop] #end position is already one-based included
        exonlist=[]
        for index in range(len(startslist)):
            exonlist.append(index+1)
        newfeature=FeatureClass_Small.SmallFeature(featurename,strand,
                            chromset,startslist,stopslist,exonlist,ranger,[])
        featslistsigova.append(newfeature)
    return featslistsigova

#===============================================================================
#----------CLASS----------------------------------------------------------------
#===============================================================================
class Exon(object):
    exonSTART=0
    exonSTOP=0
    def __init__(self,exonSTARTpos,exonSTOPpos):
        self.exonSTART=exonSTARTpos
        self.exonSTOP=exonSTOPpos
    def __str__(self):
        return "Exon "+`self.exonSTART`+"-"+`self.exonSTOP`

#===============================================================================
#----------EXECUTION------------------------------------------------------------
#===============================================================================
if __name__ == '__main__':
    #Parse config file & check that it is validly filled out
    configfile=open(sys.argv[1],'r')
    datasets={}
    outputdir=""
    use=[]
    for line in configfile:
        if line[0]!="#":
            lineaslist=line.rsplit()
            if lineaslist[0]=="broad8000:":
                datasets["broad8000"]=lineaslist[1]
            elif lineaslist[0]=="gencode:":
                datasets["gencode"]=lineaslist[1]
            elif lineaslist[0]=="sigova0:":
                datasets["sigova0"]=lineaslist[1]
            elif lineaslist[0]=="hangauerS3:":
                datasets["hangauerS3"]=lineaslist[1]
            elif lineaslist[0]=="pseudogenes:":
                datasets["pseudogenes"]=lineaslist[1]
            elif lineaslist[0]=="refgene:":
                datasets["refgene"]=lineaslist[1]
            elif lineaslist[0]=="use:":
                use=lineaslist[1:]
            elif lineaslist[0]=="output:":
                outputdir=lineaslist[1]
            else:
                pass
    assert len(datasets)>0, ("Error:"
            +" not enough information in the configuration file.")
    assert os.path.isdir(outputdir), ("Error: "
                                +outputdir+" is not an existing directory.")
    for name in ["broad8000","gencode","sigova0","hangauerS3"]:
        if name not in use:
            del datasets[name]
    for key in datasets:
        pathtofile=datasets[key]
        assert os.path.isfile(pathtofile), ("Error: "+pathtofile
                                            +" is not an existing regular file")
    #Run module
    use.sort()
    define_lincRNAs(datasets,outputdir)

