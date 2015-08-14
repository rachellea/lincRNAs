#GTFparser_general.py
#Rachel Ballantyne

"""GTF parser. GTF files use a one-based included start 
and a one-based included end for the coordinate system."""

import copy
from FeatureClass_Small import SmallFeature

def makeSmallfeatures_fromGTF(gtffile,compress,filterby):
    """Read a GTF and return a list of SmallFeature objects based on the GTF.
    PARAMETERS:
    <gtffile> is the path to a GTF file.
    <compress> is either True or False. If <compress> is True, then consecutive
        features of the same name will be put into one SmallFeature object. If
        <compress> is False, then every line of the GTF file will be made into
        a separate SmallFeature object. You almost always want <compress> to be
        True.
    <filterby> is either "ALL" or a string that specifies which gene type to
        keep. If <filterby> is not "ALL" then only lines in the GTF file with
        a gene type matching the chosen string will be kept."""
    f=open(gtffile,'r')
    featureslist=[]
    for line in f:
        if line[0]=="#":
            pass #ignore the descriptive lines at the beginning of the file
        else:
            #Extract the needed information from the line of the GTF file
            lineaslist=line.rsplit("\t")
            chromstr=lineaslist[0]
            if "chrX" in chromstr:
                chromset=23 #for X
            elif "chrY" in chromstr:
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
            startpos=int(lineaslist[3])
            stoppos=int(lineaslist[4])
            strand=lineaslist[6]
            moreinfo=lineaslist[8]
            #PARSING moreinfo
            #fname
            fname_dex=moreinfo.find("gene_name")
            if fname_dex==-1:
                fname_dex2=moreinfo.find("gene_id")
                if fname_dex2==-1:
                    fname="empty"
                else:
                    g2=moreinfo[fname_dex2:]
                    colon2=g2.find(";")
                    fname=g2[9:colon2-1] #if you can't use gene_name use gene_id
            else:
                g1=moreinfo[fname_dex:]
                colon=g1.find(";")
                fname=g1[11:colon-1]
            #exonnumber
            exonnumber_dex=moreinfo.find("exon_number")
            if exonnumber_dex==-1:
                exonnumber=0
            else:
                g2=moreinfo[exonnumber_dex:]
                colon=g2.find(";")
                if "lincRNAs" in gtffile:
                    exonnumber=int(g2[13:colon-1]) #the exon numbers have
                    #quotes around them in Broad data
                if "gencode" in gtffile:
                    exonnumber=int(g2[12:colon]) #the exon numbers do NOT have
                    #quotes around them in gencode data
            #genetype - needed for filtering
            genetype_dex=moreinfo.find("gene_type")
            if genetype_dex==-1:
                genetype="empty"
            else:
                g1=moreinfo[genetype_dex:]
                colon=g1.find(";")
                genetype=g1[11:colon-1]
            #Make a list of SmallFeature objects
            if compress and (filterby=="ALL" or genetype==filterby):
                if featureslist!=[]:
                    prevfeat=featureslist[-1]
                if ((featureslist==[]) or fname!=prevfeat.featurename
                    or strand!=prevfeat.strand
                    or chromset!=prevfeat.chromosome):
                    #if the featureslist is empty or if the current feature
                    #is not part of the previous feature, then make a new
                    #feature
                    feat=SmallFeature(fname,strand,chromset,[startpos],
                                      [stoppos],[exonnumber])
                    featureslist.append(feat)
                else: #if the current feature IS part of the previous feature,
                    #modify the previous feature accordingly
                    modyfeat=featureslist[-1]
                    modyfeat.start=modyfeat.start+[startpos]
                    modyfeat.stop=modyfeat.stop+[stoppos]
                    modyfeat.exons=modyfeat.exons+[exonnumber]
            else: #compress is False; you want every feature as its own object
                if filterby=="ALL" or genetype==filterby:
                    feat=SmallFeature(fname,strand,chromset,[startpos],
                                      [stoppos],[exonnumber])
                    featureslist.append(feat)
    return featureslist

#===============================================================================
#------------------Additional Functions Also Present in BEDparser---------------
#===============================================================================

def giveranger(featurelist):
    """Add information to the ranger of each SmallFeature object.
    PARAMETERS:
    <featurelist> is a list of SmallFeature objects; for example, it could be
        the output of makeSmallfeatures_fromGTF() which is a
    
    NOTES:
    The ranger property is normally empty, unless giveranger is called on the
    featurelist. This function extracts the lowest start position and the
    highest stop position to create a ranger for the feature that encompasses
    the entire feature. Notice how I make a deep copy of the start and stop
    lists to prevent the originals from being modified (if I sorted
    feature.start and sorted feature.stop within the function, it would also
    sort them in the featurelist, and it would basically ruin the usefulness of
    these properties completely since the start and stop positions would no
    longer correlate with each other according to order in which they were
    appended.)"""
    for feature in featurelist:
        x=feature.start
        startlist=copy.deepcopy(x)
        y=feature.stop
        stoplist=copy.deepcopy(y)
        startlist.sort()
        stoplist.sort()
        firststart=startlist[0]
        laststop=stoplist[-1]
        feature.ranger=[firststart,laststop]

def expandranger(featurelist,bases):
    """Expand the ranger of each SmallFeature in <featurelist> by <bases>.
    Specifically, <bases> is subtracted from the start in the ranger and added
    to the stop in the ranger.
    PARAMETERS:
    <featurelist> is a list of SmallFeature objects.
    <bases> is an int between 0 and 20000"""
    assert type(bases) is int, "bases must be an integer"
    assert 0<=bases<=20000, "bases must be between zero and 20000"
    for feature in featurelist:
        firststart=feature.ranger[0]
        laststop=feature.ranger[1]
        newstart=firststart-bases
        if newstart < 0:
            newstart=0
        newstop=laststop+bases
        feature.ranger=[newstart,newstop]

def sortSmallFeatsbychrom(smallfeaturelist):
    """Return a list of SmallFeature objects sorted into sub-lists by chrom.
    Specifically, take the <smallfeaturelist> and put each SmallFeature object
    into a sub-list depending on which chromosome it is on, and return a
    two-dimensional list of the form
    [[chrom1],[chrom2],[chrom3]...[chrom23=X][chrom24=Y]] where each sub-list
    contains the SmallFeature objects that are on the relevant chromosome."""
    featsbychrom=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
        [],[],[],[]]
    for feat in smallfeaturelist:
        if feat.chromosome!=0: #to ignore all those RNAs that have unknown
            #chromosomal locations and were given 0 as their chrom number
            featsbychrom[(feat.chromosome)-1].append(feat)
    return featsbychrom

def countchromifiedlist(chromifiedlist):
    """Return the number of SmallFeatures in a two-dimensional list created by
    sortSmallFeatsbychrom()."""
    counter=0
    for chromo in chromifiedlist:
        for feat in chromo:
            counter=counter+1
    return counter
