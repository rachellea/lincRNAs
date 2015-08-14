#Rachel Ballantyne
#2013

"""Parse a BED file using the correct coordinate system.
See http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for a description of BED
file format.
The overall start position provided is zero-based. The overall stop position
provided is one-based. However, the exon sizes provided are calculated so that
when you add them to the exon starts provided, you get one-based ends.
Therefore, if you want your SmallFeature to have all of its coordinates in a
1-based system, you need to do the following IN THIS ORDER:
(1) Extract the chromStart (overall start position) (ZERO BASED)
(2) Extract the blockStarts (ZERO BASED)(exon start positions, which "should be
    calculated relative to chromStart")
(3) To get exon absolute start positions (ZERO BASED), add chromStart to each
    of the blockStarts
(4) To get exon absolute end positions (ONE BASED, due to how size is
    calculated) extract the blockSizes and then add each blockSize to the
    corresponding blockStart.
(5) To get exon absolute start positions (ONE BASED), add one to every
    blockStart
(6) To get exon chromStart (ONE BASED) add one to it.
Checks: after all this, chromStart should equal the start of the first exon,
and chromEnd should equal the end of the last exon. Furthermore, since you are
now one-based for the start and one-based for the end, the size of each exon
as listed in the BED file should be calculated by using
(onebased end)-(onebased start)+1.
EVERY FEATURE OBJECT CREATED FROM A BED FILE USING THIS MODULE WILL HAVE ALL
OF ITS POSITIONS PROVIDED IN ONE-BASED COORDINATES, WHERE START IS INCLUDED IN
THE FEATURE AND END IS INCLUDED IN THE FEATURE.

Example: for Hangauer_DatasetS3_FPKM1_hg19.bed, the first two lines of the file
contain the following information:
NAME                        chromStart  chromEnd  blockSizes       blockStarts
FPKM1_group_1_tr_1  89294     120932    2335,150,105,158,    0,2796,23405,31480,
FPKM1_group_2_tr_1  160444    161525    246,213,            0,868

chromStart is zero-based, chromEnd is one-based, blockStarts are relative to
chromStart and are calculated so that when you add a blockStart to a chromStart
you get a zero-based exon Start, and blockSizes are calculated so that when a
blockSize is added to a zero-based exon Start, it produces a one-based exon End.

After parsing these two lines with makeSmallfeatures_fromBed, two SmallFeature
objects should be created with the following ranger, start, stop, and exon
properties: (EVERY COORDINATE ONE-BASED)
FPKM1_group_1_transcript_1
ranger [89295,120932]
start  [89295,92091,112700,120775]
stop   [91629,92240,112804,120932]

FPKM1_group_2_transcript_1
ranger [160445,161525]
start  [160445,161313]
stop   [160690,161525]"""

from FeatureClass_Small import SmallFeature
import copy

def makeSmallfeatures_fromBED(bedfilename):
    """Return a list of SmallFeature objects based on a BED file.
    Each line of the BED file is turned into a SmallFeature object."""
    if "/" in bedfilename: #then a full path was provided
        bedfile=open(bedfilename,'r')
    else: #then assume what directory the file is in
        bedfile=open("C:/Users/raba/MyPrograms/shared_data_and_modules/"
                     +bedfilename,'r')
    featureslistbed=[]
    for line in bedfile:
        lineaslist=line.rsplit()
        chromstr=lineaslist[0]
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
        featurename=lineaslist[3]
        strand=lineaslist[5]
        #to remove spaces from strand:
        if "." in strand:
            strand="."
        if "+" in strand:
            strand="+"
        if "-" in strand:
            strand="-"
        #Obtain list of absolute start positions for exons (ZERO BASED)
        relativestartslist=(lineaslist[11]).rsplit(",") 
        startslist_zerobased=[]
        for startpoz in relativestartslist[:-1]:
            startslist_zerobased.append(int(startpoz)+overallstart)
        #Obtain list of absolute stop positions for exons (ONE BASED--due to
        #how size is provided)
        sizeslist=(lineaslist[10]).rsplit(",")
        stopslist=[]
        for index in range(len(startslist_zerobased)):
            correspondingsize=int(sizeslist[index])
            correspondingstart=startslist_zerobased[index]
            calculatedstop=correspondingsize+correspondingstart
            stopslist.append(calculatedstop)
        #Now you can make the start positions ONE BASED
        startslist_onebased=[]
        for item in startslist_zerobased:
            startslist_onebased.append(item+1)
        #Make exon list
        exonlist=[]
        for index in range(len(startslist_onebased)):
            exonlist.append(index+1) #So if your startslist has 3 members,
            #exonlist will be [1,2,3]
        #Make list of SmallFeature objects
        newfeature=SmallFeature(featurename,strand,chromset,
                            startslist_onebased,stopslist,exonlist,ranger,[])
        #TO DO: RANGER SHOULD NOT CHANGE WHEN YOU SET IT AGAIN!!!!
        featureslistbed.append(newfeature)
    return featureslistbed


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

def makeBED_fromchromlist(chromlist,filename):
    """Given a chromified list, make a BED file with the specified name.
    PARAMETERS:
    <chromlist> (the "chromified list") is a list of lists where each sub-list
    represents a chromosome and contains SmallFeature objects that are on that
    chromosome.
    Note: the SmallFeature objects should have their positions in the one-based
    closed interval coordinate system. The positions will be manipulated
    appropriately to put them into the half-open zero-based BED file coordinate
    system. There is no header in the BED file created by this function.
    Here is some info pasted from <http://genome.ucsc.edu/FAQ/FAQformat.html>
    The first three required BED fields are:
        1.chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or
          scaffold (e.g. scaffold10671). 
        2.chromStart - The starting position of the feature in the chromosome
          or scaffold. The first base in a chromosome is numbered 0. 
        3.chromEnd - The ending position of the feature in the chromosome or
          scaffold. The chromEnd base is not included in the display of the
          feature. For example, the first 100 bases of a chromosome are defined
          as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
    The 9 additional optional BED fields are:
        4.name - Defines the name of the BED line. This label is displayed to
          the left of the BED line in the Genome Browser window when the track
          is open to full display mode or directly to the left of the item in
          pack mode. 
        5.score - A score between 0 and 1000. If the track line useScore
          attribute is set to 1 for this annotation data set, the score value
          will determine the level of gray in which this feature is displayed
          (higher numbers = darker gray). This table shows the Genome Browser's
          translation of BED score values into shades of gray
        6.strand - Defines the strand - either '+' or '-'. 
        7.thickStart - The starting position at which the feature is drawn
          thickly (for example, the start codon in gene displays). 
        8.thickEnd - The ending position at which the feature is drawn thickly
          (for example, the stop codon in gene displays). 
        9.itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track
          line itemRgb attribute is set to "On", this RBG value will determine
          the display color of the data contained in this BED line. NOTE: It is
          recommended that a simple color scheme (eight colors or less) be used
          with this attribute to avoid overwhelming the color resources of the
          Genome Browser and your Internet browser. 
        10.blockCount - The number of blocks (exons) in the BED line. 
        11.blockSizes - A comma-separated list of the block sizes. The number
          of items in this list should correspond to blockCount. 
        12.blockStarts - A comma-separated list of block starts. All of the
          blockStart positions should be calculated relative to chromStart.
          The number of items in this list should correspond to blockCount."""
    assert type(filename) is str, "Error filename must be a string"
    BEDfile=open(filename,'w')
    for chrom in chromlist:
        for feat in chrom:
            BEDfile.write("chr"+`feat.chromosome`+"\t") #chrom
            #chromStart, zero-based, start included:
            BEDfile.write(`feat.ranger[0]-1`+"\t")
            #chromEnd, zero based, end excluded:
            BEDfile.write(`feat.ranger[1]`+"\t")
            BEDfile.write(feat.featurename+"\t")        #name
            BEDfile.write("500\t")                      #score
            if feat.strand=="%" or feat.strand=="*":
                BEDfile.write(".\t")  #strand. Loss of partial strand info
                #since official BED format doesn't use my weird symbols % and *
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
