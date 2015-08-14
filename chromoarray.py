#Rachel Ballantyne
#2013

"""This module is used for creating an array of chars for a chromosome where
each char corresponds to one base of that chromosome and indicates which of the
following five categories that base belongs to:
m=mRNA only
l=lincRNA only
b=mRNA and lincRNA
i=intergenic
n=nonintergenic
All coordinates that are dealt with in this module are one-based inclusive
start, one-based inclusive stop. See document
"11-18-13_figuring_out_coordinate_systems.docx" for more information.
The module uses arrays only, which takes less memory than an implementation
that uses lists as intermediates."""

import cPickle
import array

def create_22chrom_arrays(mRNAlocslist,linclocslist,nonintlocslist,lmbases,
                          whichchroms=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                                       17,18,19,20,21,22],
                          userealsizes=True,dirtosavein=""):
    """Create an array for each autosome specifying the classification of each
    base in that chromosome.
    
    Required arguments:
        The inputs <linclocslist>, <mRNAlocslist>, and <nonintlocslist> are all
        3D lists with the following structure:
            -->>Verbal description: a list of 22 lists representing the
                autosomes; each autosome
                contains "[#,#] location  lists" that specify the intervals
                where the indicated features
                are located (e.g. where lincRNAs are located.)
            -->>Visual description:
                [[*],[*],[*],[*],[*],[*],[*],[*],[*],[*],[*],[*],[*],[*],[*],
                [*],[*],[*],[*],[*],[*],[*]]
                where an example of * is [1,2],[11,15],[22,23]
            Total size of each chromosome (for hg19/GRCh37.p13) was obtained
            from the Genome Reference Consortium:
            http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/
            human/data/index.shtml
            (NOTE: the Wikipedia numbers in the Chromosome article are wrong!)
        <lmbases> is an int >= 0 that specifies the number of bases by which
        you want to expand lincRNAs and mRNAs in either direction.
    
    Optional arguments:
    <whichchroms> is a list of ints, which can be any ints between 1 and 22,
        specifying which autosomes should have arrays created for them.
        Default is to create an array for every autosome.
    <userealsizes> is a bool. If True, then the actual (huge) chromosome sizes
        will be used. If False, fake chromosome sizes of 300 will be used
        (this is for testing purposes.)
    <dirtosavein> is a string specifying which directory to save in. It will
        be attached to the beginning of the file name so that the chromosome
        arrays will be saved inside that dir.
    
    Output: The array for each chromosome will be saved under the filename
    'chromosome#_ext<lmbases>_array.bin'
    """
    assert type(whichchroms) is list, ("Error: "
                    +"whichchroms must be a list, not "+`type(whichchroms)`)
    if userealsizes:
        chrsizedict={"chr1":249250621,"chr2":243199373,"chr3":198022430,
                "chr4":191154276,"chr5":180915260,"chr6":171115067,
                "chr7":159138663,"chr8":146364022,"chr9":141213431,
                "chr10":135534747,"chr11":135006516,"chr12":133851895,
                "chr13":115169878,"chr14":107349540,"chr15":102531392,
                "chr16":90354753,"chr17":81195210,"chr18":78077248,
                "chr19":59128983,"chr20":63025520,"chr21":48129895,
                "chr22":51304566}
    else:
        chrsizedict={"chr1":300,"chr2":300,"chr3":300,"chr4":300,
                     "chr5":300,"chr6":300,"chr7":300,"chr8":300,
                     "chr9":300,"chr10":300,"chr11":300,"chr12":300,
                     "chr13":300,"chr14":300,"chr15":300,"chr16":300,
                     "chr17":300,"chr18":300,"chr19":300,"chr20":300,
                     "chr21":300,"chr22":300}
    whichchromdexes=[]
    for number in whichchroms:
        whichchromdexes.append(number-1)
    for chromdex in whichchromdexes: #CHANGE THIS TO "IN RANGE 1" FOR TESTING
        chromnumber=chromdex+1 #add one since the output of range() starts at 0
        chromosize=chrsizedict["chr"+`chromnumber`]
        if dirtosavein != "":
            dirtosavein=dirtosavein+"/"
        makechromarray(dirtosavein+"chromosome"+`chromnumber`,chromosize,
                       mRNAlocslist[chromdex],linclocslist[chromdex],
                       nonintlocslist[chromdex],lmbases)


def makechromarray(arraynameprefix,chromsize,mRNAlist,lincRNAlist,nonlist,
                   lmbases):
    """Save an array of chars for the specified chromosome using cPickle.
    First create a list of chars. Then modify element of the list appropriately
    to correctly classify the base. Finally an array is created from the list
    and saved using cPickle.
    
    Preconditions:
    <arraynameprefix> is a string that will become the beginning of the filename
        for the pickled saved array. Example: "chromosome14"
    <chromsize> is an int>0 that indicates the size of the chromosome in bases. 
    <mRNAlist>, <lincRNAlist>, <nonlist>:  The three input lists must be
        lists of two-member lists where each two-member list is an interval
        denoting presence of the specified feature. E.g. if
        mRNAlist=[[20,30],[45,50]] that means there is an mRNA from 20 to 30
        (including 20 and 30) and another one from 45 to 50 (including 45 and
        50). Default value for a base in the chromosome is "intergenic".
        See NOTE_74356 for info about the "both" classification
    <lmbases> is an int >= 0 that specifies the number of bases by which you
        want to expand lincRNAs and mRNAs in either direction.
    
    Output: The array of chars will be saved under the file name
    '<arraynameprefix>_ext<lmbases>_array.bin'"""
    #Checks before doing stuff:
    assert type(arraynameprefix) is str,("Error: "
        +"arraynameprefix must be a string, not "+str(type(arraynameprefix)))
    assert type(chromsize) is int and chromsize>0, ("Error:"
        +" chromsize must be an int>0, not "+str(chromsize))
    for interval in mRNAlist: #check to ensure that no intervals start with 0,
        #since these intervals should all be 1-based start 1-based end inclusive
        if not(interval[0]): #if the number in the 0th position of the interval
            #is zero, that's BAD
            raise MyZeroStartError
    for interval in lincRNAlist:
        if not(interval[0]):
            raise MyZeroStartError
    for interval in nonlist:
        if not(interval[0]):
            raise MyZeroStartError
    #Doing stuff:
    chromarray=array.array('c','') #initialize empty array
    for number in xrange(chromsize): #if you use range instead of xrange you
        #consume immense amounts of memory.
        chromarray.append("i")
    assignletter(chromarray,mRNAlist,{"i":"m","m":"m"},lmbases) #NOTE_2309
    assignletter(chromarray,lincRNAlist,{"i":"l","m":"b","l":"l","b":"b"},lmbases)
    assignletter(chromarray,nonlist,{"i":"n","m":"m","b":"b","l":"l","n":"n"},0)
    #the only values that get changed on the basis of "nonlist" are i --> n
    #(lincs stay linc, mRNAs stay mRNA, etc., because
    #lincs and mRNAs are actually subsets of "nonintergenic" and we don't
    #want to classify everything as nonintergenic)
    cPickle.dump(chromarray,open(arraynameprefix+"_ext"+`lmbases`+"_array.bin",'wb'),2)
    print "Done with makechromarray() for "+arraynameprefix


def assignletter(bigarray,intervallist,letterdict,basez):
    """Use the <intervallist> to assign a char specified by <letterdict> to the
    appropriate intervals in <bigarray>.
    This is a helper function for makechromarray()
    Arguments:
    <bigarray> is an array of chars
    <intervallist> is a two-dimensional list of intervals as described in
    makechromarray() above The intervallist must not specify any interval
    endpoint that is greater than the length of bigarray.
    <letterdict> is a dictionary describing a mapping of one char to another
    char. The keys and values are both chars. There must be exactly one key
    for every char that appears in bigarray, and each key must have exactly
    one value (also a char). The <intervallist> specifies the intervals
    within <bigarray> where the char in bigarray
    must be changed to a different char as specified by <letterdict>.
    The number 1 in an interval corresponds to the 0th position of bigarray.
    E.g. if bigarray is ['a','b','a','a','b','b','c','a'],
    intervallist is [[1,2],[5,7]],
    and letterdict is {a:z,b:y,c:x} then assignletter() will modify
    <bigarray> so that it becomes ['z','y','a','a','y','y','x','a']
    <basez> is an int >= 0, the number of bases by which the intervals in the
    intervallist will be expanded in each direction (<basez> will be subtracted
    from the start of the interval and will be added to the end of the
    interval.)"""
    #THIS IS TRICKY YOU HAVE TO PRETEND THAT LISTS ARE NUMBERED STARTING AT ONE!
    #SO THE ZEROTH POSITION OF THE BIGARRAY IS SPECIFIED BY THE NUMBER 1
    #IN AN INTERVAL
    assert type(bigarray) is array.array and len(bigarray)>0, ("Error:"
                                +" bigarray must be an array of length > 0")
    assert type(basez) is int and basez >= 0, ("Error: "
        +"basez must be an int >= 0, not "+`basez`)
    for interval in intervallist:
        for index in xrange(interval[1]-interval[0]+1+(2*basez)): #NOTE_5563
            pos=max(interval[0]+index-1-basez,0) #NOTE_5564
            if basez>0 and pos > len(bigarray)-1: #NOTE_5565
                return
            charr=bigarray[pos]
            #modify bigarray according to letterdict:
            bigarray[pos]=letterdict[charr] 

#----------TESTING--------------------------------------------------------------
def _test_assignletter():
    """Test the function assignletter()"""
    bigarray1=array.array('c','aaaa')
    intervallist1=[[1,2]]
    letterdict1={"a":"x"}
    assignletter(bigarray1,intervallist1,letterdict1,0)
    assert bigarray1==array.array('c','xxaa'), "Error: Test1 failed"
    bigarray2=array.array('c','ixzp')
    intervallist2=[[1,1],[2,2],[3,3],[4,4]]
    letterdict2={'i':'I','x':'X','z':'Z','p':'P'}
    assignletter(bigarray2,intervallist2,letterdict2,0)
    assert bigarray2==array.array('c','IXZP'), ("Error:"
                +" Test2 failed. bigarray was "+`bigarray`+" instead of 'IXZP'")
    bigarray3=array.array('c','b')
    intervallist3=[[1,1]]
    letterdict3={'b':'q'}
    assignletter(bigarray3,intervallist3,letterdict3,0)
    assert bigarray3==array.array('c','q'), ("Error: "
                +"Test3 failed. bigarray was "+`bigarray`+" instead of 'q'")
    bigarray4=array.array('c','hello')
    intervallist4=[[1,3],[5,5]]
    letterdict4={'h':'w','e':'o','l':'r','o':'d'}
    assignletter(bigarray4,intervallist4,letterdict4,0)
    assert bigarray4==array.array('c','world'), ("Error: "
        +"Test4 failed. bigarray was "+`bigarray`
        +" instead of ['w','o','r','l','d']")
    print "Testing complete: assignletter() works properly"


def _test_makechromarray():
    """Test the function makechromarray()"""
    #CHECKING THAT ACTUAL OUTPUT MATCHES EXPECTED OUTPUT, for bases=0
    chromname1="chromtest1" #TEST1
    chromsize1=10
    mRNAlist1=[[1,1],[5,5]]
    lincRNAlist1=[[2,3]]
    nonlist1=[[9,9]]
    assert _is_makechromarray_right(chromname1,chromsize1,mRNAlist1,
                        lincRNAlist1,nonlist1,array.array('c','mllimiiini'),0)
    print "done with test1"
    assert _is_makechromarray_right("chromtest2",10,[[3,5],[9,10]],
                [[4,7]],[[5,5],[6,6],[7,8]],array.array('c','iimbbllnmm'),0)
    print "done with test2"
    assert _is_makechromarray_right("chromtest3",20,[[1,3],[6,8],[14,16]],
            [[3,4],[10,11]],[[4,4],[19,20]],
            array.array('c','mmblimmmilliimmmiinn'),0)
    print "done with test3"
    #TEST4: only mRNA list has entries
    assert _is_makechromarray_right("chromtest4",5,[[2,3]],[],[],
                                    array.array('c','immii'),0)
    print "done with test4"
    #TEST5: only lincRNA list has entries
    assert _is_makechromarray_right("chromtest5",5,[],[[3,4]],[],
                                    array.array('c','iilli'),0)
    print "done with test5"
    #TEST6: only nonlist has entries
    assert _is_makechromarray_right("chromtest6",5,[],[],[[1,2]],
                                    array.array('c','nniii'),0)
    print "done with test6"
    #TEST7: length of chromosome is 1
    assert _is_makechromarray_right("chromtest7",1,[[1,1]],[[1,1]],[[1,1]],
                                    array.array('c','b'),0)
    print "done with test7"
    #CHECKING THAT ERRORS ARISE WHEN THEY ARE SUPPOSED TO
    #TEST8: error if a starting index is zero for mRNA or lincRNA or nonlist
    try:
        _is_makechromarray_right("chromtest8a",5,[[0,3]],[],[],
                                 array.array('c','iiiii'),0)
        raise MyErrorNotCaughtError
    except MyZeroStartError:
        print ("test8a: "
               +"MyZeroStartError was caught for mRNA intervals list. Good.")
    try:
        _is_makechromarray_right("chromtest8b",5,[],[[0,2]],[],
                                 array.array('c','iiiii'),0)
        raise MyErrorNotCaughtError
    except MyZeroStartError:
        print ("test8b: "
               +"MyZeroStartError was caught for lincRNA intervals list. Good.")
    try:
        _is_makechromarray_right("chromtest8c",5,[],[],[[0,1]],
                                 array.array('c','iiiii'),0)
        raise MyErrorNotCaughtError
    except MyZeroStartError:
        print ("test8c: "
        +"MyZeroStartError was caught for nonintergenic intervals list. Good.")
    #TEST9: error if a stopping index is greater than the length of the chrom
    try:
        _is_makechromarray_right("chromtest9a",5,[[4,6]],[],[],
                                 array.array('c','iiiii'),0)
        raise MyErrorNotCaughtError
    except IndexError:
        print "test9a: IndexError was caught for mRNA intervals list. Good."
    try:
        _is_makechromarray_right("chromtest9b",5,[],[[3,7]],[],
                                 array.array('c','iiiii'),0)
        raise MyErrorNotCaughtError
    except IndexError:
        print "test9b: IndexError was caught for lincRNA intervals list. Good."
    try:
        _is_makechromarray_right("chromtest9c",5,[],[],[[4,8]],
                                 array.array('c','iiiii'),0)
        raise MyErrorNotCaughtError
    except IndexError:
        print ("test9c: "
            +"IndexError was caught for nonintergenic intervals list. Good.")
    #CHECKING THAT ACTUAL OUTPUT MATCHES EXPECTED OUTPUT, for various
    #bases parameters > 0
    assert _is_makechromarray_right("chromtest10",20,[[1,2],[19,20]],[[9,9]],
                        [[13,14]],array.array('c','mmmmmlllllllnnimmmmm'),3)
    print "done with test10"
    assert _is_makechromarray_right("chromtest11",21,[[9,11],[21,21]],[[1,2]],
                    [[1,2],[9,10]],array.array('c','lllliimmmmmmmiiiiimmm'),2)
    print "done with test11"
    assert _is_makechromarray_right("chromtest12",15,[[1,2],[7,8],[15,15]],
                    [[6,7]],[[13,13]],array.array('c','mmmilbbbmiiinmm'),1)
    print "done with test12"
    print "ALL TESTING COMPLETE"


def _is_makechromarray_right(chrnm,chrsiz,mRNAL,lincL,nonL,expectedoutput,
                             lmbayses):
    """Return True if the actual output (based on the specified parameters)
    matches the expected output"""
    makechromarray(chrnm,chrsiz,mRNAL,lincL,nonL,lmbayses)
    result=cPickle.load(open(chrnm+"_ext"+`lmbayses`+"_array.bin",'rb'))
    #print "result was "+`result`+" and the expected output was "+
    #`array.array('c',expectedoutputlist)`
    return result==expectedoutput

def _test_create_22chrom_arrays():
    """Test the function create_22chrom_arrays().
    NOTE: you must make one temporary change to create_22chrom_arrays() before
    running this test: Comment out chrsizedict and replace it with this line:
    chrsizedict={"chr1":10,"chr2":11,"chr3":12,"chr4":9}
    """
    print ("REMEMBER: you have to make two changes to create_22chrom_arrays()"
        +"before running this test.")
    mRNAlklist=[[[2,5],[9,9]],[[1,2],[10,11]],[[1,1],[5,5],[6,6]],[[1,9]],[],[],
        [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    linclklist=[[[4,6]],[[2,3],[5,5]],[[10,12]],[[1,9]],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[]]
    nonintlklist=[[[7,7],[8,8]],[[8,10]],[[3,4],[7,8]],[[1,9]],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],[],[]]
    create_22chrom_arrays(mRNAlklist,linclklist,nonintlklist,0,
                          whichchroms=[1,2,3,4])
    assert cPickle.load(
        open("chromosome1_ext0_array.bin",'rb'))==array.array('c','immbblnnmi')
    assert cPickle.load(
        open("chromosome2_ext0_array.bin",'rb'))==array.array('c','mbliliinnmm')
    assert cPickle.load(
        open("chromosome3_ext0_array.bin",'rb'))==array.array('c','minnmmnnilll')
    assert cPickle.load(
        open("chromosome4_ext0_array.bin",'rb'))==array.array('c','bbbbbbbbb')
    print ("PASSED TESTING"
        +" -- don't forget to reverse the changes you made "
        +"before running this function!!!")

#----------Custom Error Classes-------------------------------------------------
class MyZeroStartError(Exception):
    pass

class MyErrorNotCaughtError(Exception):
    """This is an error that can be raised during testing where the point was
    for Python to catch an error and yet somehow Python missed catching that
    error (thereby committing a 'MyErrorNotCaughtError')"""
    pass

