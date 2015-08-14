#Rachel Ballantyne
#2013

class SmallFeature(object):
    """A SmallFeature Object has the following properties:
        featurename: string indicating the name or ID of the feature
        strand: string. Can be "+","-", or "." Initialized to "." (no strand
            info)
        chromosome: integer indicating the chromosome number, 1 - 22 (for
            autosome), 23 (for X), 24 (for Y)
        start: list of genomic start locations, which are ints (1-based)
        stop: list of genomic stop locations, which are ints (1-based)
        exons: list of ints representing exon numbers. These correspond to the
            start/stop location provided in the corresponding position of the
            start/stop lists
        ranger: This is based off the start and stop properties.
            When parsing GTF files using GTFparser_general() with compress set
            to True, a SmallFeature contains many different start and stop
            locations for various components of that SmallFeature (e.g. exons).
            The giveranger() function can be used to get the earliest start
            position and the latest stop position and save that as a 2-member
            list in the ranger property for that feature. It's only called
            ranger instead of range because "range" is a Python keyword.
            ranger is initialized to [0,0]
        misc: a field that can contain whatever you want (useful in various
            functions). Initialized to []
        source: a string that indicates where this SmallFeature was derived
            from."""
    #FIELDS
    _featurename = ""
    _strand = "."
    _chromosome = 0
    _start = [0]
    _stop = [0]
    _exons = [0]
    _ranger = [0,0]
    _misc = []
    _source = ""
    
    @property
    def featurename(self):
        return self._featurename
    @featurename.setter
    def featurename(self,value):
        assert type(value) is str, "Error: cannot set featurename to "+`value`
        self._featurename=value
    
    @property
    def strand(self):
        return self._strand
    @strand.setter
    def strand(self,value):
        assert type(value) is str
        assert (value=="." or value=="+"
                or value=="-" or value=="*"
                or value=="%"), "Error: cannot set strand to "+value
        self._strand=value
    
    @property
    def chromosome(self):
        return self._chromosome
    @chromosome.setter
    def chromosome(self, value):
        assert type(value) is int, ("Error: "
                            +"chromosome must be an int; cannot be "+`value`)
        self._chromosome = value
    
    @property
    def start(self):
        return self._start
    @start.setter
    def start(self,value):
        assert type(value) is list
        for item in value:
            assert type(item) is int
        self._start=value
    
    @property
    def stop(self):
        return self._stop
    @stop.setter
    def stop(self,value):
        assert type(value) is list
        for item in value:
            assert type(item) is int
        self._stop=value
    
    @property
    def exons(self):
        return self._exons
    @exons.setter
    def exons(self,value):
        assert type(value) is list
        for item in value:
            assert type(item) is int
        self._exons=value
    
    @property
    def ranger(self):
        return self._ranger
    @ranger.setter
    def ranger(self,value):
        assert type(value) is list and len(value)==2, ("Error:"
                                        +"ranger must be a two-member list")
        self._ranger=value
    
    @property
    def misc(self):
        return self._misc
    @misc.setter
    def misc(self,value):
        self._misc=value
    
    @property
    def source(self):
        return self._source
    @source.setter
    def source(self,value):
        assert type(value) is str, "source must be a string"
        self._source=value
    
    #METHODS
    def __init__(self, featurename,strand,chromosome,start,stop,exons,
                 ranger=[0,0], misc=[], source=""):
        self.featurename=featurename
        self.strand=strand
        self.chromosome=chromosome
        self.start=start
        self.stop=stop
        self.exons=exons
        self.ranger=ranger
        self.misc=misc
        self.source=source

    def __str__(self):
        return ("\n\n"+self.featurename+"\nchr "+`self.chromosome`+", "
                +`self.ranger[0]`+"-"+`self.ranger[1]`+", strand: "+self.strand
                +"\nstart: "+`self.start`+"\nstop: "+`self.stop`+""
                +"\nexons: "+`self.exons`+"\nmisc: "+`self.misc`+"\nsource: "
                +self.source)
    
    def __eq__(self,otherobject):
        return self.__dict__ == otherobject.__dict__
