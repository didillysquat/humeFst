class sample(object):

    def __init__(self, name, compComplement, hostTaxon = None, region = None, reef = None, totalSeqs = None):
        self.name = name
        self.compComplement = compComplement
        self.hostTaxon = hostTaxon
        self.region = region
        self.reef = reef
        self.cladeCollectionList = []
        self.totalSeqs = totalSeqs
        self.finalTypeCladeCollectionList = []
        self.intraAbundanceDict = self.createIntraAbundanceDict()


    def createIntraAbundanceDict(self):
        abundanceDict = {}
        for OCCURENCE in self.compComplement.listOfits2SequenceOccurances:
            abundanceDict[OCCURENCE.name] = OCCURENCE.abundance
        return abundanceDict

    def __str__(self):
        return self.name

    def addCladeCollection(self, cladeCollection):
        self.cladeCollectionList.append(cladeCollection)


class its2SequenceOccurance:

    def __init__(self, name, sequence, abundance, clade):
        self.name = name
        self.sequence = sequence
        self.abundance = abundance
        self.clade = clade


class cladeCollection:

    def __init__(self, clade, cutoff, listOfits2SequenceOccurances, foundWithinSample, proportion):
        self.clade = clade
        self.cutoff = cutoff
        self.listOfits2SequenceOccurances = listOfits2SequenceOccurances
        self.initialType = None
        self.footPrint = frozenset([a.name for a in listOfits2SequenceOccurances])
        self.maj = listOfits2SequenceOccurances[0].name
        self.foundWithinSample = foundWithinSample
        self.cladalProportion = proportion


    def addInitialType(self, symbiodiniumType):
        self.initialType = symbiodiniumType


class completeComplement:

    def __init__(self, listOfits2SequenceOccurances):
        self.listOfits2SequenceOccurances = listOfits2SequenceOccurances