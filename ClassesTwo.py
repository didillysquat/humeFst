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
        self.cladalProportions = self.initCladProps()


    def initCladProps(self):
        cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        counterDict = {letter: 0 for letter in cladeList}
        for occur in self.compComplement:
            counterDict[occur.clade] += occur.abundance

        propDict = {letter: 0 for letter in cladeList}
        for clade in cladeList:
            propDict[clade] = counterDict[clade]/self.totalSeqs

        return propDict

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

    def __init__(self, clade, cutoff, listofseqsabovecutoff, foundwithinsample, cladalproportion):
        self.clade = clade # The clade in question that all seqs are from
        self.cutoff = cutoff # The cutoff percentage that the sequence abundance must be higher than to make the listOfSeqsAboveCutOff list and to be considered as the footprint
        self.listOfSeqsAboveCutOff = listofseqsabovecutoff # List of the sequences from the sample that were found as a proportion of the total sequences from the given clade above the cutoff
        self.initialType = None
        self.footPrint = frozenset([a.name for a in listofseqsabovecutoff]) # The footprint of this cladeCollection as defined by the presence of sequences above the cutoff
        try:
            self.maj = listofseqsabovecutoff[0].name # The predominant ITS2 sequence of the given clade
        except:
            self.maj = 'None'
        self.foundWithinSample = foundwithinsample # Which sample this clade collection is from
        self.cladalProportion = cladalproportion # The proportion of the sample's total seqs that are of the clade in question


    def addInitialType(self, symbiodiniumType):
        self.initialType = symbiodiniumType


class completeComplement:

    def __init__(self, listOfits2SequenceOccurances):
        self.listOfits2SequenceOccurances = listOfits2SequenceOccurances