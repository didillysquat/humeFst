import HumeFst
import config


class symbiodiniumType:
# when creating the final types try to pass in the total seqs from sample and the list of occurences from sample.compcomplement.listofits2collection
    def __init__(self, footPrint, clade, maj, typeOfType, coDom, listofoccurences = None, listofSamples=None,  supportedType = None, coDomDict = None,  totalseqs=None,  typeSupport=None, totalSeqs=None):
        self.coDom = coDom
        # If coDom then cannot have a singlemaj
        # We will condition off the setting of the self.maj property conditional on the type being not coDom and see where we fail.
        # In it's place we will make a dictionary of which samples (of this initial type) have which Maj
        # We will see if the machine breaks if the final type doesn't have a maj
        if coDom and typeOfType == 'INITIAL':
            self.maj = self.makeCoDomDict(clade = clade, listofsamples = listofSamples)
        else:
            self.maj = maj
        self.supportedType = supportedType
        self.listOfSamples = listofSamples
        self.footPrint = footPrint # This should be a frozenset
        self.coDomDict = coDomDict
        self.clade = clade
        self.typeOfType = typeOfType
        if self.typeOfType == 'INITIAL':
            self.name, self.sortedDefiningIts2Occurances = self.createSymbiodiniumTypeNameClass()
        elif self.typeOfType == 'FINAL':
            self.name, self.sortedDefiningIts2Occurances = self.createSymbiodiniumTypeNameClass(listOfOccurences=listofoccurences, totalSeqs=totalSeqs, typeSupport=typeSupport)


    def createSymbiodiniumTypeNameClass(self,  listOfOccurences = None, totalSeqs = None, typeSupport = None): # Outputname

        if self.typeOfType == 'INITIAL': # create an abundance dictionary based on all occurances of the types in all sample
            ProportionDict = {seq: 0 for seq in self.footPrint}
            for SAMPLE in config.abundanceList:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.footPrint == self.footPrint:
                        for occur in CLADECOLLECTION.listOfits2SequenceOccurances:
                            ProportionDict[occur.name] = ProportionDict[occur.name] + occur.abundance
            sortedList = [(a[0], a[1]) for a in sorted(ProportionDict.items(), key=lambda x: x[1], reverse=True)]
            copyOfSortedList = list(sortedList)
            sortedList = [a[0] for a in sortedList]


        else: #self.typeOfType == 'FINAL': # create an abundance dictionary based on only occurances from the sample in question
            ProportionDict = {}
            cladalProportion = sum([occurence.abundance for occurence in listOfOccurences if occurence.clade == self.clade])/totalSeqs
            for OCCURENCE in listOfOccurences:
                if OCCURENCE.clade == self.clade:
                    ProportionDict[OCCURENCE.name] = OCCURENCE.abundance/(cladalProportion*totalSeqs) # This is the decimal percentage of the proportion of that clades sequences in that sample
            sortedList = [(a[0], a[1]) for a in sorted(ProportionDict.items(), key=lambda x: x[1], reverse=True) if a[0] in self.footPrint]
            copyOfSortedList = list(sortedList)
            sortedList = [a[0] for a in sortedList]
            self.typeTotalProportion = sum([ProportionDict[intra] for intra in self.footPrint])

        footprintList = list(self.footPrint)
        #convert  the footPrint Intras, the ProportionDict and the coDomdict to LaJeunesse names
        i = 0
        while i < len(footprintList):
            if footprintList[i] in config.oursToLaJDict.keys():
                if self.coDom:
                    if footprintList[i] in self.coDomDict.keys():
                        self.coDomDict[config.oursToLaJDict[footprintList[i]]] = self.coDomDict[footprintList[i]]
                        del self.coDomDict[footprintList[i]]
                if footprintList[i] in sortedList:
                    j = 0
                    while j < len(sortedList):
                        if sortedList[j] == footprintList[i]:
                            sortedList[j] = config.oursToLaJDict[footprintList[i]]
                            break
                        j += 1
                footprintList[i] = config.oursToLaJDict[footprintList[i]]
            i += 1

        # typeName
        added = []
        if self.coDom:
            sortedcoDomList = [item for item in sortedList if item in self.coDomDict.keys()] # Need this as cant pass directly through the coDomDict.keys() as these are not sorted and can't parse through sortedList directly as some of them may not be in the coDOmdict and so will through error at the conditional
            namePart1 = '/'.join([codomintra for codomintra in sortedcoDomList if self.coDomDict[codomintra]]) # Add any coDom intras first
            added.extend([codomintra for codomintra in sortedcoDomList if self.coDomDict[codomintra]])
        namePart2 = '-'.join([noncoDomIntras  for noncoDomIntras in sortedList if noncoDomIntras not in added]) # If it isn't already in the name because it is a codom then add in order of abundance within the sample
        if self.coDom:
            if len(namePart2) > 0: # Only if there is something in name Part 2
                typeName = '-'.join([namePart1, namePart2])
            else:
                typeName = namePart1
        else:
            typeName = namePart2


        return typeName, copyOfSortedList

    def makeCoDomDict(self, clade, listofsamples):
        sampleToMajDict = {}
        for SAMPLE in config.abundanceList:
            if SAMPLE.name in listofsamples:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                        sampleToMajDict[SAMPLE] =  sorted(CLADECOLLECTION.listOfits2SequenceOccurances, key=lambda x: x.abundance, reverse=True)[0].name

        # majToSampleListDict = {}
        # for MAJ in set(sampleToMajDict.values()):
        #     tempList = []
        #     for KEY in sampleToMajDict.keys():
        #         if sampleToMajDict[KEY] == MAJ:
        #             tempList.append(KEY)
        #     majToSampleListDict[MAJ] = tempList

        return sampleToMajDict

class finalTypeCladeCollection:

    def __init__(self, foundWithinSample, clade, cutoff, listOfFinalTypes):
        self.foundWithinSample = foundWithinSample
        self.clade = clade
        self.cutoff = cutoff
        self.listOfFinalTypes = self.sortedListOfFinalTypes(listOfFinalTypes) # This should be sorted according to the abundance of the Type in the sample
        # self.typeBasedCompCollection = [] # A list of its2SequenceOccurances that only contain intras that have been found in supported types identified in this sample
        self.identified = False # True if all of the intras of the given clade above the given cutoff are found in the types' defining footprints
        self.maj  = self.listOfFinalTypes[0].maj# The most abundant intra in the most abundant type
        self.mostAbundantType = self.listOfFinalTypes[0].name



    def sortedListOfFinalTypes(self, listOfFinalTypes): # I think this might be wrong, we need a list where the most abundant type is first, i.e. addition of all it's
        compCompDict = {}
        for SAMPLE in config.abundanceList:
            if SAMPLE.name == self.foundWithinSample:
                for OCCURENCE in SAMPLE.compComplement.listOfits2SequenceOccurances:
                    compCompDict[OCCURENCE.name] = OCCURENCE.abundance
        # Here we have the abundance dict where key = intra in sample and value = abundance
        typeAbundDict = {}
        for TYPE in listOfFinalTypes:
            for INTRANAME in TYPE.footPrint:
                if TYPE in typeAbundDict.keys():
                    typeAbundDict[TYPE] = typeAbundDict[TYPE] + compCompDict[INTRANAME]
                else:
                    typeAbundDict[TYPE] = compCompDict[INTRANAME]
        orderedTypeList = [a[0] for a in sorted(typeAbundDict.items(), key=lambda x: x[1], reverse=True)]

        return orderedTypeList

    def containsCoDom(self):
        containsCoDom = False
        for TYPE in self.listOfFinalTypes:
            if TYPE.coDom:
                containsCoDom = True
        return containsCoDom

    def typeBasedCompCollection(self):
        tempList = []
        for TYPE in self.listOfFinalTypes:
            for INTRA in TYPE.footPrint:
                tempList.append(INTRA)

        return tempList





