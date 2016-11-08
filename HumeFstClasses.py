import HumeFst
import config

# We should consider re-naming this a sampleBasedSymbiodiniumType or something similar vs a type in the symbiodiniumTypeDB
class symbiodiniumType:
# when creating the final types try to pass in the total seqs from sample and the list of occurences from sample.compcomplement.listofits2collection
    def __init__(self, footPrint, clade, maj, typeOfType, coDom, listofoccurences = None, name =None, listofSamples=None, listofcodommajs = None, majList = None,  totalseqs=None,  typeSupport=None, totalSeqs=None, permute=None, abundancelist = None):
        self.coDom = coDom
        # If coDom then cannot have a singlemaj
        # We will condition off the setting of the self.maj property conditional on the type being not coDom and see where we fail.
        # In it's place we will make a dictionary of which samples (of this initial type) have which Maj
        # We will see if the machine breaks if the final type doesn't have a maj
        if coDom and typeOfType == 'INITIAL':
            self.maj = {listofSamples[i]: majList[i] for i in range(len(listofSamples))}

        else:
            self.maj = maj

        # self.supportedType = supportedType
        self.majList = majList
        self.listOfSamples = listofSamples
        self.footPrint = footPrint # This should be a frozenset
        self.listofcodommajs = listofcodommajs
        self.clade = clade
        self.typeOfType = typeOfType
        # Here we have a slight problem as the names of types might change slightly between initial and final type creations due to the abundances of the defining intras changing.
        # This causes us problems futher on in the software where for example we are trying to look up a final type in both the initial support and final support dictionaries
        # A1-Otu36794-Otu37876 may have become A1-Otu37876-Otu36794
        # I think one solution to this would be to simply stick with the names as created when the inital type was created.
        # I have implemented this below
        if self.typeOfType == 'INITIAL':
            self.name = self.createSymbiodiniumTypeNameClass()
            # This will be a lists of tuples giving the proportion that the defining sequences of the type represent as part of the sequences from their clade within the sample.
            # Same as for the FINAL types. Only difference being that this will be initialized when the newly identified types are added to the sequences they are found in.
            self.sortedDefiningIts2Occurances = None
        #elif self.typeOfType == 'FINAL':
        #    self.name, self.sortedDefiningIts2Occurances = self.createSymbiodiniumTypeNameClass(listOfOccurences=listofoccurences, totalSeqs=totalSeqs, typeSupport=typeSupport)
        elif self.typeOfType == 'FINAL':
            self.name = name
            # These sortedDefininITS2Occurances lists are lists of tuples giving the proportion that the defining sequences of the type represent as part of the sequences from their clade within the sample.
            # e.g. if there are 600 C and 400 D sequences and the footprint for the type in question is made up of two D intras of 10 and 80 % proportion of the clade D sequences then their proportions would be .1 and .8.
            self.sortedDefiningIts2Occurances = self.createSymbiodiniumTypeNameClass(listOfOccurences=listofoccurences, totalSeqs=totalSeqs)
        if coDom and typeOfType == 'FINAL':
            self.maj = self.sortedDefiningIts2Occurances[0][0]


    def createSortedDefiningIts2Occurances(self, listofoccurences, totalseqs):
        proportionDict = {}
        cladalProportion = sum([occurence.abundance for occurence in listofoccurences if occurence.clade == self.clade]) / totalseqs  # This is the decimal percentage of the proportion of that clades sequences in that sample
        for OCCURENCE in listofoccurences:
            if OCCURENCE.clade == self.clade:
                proportionDict[OCCURENCE.name] = OCCURENCE.abundance / (cladalProportion * totalseqs)  # This is the decimal percentage of this sequence as a proportion of it's clades sequences within the sample
        sortedList = [(a[0], a[1]) for a in sorted(proportionDict.items(), key=lambda x: x[1], reverse=True) if a[0] in self.footPrint]
        return sortedList, proportionDict

    def createSymbiodiniumTypeNameClass(self,  listOfOccurences = None, totalSeqs = None): # Outputname

        # Create an abundance dictionary for the sequences in the footprint.
        # This will count how many times the sequences in the footprint are found in all samples that contain the footprint
        if self.typeOfType == 'INITIAL':
            ProportionDict = {seq: 0 for seq in self.footPrint}
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.footPrint == self.footPrint:
                        for occur in CLADECOLLECTION.listOfSeqsAboveCutOff:
                            ProportionDict[occur.name] = ProportionDict[occur.name] + occur.abundance
            # A sorted list of items with the highest abundnace first
            sortedList = [(a[0], a[1]) for a in sorted(ProportionDict.items(), key=lambda x: x[1], reverse=True)]
            copyOfSortedList = list(sortedList)
            # A sorted list of the names of the seqs by highest abundance
            sortedList = [a[0] for a in sortedList]


        else: #self.typeOfType == 'FINAL': # create an abundance dictionary based on only occurances from the sample in question
            sortedList, ProportionDict = self.createSortedDefiningIts2Occurances(listOfOccurences, totalSeqs)
            copyOfSortedList = list(sortedList)
            sortedList = [a[0] for a in sortedList]
            self.typeTotalProportion = sum([ProportionDict[intra] for intra in self.footPrint])


        footprintList = list(self.footPrint)
        #convert  the footPrint Intras, the ProportionDict and the coDomdict to LaJeunesse names
        i = 0
        while i < len(footprintList):
            if footprintList[i] in config.oursToLaJDict.keys():
                if self.coDom:
                    if footprintList[i] in self.listofcodommajs:
                        for n, seq in enumerate(self.listofcodommajs):
                            if seq == footprintList[i]:
                                self.listofcodommajs[n] = config.oursToLaJDict[footprintList[i]]
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
            sortedcoDomList = [item for item in sortedList if item in self.listofcodommajs] # Need this as cant pass directly through the coDomDict.keys() as these are not sorted and can't parse through sortedList directly as some of them may not be in the coDOmdict and so will through error at the conditional
            namePart1 = '/'.join([codomintra for codomintra in sortedcoDomList]) # Add any coDom intras first
            added.extend([codomintra for codomintra in sortedcoDomList])
        namePart2 = '-'.join([noncoDomIntras for noncoDomIntras in sortedList if noncoDomIntras not in added]) # If it isn't already in the name because it is a codom then add in order of abundance within the sample
        if self.coDom:
            if len(namePart2) > 0: # Only if there is something in name Part 2
                typeName = '-'.join([namePart1, namePart2])
            else:
                typeName = namePart1
        else:
            typeName = namePart2
        #Delete these two lines if you want to go back to creating names for the final and inital types instead of assigning final types the inital types name without taking
        # into consideration any changes in the intra abundances
        if self.typeOfType == 'FINAL':
            return copyOfSortedList
        else:
            return typeName

    def makeCoDomDict(self, clade, listofsamples):
        sampleToMajDict = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name in listofsamples:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                        sampleToMajDict[SAMPLE.name] =  sorted(CLADECOLLECTION.listOfSeqsAboveCutOff, key=lambda x: x.abundance, reverse=True)[0].name

        return sampleToMajDict

class finalTypeCladeCollection:

    def __init__(self, foundWithinSample, clade, cutoff, listOfFinalTypes):
        self.foundWithinSample = foundWithinSample
        self.clade = clade
        self.cutoff = cutoff
        self.listOfFinalTypes = self.sortedListOfFinalTypes(listOfFinalTypes) # This should be sorted according to the abundance of the Type in the sample

        self.identified = False # True if all of the intras of the given clade above the given cutoff are found in the types' defining footprints
        self.isMixedIdentification = False # If we have two final footprints that share intras but one is not a subset of the other this is TRUE
        self.maj  = self.listOfFinalTypes[0].maj # The most abundant intra in the most abundant type
        self.mostAbundantType = self.listOfFinalTypes[0].name



    def sortedListOfFinalTypes(self, listOfFinalTypes): # I think this might be wrong, we need a list where the most abundant type is first, i.e. addition of all it's
        compCompDict = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
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

    # Return a list of seq names that are found in all types in the listOfFinalTypes
    def typeBasedCompCollection(self):
        tempList = []
        for TYPE in self.listOfFinalTypes:
            for INTRA in TYPE.footPrint:
                tempList.append(INTRA)

        return tempList

class symbiodiniumTypeDB(dict):

    def __init__(self, *arg, **kwargs):
        super(symbiodiniumTypeDB, self).__init__(*arg, **kwargs)

    def addType(self, symbiodiniumType):
        self[symbiodiniumType.name] = symboidiniumDBTypeEntry(name=symbiodiniumType.name, clade=symbiodiniumType.clade,
                                                              samplename=symbiodiniumType.listOfSamples,
                                                              codom=symbiodiniumType.coDom,
                                                              typeOfType=symbiodiniumType.typeOfType,
                                                              maj=symbiodiniumType.majList, majList = symbiodiniumType.majList)

    def initialiseFromAbundanceList(self, abundanceList):
        for SAMPLE in abundanceList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                # Type instance
                TI = CLADECOLLECTION.initialType
                if TI.name not in self.keys():
                    # If this is the first time we come across this type we intialise it to the type DB with an inital type instance
                    self[TI.name] = symboidiniumDBTypeEntry(name=TI.name,
                                                                    clade=TI.clade,
                                                                    codom=TI.coDom,
                                                                    samplename=SAMPLE.name,
                                                                    typeOfType=TI.typeOfType,
                                                                    maj=TI.maj, majList=TI.name.majList)
                elif TI.name in self.keys():
                    self[TI.name].update(typeoftype=TI.typeOfType, samplename=SAMPLE.name, coDom=TI.coDom,
                                                  listofdefiningintras=TI.sortedDefiningIts2Occurances, maj=TI.maj)
                    #Here we need to upgrade the type entry rather than start a new one
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                for FINALTYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
                    if FINALTYPE.name not in self.keys():
                        self[FINALTYPE.name] = symboidiniumDBTypeEntry(name=FINALTYPE.name,
                                                                         clade=FINALTYPE.clade,
                                                                         codom=FINALTYPE.coDom,
                                                                         samplename=SAMPLE.name, maj=FINALTYPE.maj,
                                                                         typeOfType=FINALTYPE.typeOfType, majList=FINALTYPE.name.majList)
                        # Here we need to initialise an entry into the DB using a final type to intialise
                    else:
                        # Here we need to upgrade the type with a FINAL type info, i.e. add to which samples
                        # found in (final) but not add anything to the defining intras
                        self[FINALTYPE.name].update(typeoftype=FINALTYPE.typeOfType, samplename=SAMPLE.name)
        # This method will initialise the database from a given abundanceList
        # Alternatively the addType and __init__ can be used to add types to the DB




class symboidiniumDBTypeEntry:

    '''This creates an instance of the typeEntry probably with only a single instance of the type found in a sample
    We will continue to update the information as we come across instances of the type within the samples'''
    def __init__(self, name, clade, codom, samplename, typeOfType, maj, majList):
        ''' If type is final then we take no listofdefiningintras as we only want this info from intial cases'''
        print('Initialising type: {0}'.format(name))
        self.name = name
        self.coDom = codom
        self.clade = clade
        self.majList = majList
        # If listofdefiningintras == None then set to empty
        self.definingIntras = self.calculateDefiningIntrasInfo()
        if self.coDom == False:
            self.majDict = {maj[0]: 1}
        else:
            if typeOfType == 'INITIAL':
                self.majDict = {samplename[i]: maj[i] for i in range(len(samplename))}
            # We don't count the majs for final type allocations
            elif typeOfType == 'FINAL':
                self.majDict = {}
        if typeOfType == 'FINAL':
            self.samplesFoundInAsFinal = [samplename]
            self.samplesFoundInAsInitial = []
        elif typeOfType == 'INITIAL':
            self.samplesFoundInAsFinal = []
            if len(samplename) > 1:
                self.samplesFoundInAsInitial = samplename
            else:
                self.samplesFoundInAsInitial = [samplename]


    def calculateDefiningIntrasInfo(self):
        a=5

    def update(self, typeoftype, samplename, coDom = None, listofdefiningintras = None, maj = None):
        ''' We only take a maj argument if this is coDom and initial
        we only take listofdefiningintras argument if this typeoftype is inital'''
        if typeoftype == 'INITIAL':
            self.samplesFoundInAsInitial.append(samplename)
            self.definingIntras.update(anotherlistofdefiningintras=listofdefiningintras)
            if coDom:
                self.majDict[samplename] = maj[samplename]

        # We don't count the majs for final type allocation
        elif typeoftype == 'FINAL':
            self.samplesFoundInAsFinal.append(samplename)




    def initialiseSymTypeEntry(self):
        # TODO write this method that will cycle through config.abundance to get the above info
        a = 5




class definingIntraSet:
    '''This will be a set of defiing intras for a given type. We will only consider abundance of
    intras within types that are initial (not final). We will use the abundances, and standard deviations from
    these abundances as well as the ratio between the intras to determine whether types identified in final types
    are indeed representative of this given type'''
    def __init__(self, listofdefiningintras):
        if listofdefiningintras != None:
            self.definingIntrasSet = [listofdefiningintras]
        else:
            self.definingIntrasSet = []


    def definingInfo(self):
        a = 5
        # This will give us a data set of mean abundances and s.d.s of each of the defining intas and possibly all seqs
        # found in all the samples containing this type as an initial type. Used in LDA linear discriminant anlysis

    def update(self, anotherlistofdefiningintras):
        self.definingIntrasSet.append(anotherlistofdefiningintras)

class definingIntra:
    ''' This is a given intragenomic sequence within the context of a type
    we are only considering its abundance within initial (not final) instances of the type
    '''
    def __init__(self, typeName, intraName):
        self.listOfIntraInstances


