
import config
import statistics
import re
# We should consider re-naming this a sampleBasedSymbiodiniumType or something similar vs a type in the symbiodiniumTypeDB
class symbiodiniumType:
# when creating the final types try to pass in the total seqs from sample and the list of occurences from sample.compcomplement.listofits2collection
    def __init__(self, footPrint, clade, maj, typeOfType, coDom, listofoccurences = None, name =None, sorteddefiningits2occurances=None, listofSamples=None, majList = None,  totalseqs=None,  typeSupport=None, totalSeqs=None, permute=None, abundancelist = None):
        self.coDom = coDom
        if coDom :
            self.maj = {listofSamples[i]: majList[i] for i in range(len(listofSamples))}
        else:
            self.maj = maj
        self.majList = majList
        self.listOfSamples = listofSamples
        self.footPrint = footPrint # This should be a frozenset
        self.clade = clade
        self.typeOfType = typeOfType
        if name == None and sorteddefiningits2occurances == None:
            self.name, self.sortedDefiningIts2Occurances = self.createSymbiodiniumTypeName()
        else:
            self.name = name
            self.sortedDefiningIts2Occurances = sorteddefiningits2occurances
    def __str__(self):
        return self.name

    def createSortedDefiningIts2Occurances(self, listofoccurences, totalseqs):
        proportionDict = {}
        cladalProportion = sum([occurence.abundance for occurence in listofoccurences if occurence.clade == self.clade]) / totalseqs  # This is the decimal percentage of the proportion of that clades sequences in that sample
        for OCCURENCE in listofoccurences:
            if OCCURENCE.clade == self.clade:
                proportionDict[OCCURENCE.name] = OCCURENCE.abundance / (cladalProportion * totalseqs)  # This is the decimal percentage of this sequence as a proportion of it's clades sequences within the sample
        sortedList = [(a[0], a[1]) for a in sorted(proportionDict.items(), key=lambda x: x[1], reverse=True) if a[0] in self.footPrint]
        return sortedList, proportionDict

    def createSymbiodiniumTypeName(self,  listOfOccurences = None, totalSeqs = None): # Outputname
        # Create an abundance dictionary for the sequences in the footprint across all samples that contain that type
        proportionDict = {seq: 0 for seq in self.footPrint}
        for SAMPLEKEY in self.listOfSamples:
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                for occurance in [occur for occur in CLADECOLLECTION.listOfSeqsAboveCutOff if occur.name in self.footPrint]:
                    proportionDict[occurance.name] = proportionDict[occurance.name] + occurance.abundance


        # A sorted list of the intra abundances as calculated across all samples containing the type
        sortedList = [(a[0], a[1]) for a in sorted(proportionDict.items(), key=lambda x: x[1], reverse=True)]
        copyOfSortedList = list(sortedList)
        # A sorted list of the names of the seqs by highest abundance
        sortedList = [a[0] for a in sortedList]


        # typeName
        added = []
        if self.coDom:
            sortedcoDomList = [item for item in sortedList if item in set(self.majList)] # Need this as cant pass directly through the coDomDict.keys() as these are not sorted and can't parse through sortedList directly as some of them may not be in the coDOmdict and so will through error at the conditional
            namePart1 = '/'.join([CLJ(codomintra) for codomintra in sortedcoDomList]) # Add any coDom intras first
            added.extend([codomintra for codomintra in sortedcoDomList])
            namePart2 = '-'.join([CLJ(noncoDomIntras) for noncoDomIntras in sortedList if noncoDomIntras not in added]) # If it isn't already in the name because it is a codom then add in order of abundance within the sample
            if len(namePart2) > 0: # Only if there is something in name Part 2
                typeName = '-'.join([namePart1, namePart2])
            else:
                typeName = namePart1
        else:
            typeName = '-'.join([CLJ(noncoDomIntras) for noncoDomIntras in sortedList if noncoDomIntras not in added])
        #Delete these two lines if you want to go back to creating names for the final and inital types instead of assigning final types the inital types name without taking
        # into consideration any changes in the intra abundances
        return typeName, copyOfSortedList



    def makeCoDomDict(self, clade, listofsamples):
        sampleToMajDict = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name in listofsamples:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                        sampleToMajDict[SAMPLE.name] =  sorted(CLADECOLLECTION.listOfSeqsAboveCutOff, key=lambda x: x.abundance, reverse=True)[0].name

        return sampleToMajDict
#Check LaJ Covnertion
def CLJ(VAR):
    try:
        return config.oursToLaJDict[VAR]
    except:
        return VAR

class finalTypeCladeCollection:

    def __init__(self, foundWithinSample, clade, cutoff, listOfFinalTypes):
        self.foundWithinSample = foundWithinSample
        self.clade = clade
        self.cutoff = cutoff
        self.sortedListOfFinalTypes = self.generateSortedListOfFinalTypes(listOfFinalTypes, foundWithinSample) # This should be sorted according to the abundance of the Type in the sample

    def generateSortedListOfFinalTypes(self, listOfFinalTypes, foundWithinSample): # I think this might be wrong, we need a list where the most abundant type is first, i.e. addition of all it's
        compCompDict = {}
        SAMPLE = config.abundanceList[foundWithinSample]
        for OCCURENCE in SAMPLE.compComplement.listOfits2SequenceOccurances:
            compCompDict[OCCURENCE.name] = OCCURENCE.abundance
        # Here we have the abundance dict where key = intra in sample and value = abundance
        typeAbundDict = {}
        for TYPE in listOfFinalTypes:
            for INTRANAME in config.typeDB[TYPE].footPrint:
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

    def generateIntrasInfoFinalForAllTypes(self):
        for types in self.keys():
            self[types].generateIntrasInfoFinal()

    def addType(self, symbiodiniumType):
        self[symbiodiniumType.name] = symboidiniumDBTypeEntry(name=symbiodiniumType.name, clade=symbiodiniumType.clade,
                                                              samplename=symbiodiniumType.listOfSamples,
                                                              codom=symbiodiniumType.coDom,
                                                              typeOfType=symbiodiniumType.typeOfType,
                                                              maj=symbiodiniumType.majList, majList = symbiodiniumType.majList, footprint=symbiodiniumType.footPrint, sorteddefiningits2occurances=symbiodiniumType.sortedDefiningIts2Occurances)

    # def initialiseFromAbundanceList(self, abundanceList):
    #     for SAMPLE in abundanceList:
    #         for CLADECOLLECTION in SAMPLE.cladeCollectionList:
    #             # Type instance
    #             TI = CLADECOLLECTION.initialType
    #             if TI.name not in self.keys():
    #                 # If this is the first time we come across this type we intialise it to the type DB with an inital type instance
    #                 self[TI.name] = symboidiniumDBTypeEntry(name=TI.name,
    #                                                                 clade=TI.clade,
    #                                                                 codom=TI.coDom,
    #                                                                 samplename=SAMPLE.name,
    #                                                                 typeOfType=TI.typeOfType,
    #                                                                 maj=TI.maj, majList=TI.name.majList)
    #             elif TI.name in self.keys():
    #                 self[TI.name].update(typeoftype=TI.typeOfType, samplename=SAMPLE.name, coDom=TI.coDom,
    #                                               listofdefiningintras=TI.sortedDefiningIts2Occurances, maj=TI.maj)
    #                 #Here we need to upgrade the type entry rather than start a new one
    #         for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
    #             for FINALTYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
    #                 if FINALTYPE.name not in self.keys():
    #                     self[FINALTYPE.name] = symboidiniumDBTypeEntry(name=FINALTYPE.name,
    #                                                                      clade=FINALTYPE.clade,
    #                                                                      codom=FINALTYPE.coDom,
    #                                                                      samplename=SAMPLE.name, maj=FINALTYPE.maj,
    #                                                                      typeOfType=FINALTYPE.typeOfType, majList=FINALTYPE.name.majList)
    #                     # Here we need to initialise an entry into the DB using a final type to intialise
    #                 else:
    #                     # Here we need to upgrade the type with a FINAL type info, i.e. add to which samples
    #                     # found in (final) but not add anything to the defining intras
    #                     self[FINALTYPE.name].update(typeoftype=FINALTYPE.typeOfType, samplename=SAMPLE.name)
    #     # This method will initialise the database from a given abundanceList
    #     # Alternatively the addType and __init__ can be used to add types to the DB

class symboidiniumDBTypeEntry:

    '''This creates an instance of the typeEntry probably with only a single instance of the type found in a sample
    We will continue to update the information as we come across instances of the type within the samples'''
    def __init__(self, name, clade, codom, samplename, typeOfType, maj, majList, footprint, sorteddefiningits2occurances):
        ''' If type is final then we take no listofdefiningintras as we only want this info from intial cases'''
        print('Initialising type: {0}'.format(name))
        self.name = name
        self.coDom = codom
        self.clade = clade
        self.majList = majList
        self.footPrint = footprint
        self.sortedDefiningIts2Occurances = sorteddefiningits2occurances


        if self.coDom == False:
            self.majDict = {maj[0]: len(maj)}
        else:

            self.majDict = {samplename[i]: maj[i] for i in range(len(samplename))}




        self.samplesFoundInAsFinal = []
        self.updateSamplesFoundIn(samplename) # Initialises self.samplesFoundInAsFinal and self.definingIntras

    def __str__(self):
        return self.name

    def generateIntrasInfoFinal(self):
        '''
        This will create a self.intrasInfoFinal variable that will be exactly the same as the self.definingIntrasInfo
        but created form the intra info found in the samples in the self.samplesFoundInAsFinal list
        :return: This should initiate a self.intrasInfoFinal
        '''
        listOfIntrasInOrderOfAbun = [a[0] for a in self.sortedDefiningIts2Occurances]
        self.intrasInfoFinal = [[[], 0, 0 , []] for intra in listOfIntrasInOrderOfAbun]


        if len(self.footPrint) > 1:
            for SAMPLEKEY in self.samplesFoundInAsFinal:
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == self.clade:
                        totSeqs = sum([SAMPLE.intraAbundanceDict[intra] for intra in listOfIntrasInOrderOfAbun])
                        for i in range(len(self.footPrint)):
                            self.intrasInfoFinal[i][0].append(SAMPLE.intraAbundanceDict[listOfIntrasInOrderOfAbun[i]]/totSeqs)
                        # Here we add the ratios info
                        # In the second list of the self.definingIntrasInfo, we divide the abundance of the given intra
                        # by the abundance of the most abundant intra so as to get a ratio
                        # The ratios of the first intra will always be 1 as we are dividing by itself
                        for i in range(len(self.footPrint)):
                            self.intrasInfoFinal[i][3].append(self.intrasInfoFinal[i][0][-1]/self.intrasInfoFinal[0][0][-1])

            for i in range(len(listOfIntrasInOrderOfAbun)):
                self.intrasInfoFinal[i][1] = sum(self.intrasInfoFinal[i][0])/len(self.intrasInfoFinal[i][0])
                self.intrasInfoFinal[i][2] = statistics.stdev(self.intrasInfoFinal[i][0])


    def updateSamplesFoundIn(self, additionalSamplesFoundIn):
        #Here we will update the self.definingIntras parameter so that we
        # check through the current list of samples and work out the average abundance and s.d. for each
        # of the intras in the type
        listOfIntrasInOrderOfAbun = [a[0] for a in self.sortedDefiningIts2Occurances]
        try:
            self.samplesFoundInAsInitial = self.samplesFoundInAsInitial + additionalSamplesFoundIn
        except:
            self.definingIntrasInfo = [[[], 0, 0 , []] for intra in listOfIntrasInOrderOfAbun]  # A list of tuples where each tuple holds the mean and S.d. for the intra (represented by position)

            self.samplesFoundInAsInitial = additionalSamplesFoundIn

        if len(self.footPrint) > 1:
            for SAMPLEKEY in additionalSamplesFoundIn:
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == self.clade:
                        totSeqs = sum([SAMPLE.intraAbundanceDict[intra] for intra in listOfIntrasInOrderOfAbun])
                        for i in range(len(self.footPrint)):
                            self.definingIntrasInfo[i][0].append(SAMPLE.intraAbundanceDict[listOfIntrasInOrderOfAbun[i]]/totSeqs)
                        # Here we add the ratios info
                        # In the second list of the self.definingIntrasInfo, we divide the abundance of the given intra
                        # by the abundance of the most abundant intra so as to get a ratio
                        # The ratios of the first intra will always be 1 as we are dividing by itself
                        for i in range(len(self.footPrint)):
                            self.definingIntrasInfo[i][3].append(self.definingIntrasInfo[i][0][-1]/self.definingIntrasInfo[0][0][-1])

            for i in range(len(listOfIntrasInOrderOfAbun)):
                self.definingIntrasInfo[i][1] = sum(self.definingIntrasInfo[i][0])/len(self.definingIntrasInfo[i][0])
                self.definingIntrasInfo[i][2] = statistics.stdev(self.definingIntrasInfo[i][0])

        return



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


