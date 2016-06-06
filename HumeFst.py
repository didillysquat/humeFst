# Imports
import config
import numpy as np
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.filedialog import askdirectory
import re
from tkinter import messagebox
import math
import itertools
import os
import multiprocessing
from HumeFstClasses import *
from ClassesTwo import *
import pickle
import config
from jinja2 import Environment, FileSystemLoader
from subprocess import call
import time


## THESE ARE ALL THE SUBFUNCTIONS ##
# read a file line by line with each line containing multiple items separated by a ','
def write2DListToFile(initialFile, defaultExtention, twoDlistToWrite):
    fileName = asksaveasfilename(initialfile=initialFile, defaultextension=defaultExtention)
    with open(fileName, mode='w+') as writer:
        i = 0
        while i < len(twoDlistToWrite):
            if i == len(twoDlistToWrite)-1:
                writer.write(','.join([str(a) for a in twoDlistToWrite[i]]))
            else:
                writer.write(','.join([str(a) for a in twoDlistToWrite[i]]) + '\n')
            i += 1

def write2DListToDestination(destination, twoDlistToWrite):
    with open(destination, mode='w+') as writer:
        i = 0
        while i < len(twoDlistToWrite):
            if i == len(twoDlistToWrite)-1:
                writer.write(','.join([str(a) for a in twoDlistToWrite[i]]))
            else:
                writer.write(','.join([str(a) for a in twoDlistToWrite[i]]) + '\n')
            i += 1

def readDefinedFileToList(filename):
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

# Write a list to a text file at a given location and name with each line being an item
def writeListToDestination(destination, listToWrite):
    print('Writing list to ' + destination);
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite)-1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite)-1:
                writer.write(listToWrite[i])
            i += 1

# Convert an interleaved fasta file into a normal fasta file. T
# The output of the Trex alignment is interleaved so needs converting
def convertInterleavedFasta(fastaFile):
    fastaout = []
    titleString = ''
    seqString = ''
    i = 0
    while i < len(fastaFile): # For each line of the interleaved fastafile
        if len(fastaFile[i]) > 0:
            if fastaFile[i][0] == '>': # Then this is the start of a new sequence
                if len(seqString) > 0:
                    fastaout.append(titleString)
                    fastaout.append(seqString)
                titleString = fastaFile[i]
                seqString = ''
            elif len(fastaFile[i]) > 0: # If a following line contains something then add to seqString, if empty then simply ignore
                seqString = seqString + fastaFile[i]
        if i == len(fastaFile)-1: # If we have just processed the last line then we need to add to fastaout
            fastaout.append(titleString)
            fastaout.append(seqString)
        i += 1
    return fastaout

# Converts a nucleotide matrix back to a fasta
def convertMatrixToFasta(fastaFile, matrix):
    i = 0
    j = 0
    while i < len(fastaFile):
        fastaFile[i+1] = ''.join(matrix[j])
        i += 2
        j += 1
    return fastaFile

def create2DmatrixfromFasta(afastafile):
    onlyseqs = []
    nucleotideMatrix = []
    # Create a matrix of nuleotides in form of 2D list
    i = 0
    while i < len(afastafile):
        onlyseqs.append(afastafile[i + 1])
        i += 2
    i = 0
    while i < len(onlyseqs):  # For each sequence
        templist = []
        j = 0
        while j < len(onlyseqs[i]):  # For each character in the sequence
            templist.append(onlyseqs[i][j])
            j += 1
        nucleotideMatrix.append(templist)
        i += 1
    return nucleotideMatrix

# create a new fasta alignment where any gap columns have been deleted
def removeGapsInFastaAlignment(gapfasta):
    fastamatrix = create2DmatrixfromFasta(gapfasta)
    fastamatrix = np.array(fastamatrix)
    seq = 0
    nuc = 0
    delete = True
    while nuc < (len(fastamatrix[seq])):  # For each column (nucleotide position)
        tempList = []
        delete = True
        while seq < len(fastamatrix):  # Go through each sequence
            tempList.append(fastamatrix[seq][nuc])
            seq +=1
        if all(p == '-' for p in tempList):
            fastamatrix = np.delete(fastamatrix, nuc, 1)
            nuc -= 1
            seq = 0
        nuc += 1
        seq = 0
    fastamatrix = fastamatrix.tolist()
    return convertMatrixToFasta(gapfasta, fastamatrix)

def CalculateSinglePairwiseFst(sampleoneits2occurances, sampletwoits2occurances, seqdistances):

    between = calcAvBtwGenDist(sampleoneits2occurances, sampletwoits2occurances, seqdistances)
    within = calcAvWthGenDist(sampleoneits2occurances, sampletwoits2occurances, seqdistances)
    if between != 0:
        return (between-within)/between
    else:
        return 0

def calcAvWthGenDist(sampleoneits2occurances, sampletwoits2occurances, seqdistances):

    runningTotalOne = 0.00
    runningTotalTwo = 0.00

    # When working out within seqs we can ignore comparisons between the same seqs so all we need to concentrate
    # on is comparison between different seqs. For each pair of seqs there will be N x M comparisons where N is
    # the abundance of the first seq and M is the abundance of the second seq so we merely need to multiply the
    # two abundances to gether and then multiply that result with the distance for that pair.

    # Similar to the between calculations we need to make sure that the logged abundances are re-normalised to 1000 for each sample
    # First calculate average for sampleone

    loggedTotalOne = sum([math.log(occurance.abundance) for occurance in sampleoneits2occurances])
    loggedTotalTwo = sum([math.log(occurance.abundance) for occurance in sampletwoits2occurances])

    for occuranceOne, occuranceTwo in itertools.combinations(sampleoneits2occurances, 2): # This gives us the combos that are possible. if only one type then this may fail but we can sort that easily
        runningTotalOne += (float((math.log(occuranceOne.abundance)/loggedTotalOne)*1000))*(float((math.log(occuranceTwo.abundance)/loggedTotalTwo)*1000))*float(seqdistances[frozenset({occuranceOne.name, occuranceTwo.name})])
    # Second calculate average for sampletwo
    for occuranceOne, occuranceTwo in itertools.combinations(sampletwoits2occurances, 2): # This gives us the combos that are possible. if only one type then this may fail but we can sort that easily
        runningTotalOne += (float((math.log(occuranceOne.abundance)/loggedTotalOne)*1000))*(float((math.log(occuranceTwo.abundance)/loggedTotalTwo)*1000))*float(seqdistances[frozenset({occuranceOne.name, occuranceTwo.name})])


    permutationsOne = 499500 # (n(n-1))/2 number of permutations
    permutationsTwo = 499500 # (n(n-1))/2 number of permutations
    return ((runningTotalOne/permutationsOne)+(runningTotalTwo/permutationsTwo))/2

def calcAvBtwGenDist(sampleoneits2occurances, sampletwoits2occurances, seqdistances):

    # When calculating the between seqs we will just use all combinations of uniqe sequences
    # i.e. the sampleonseqs and sampletwoseqs lists. we multiply the distance between the pair of seqs by the
    # result of multiplying the abundance of each of the seqs in each sample.

    # We need to consider how we will work with the abundances here. We want to log transform them so that the less abundant sequences have more power when it comes to calculating the Fst but at the same time we don't want to have different abundances in given samples just because of the different effects of logging
    # e.g. if we have a sample that has 1000 A1 that will log to 3, if we have one with 100 and 900 that will log to 2 + 2.954. It will then be as if the second sample has more abundance. So what we need to do is used the logged abundances to work out proportions of the original 1000. i.e. re-normalizeafter logging
    # so in the case of the second sequence, abundances would be (2/(2 + 2.954))*1000 and (2.954/(2 + 2.954))*1000. This way all samples will still have a total number of seqs that is 1000. THe only effect of the logging will be to have shifted the realtive proportions of the contained intras
    # In order to do this we should make a copy of the occurances lists and work out the new logs before using them in the calculations
    # In order to work out the proportional logged abundances we need to know the logged total for each of the occurances, then whenever we calculate abundance we can divide the logged abundnace by this and multiply by 1000
    loggedTotalOne = sum([math.log(occurance.abundance) for occurance in sampleoneits2occurances])
    loggedTotalTwo = sum([math.log(occurance.abundance) for occurance in sampletwoits2occurances])

    runningTotal = 0.00
    for oneOccurs in sampleoneits2occurances:
        for twoOccurs in sampletwoits2occurances:
            if oneOccurs.name != twoOccurs.name:#we only need to calculate distances if they are not the same
                runningTotal += seqdistances[frozenset({oneOccurs.name, twoOccurs.name})]*(float((math.log(oneOccurs.abundance)/loggedTotalOne)*1000))*(float((math.log(twoOccurs.abundance)/loggedTotalTwo)*1000))

    totalSeqCombos = 1000000

    return runningTotal/totalSeqCombos

def createMatrixFromColDists(Fstcoldistdict):
    FstMatrixCollection = []
    for CLADE in config.args.cladeList:
        listOfSamples = [] # [[cladeCollection.foundWithinSample for cladeCollection in sample.cladeCollectionList if cladeCollection.clade == CLADE] for sample in abundanceList] # All samples whether they have a type defined or not, then we can display the unsupported in R
        for SAMPLE in config.abundanceList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    listOfSamples.append(CLADECOLLECTION.foundWithinSample)
        # Create empty matrix that is of dimensions [rows(n samples + 1)][columns(number of samples + 1)]
        FstMatrix = [list(listOfSamples)]
        FstMatrix[0].insert(0, 'Samples')
        for sample in listOfSamples:
            t = ['N/A']*len(listOfSamples)
            t.insert(0, sample)
            FstMatrix.append(t)



        #Matrix is FstMatrix[row][column]
        # Fill matrix from the column Fsts in FstColDist
        row = 1
        while row < len(FstMatrix): # For each row
            col = 1
            while col < len(FstMatrix[row]): # For each col
                sampleOne = FstMatrix[0][col]
                sampleTwo = FstMatrix[row][0]
                if sampleOne == sampleTwo:
                    FstMatrix[row][col] = 0.00
                else:
                        FstMatrix[row][col] = FstMatrix[col][row] = Fstcoldistdict[config.args.cladeList.index(CLADE)][frozenset({sampleOne, sampleTwo})]
                col += 1
            row += 1
        FstMatrixCollection.append(FstMatrix)
    return FstMatrixCollection

def createFinalMatrixFromColDists(Fstcoldistdictsecond):
    FstMatrixCollection = []
    for CLADE in config.args.cladeList:
        listOfSamples = [] # [[cladeCollection.foundWithinSample for cladeCollection in sample.cladeCollectionList if cladeCollection.clade == CLADE] for sample in abundanceList] # All samples whether they have a type defined or not, then we can display the unsupported in R
        for SAMPLE in config.abundanceList:
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == CLADE and FINALTYPECLADECOLLECTION.identified:
                    listOfSamples.append(FINALTYPECLADECOLLECTION.foundWithinSample)
        # Create empty matrix that is of dimensions [rows(n samples + 1)][columns(number of samples + 1)]
        FstMatrix = [list(listOfSamples)]
        FstMatrix[0].insert(0, 'Samples')
        for sample in listOfSamples:
            t = ['N/A']*len(listOfSamples)
            t.insert(0, sample)
            FstMatrix.append(t)



        #Matrix is FstMatrix[row][column]
        # Fill matrix from the column Fsts in FstColDist
        row = 1
        while row < len(FstMatrix): # For each row
            col = 1
            while col < len(FstMatrix[row]): # For each col
                sampleOne = FstMatrix[0][col]
                sampleTwo = FstMatrix[row][0]
                if sampleOne == sampleTwo:
                    FstMatrix[row][col] = 0.00
                else:
                        FstMatrix[row][col] = FstMatrix[col][row] = Fstcoldistdictsecond[config.args.cladeList.index(CLADE)][frozenset({sampleOne, sampleTwo})]

                col += 1
            row += 1
        FstMatrixCollection.append(FstMatrix)
    return FstMatrixCollection

def writeByteObjectToDefinedDirectory(directory, objectString, object):
    f = open(directory + '\\' + objectString, 'wb+')
    pickle.dump(object, f)

def readByteObjectFromDefinedDirectory(directory, objectname):
    f = open(directory + '\\' + objectname, 'rb')
    return pickle.load(f)

def createFstColDists(cladeCollectionList, masterSeqDistancesDict):
    # Create a column formatted distance item for between samples
    FstColDist = []
    i = 0
    while i < len(cladeCollectionList): # For each clade
        counter = 0
        tempColDist = []
        seqProportionDict = {a[0]: a[1:] for a in cladeCollectionList[i]} # Create a dictionary key = sample name, value = sequence proportions
        for a,b in itertools.combinations([a[0] for a in cladeCollectionList[i]], 2): # For each paired combination of sample
            if a == b or seqProportionDict[a] == seqProportionDict[b]:
                tempColDist.append('\t'.join([str(a), str(b), str(0.00)]))
            elif a != b:
                tempColDist.append('\t'.join([str(a),str(b),str(CalculateSinglePairwiseFst(seqProportionDict[a], seqProportionDict[b], masterSeqDistancesDict))]))
            print(str(counter))
            counter += 1

        FstColDist.append(tempColDist)
        i += 1
    return FstColDist

def writeTypeBasedOutput(): # Produce all of the info such as number of Maj ITS2's per clade, number of clades, number of codoms, number of predefined types

    # THis had a problem in that it might show when there is increased support for a type but it doesn't work out when support has decreased for a type. For example, a lot of th
    # samples that had initial type A1 will have alternative final types. so the final support for A1 will drop when there is an alternative
    # In order to reflect this we will make an finalSupportOfInitialTypeDict which will keep track of the changes in support.
    # If a finaltypelist contains a type then this causes an increase of 1 in the support for the initial type.
    # If a finalType list has any types in it then the inital type of that cladecollection will get a -1 support as alterantives have been found.
    # If the finalType list is empty then that means that the support will be unchanged.
    # Unfortunately for the types that have only one defining intra, for example A1 they cannot gain support as you cannot have a finaltype of A1,
    # As your initial type would have had to have been A1 and therefore it wouldn't be in the final types list.

    # Create and populate the finalSupportOfInitialTypeDict
    # The values from this dict will then be used by adding them (some values will be negative) to the intial support to get the final support
    changeInFinalSupportOfInitialTypeDict = {}
    initialSupportOfInitialTypesDict = {} # This will be used in the next method key = footprint, value = initial support
    for SAMPLE in config.abundanceList:
        for CLADE in config.args.cladeList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.initialType.footPrint not in initialSupportOfInitialTypesDict:
                        initialSupportOfInitialTypesDict[CLADECOLLECTION.initialType.footPrint] = len(CLADECOLLECTION.initialType.listOfSamples)
                    for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                        if FINALTYPECLADECOLLECTION.clade == CLADE:
                            if len(FINALTYPECLADECOLLECTION.listOfFinalTypes) > 0: # If there is nothing in the listoffinaltypes then the final support for the initial doesn't change but we need to have this in the dictionary
                                if CLADECOLLECTION.initialType.footPrint in changeInFinalSupportOfInitialTypeDict.keys():
                                    changeInFinalSupportOfInitialTypeDict[CLADECOLLECTION.initialType.footPrint] = changeInFinalSupportOfInitialTypeDict[CLADECOLLECTION.initialType.footPrint] - 1
                                else:
                                    changeInFinalSupportOfInitialTypeDict[CLADECOLLECTION.initialType.footPrint] = - 1
                            else:
                                if CLADECOLLECTION.initialType.footPrint not in changeInFinalSupportOfInitialTypeDict.keys():
                                    changeInFinalSupportOfInitialTypeDict[CLADECOLLECTION.initialType.footPrint] = 0
                            for TYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
                                if TYPE.footPrint in changeInFinalSupportOfInitialTypeDict.keys():
                                    changeInFinalSupportOfInitialTypeDict[TYPE.footPrint] = changeInFinalSupportOfInitialTypeDict[TYPE.footPrint] + 1
                                else:
                                    changeInFinalSupportOfInitialTypeDict[TYPE.footPrint] = 1

    # Collect data to write out the type info and to make the strucuture. I.e. we will write out the types within the Maj ITS2 categories within the clades
    outPutDoc = []
    for CLADE in config.args.cladeList:
        collectionOfNonCoDomTypeMajs = []
        coDomTypes = []
        nonCoDomTypes = []
        numSamples = 0

        for SAMPLE in config.abundanceList:
            for cladeCollection in SAMPLE.cladeCollectionList:
                if cladeCollection.clade == CLADE: # Then this is a sample with a cladeCollection within the given clade
                    numSamples += 1
                    if cladeCollection.initialType.coDom: # Then this is a cladeCollection of the clade in question that is a supported type and is not a coDom
                        coDomTypes.append(cladeCollection.initialType)
                    elif not cladeCollection.initialType.coDom:
                        if cladeCollection.initialType.supportedType:
                            nonCoDomTypes.append(cladeCollection.initialType.name)
                            collectionOfNonCoDomTypeMajs.append(cladeCollection.initialType.maj)

        # ORDERED MAJs
        orderedCollectionOfNonCoDomSupportedTypeMajs = [a[0] for a in sorted({unique: collectionOfNonCoDomTypeMajs.count(unique) for unique in set(collectionOfNonCoDomTypeMajs)}.items(), key=lambda x: x[1], reverse=True)]
        # ORDERED CODOMs
        # This is now an ordered list of the coDomTypes
        orderedCoDomTypes = [a[0] for a in sorted({unique: coDomTypes.count(unique) for unique in set(coDomTypes)}.items(), key=lambda x: x[1], reverse=True)]
        if len(orderedCollectionOfNonCoDomSupportedTypeMajs) > 0 or len(orderedCoDomTypes) > 0:
            outPutDoc.append(('Clade ' + CLADE + ': Number of non-CoDom Maj = ' + str(len(set(collectionOfNonCoDomTypeMajs))) + ' / Number of CoDoms = ' + str(len(set(coDomTypes))) + ' / Number of Types = ' + str(len(set(coDomTypes)) + len(set(nonCoDomTypes))) + ' / Samples with supported Type = ' + str(len(coDomTypes) + len(nonCoDomTypes)) + '/' + str(numSamples), 'cladeHeader'))
            # INSERT FIGURE HERE
            # If the plots that we want to add exists
            fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/Clade' + CLADE + '_FstPlot.svg'
            paretoPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/Clade' + CLADE + '_pareto.svg'
            if os.path.isfile(fstPlotDir) and os.path.isfile(paretoPlotDir):
                outPutDoc.append(([fstPlotDir, paretoPlotDir], 'cladalPlot')) # <img src={{cell}} alt='Insufficient samples to warrant plot'>
            else:
                outPutDoc.append((config.args.saveLocation + '/html templates/image.jpg', 'plotError'))
            #outPutDoc.append(([None, 'Type', 'Initial:Final Support', 'Defining ITS2 Seqs'], 'subHeader'))
            for MAJ in orderedCollectionOfNonCoDomSupportedTypeMajs: # For each supported Maj that contains types for the given clade in order of abundance
                # Num of Types within maj
                listOfNonCoDomTypesWithinMaj = []
                for SAMPLE in config.abundanceList:
                    for cladeCollection in SAMPLE.cladeCollectionList:
                        if not cladeCollection.initialType.coDom and cladeCollection.initialType.maj == MAJ and cladeCollection.initialType.supportedType:
                            listOfNonCoDomTypesWithinMaj.append(cladeCollection.initialType)

                # Will now produce list of types which are in order of their abundance within the MAJ
                #orderedlistOfNonCoDomTypesWithinMaj = [a[0] for a in sorted({unique: listOfNonCoDomTypesWithinMaj.count(unique) for unique in set(listOfNonCoDomTypesWithinMaj)}.items())]
                orderedlistOfNonCoDomTypesWithinMaj = [a[0] for a in sorted({unique: listOfNonCoDomTypesWithinMaj.count(unique) for unique in set(listOfNonCoDomTypesWithinMaj)}.items(), key=lambda x: x[1], reverse=True)]

                if MAJ in config.oursToLaJDict.keys():
                    outPutDoc.append(('MajITS2Seq ' + config.oursToLaJDict[MAJ] + ': Number of Types = ' + str(len(orderedlistOfNonCoDomTypesWithinMaj)) + ' / Number of Samples = ' + str(len(listOfNonCoDomTypesWithinMaj)), 'majHeader'))
                else:
                    outPutDoc.append(('MajITS2Seq ' + MAJ[3:] + ': Number of Types = ' + str(len(orderedlistOfNonCoDomTypesWithinMaj)) + ' / Number of Samples = ' + str(len(listOfNonCoDomTypesWithinMaj)), 'majHeader'))
                #INSERT PLOT HERE
                if MAJ in config.oursToLaJDict.keys():
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + config.oursToLaJDict[MAJ] + '/' + config.oursToLaJDict[MAJ] + '_FstPlot.svg'
                else:
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + MAJ + '/Clade C ' + MAJ + '_FstPlot.svg'
                if os.path.isfile(fstPlotDir):
                    outPutDoc.append((fstPlotDir, 'majPlot'))
                else:
                    outPutDoc.append((config.args.saveLocation + '/html templates/image.jpg', 'plotError'))
                outPutDoc.append((['Type', 'Initial:Final Support', 'Defining ITS2 Seqs'], 'subHeader'))
                for TYPE in orderedlistOfNonCoDomTypesWithinMaj: # For each type within the given maj within the given clade
                    # try:
                    #     supportChangeValue = changeInFinalSupportOfInitialTypeDict[TYPE.footPrint]
                    #     absoluteFinalSupportOfInitialTypeDict[TYPE.footPrint] = listOfNonCoDomTypesWithinMaj.count(TYPE) + supportChangeValue
                    # except:
                    #     supportChangeValue = 0
                    #     absoluteFinalSupportOfInitialTypeDict[TYPE.footPrint] = listOfNonCoDomTypesWithinMaj.count(TYPE)
                    outPutDoc.append((['Type ' + TYPE.name, str(listOfNonCoDomTypesWithinMaj.count(TYPE)) + ':' + str(listOfNonCoDomTypesWithinMaj.count(TYPE) + changeInFinalSupportOfInitialTypeDict[TYPE.footPrint]),  TYPE.name.split('-')[0]], 'typeInfo'))
                    for INTRA in TYPE.name.split('-')[1:]:#Add the remaining intras on a new line
                        outPutDoc.append(([ None, None, INTRA], 'typeInfo'))

            #Now add coDom types
            if len(orderedCoDomTypes) > 0:
                outPutDoc.append(('Co-Dominant Types = ' + str(len(orderedCoDomTypes)), 'majHeader'))
            # outPutDoc.append('\tType' + '\tInitial Support' + '\tFinalSupport' + '\tDefining ITS2 Seqs')
            # Need an ordered type

            for CODOMTYPE in orderedCoDomTypes: # For each coDom type in order of abundance
                # try:
                #     supportChangeValue = changeInFinalSupportOfInitialTypeDict[CODOMTYPE.footPrint]
                #     absoluteFinalSupportOfInitialTypeDict[CODOMTYPE.footPrint] = coDomTypes.count(CODOMTYPE) + supportChangeValue
                # except:
                #     supportChangeValue = 0
                #     absoluteFinalSupportOfInitialTypeDict[CODOMTYPE.footPrint] = coDomTypes.count(CODOMTYPE)
                outPutDoc.append(([CODOMTYPE.name, str(coDomTypes.count(CODOMTYPE)) + ':' + str(coDomTypes.count(CODOMTYPE) + changeInFinalSupportOfInitialTypeDict[CODOMTYPE.footPrint]), re.split(r'/|-', CODOMTYPE.name)[0]], 'typeInfo'))
                for INTRA in re.split(r'/|-', CODOMTYPE.name)[1:]:
                    outPutDoc.append(([ None, None, INTRA], 'typeInfo'))
            outPutDoc.extend([(None, 'blankRow'), (None, 'blankRow'), (None, 'blankRow'), (None, 'blankRow')])

    jinjaEnvironment = Environment(loader=FileSystemLoader(config.args.saveLocation + r'\html templates'), trim_blocks=True)
    jinjaTemplate = jinjaEnvironment.get_template('typeCharacterisation__TEMPLATE.html')
    stringThing = jinjaTemplate.render(outPut=outPutDoc)
    htmlString = [a for a in stringThing.split('\n')]
    if not os.path.exists(config.args.saveLocation + '/html outputs'):
        os.makedirs(config.args.saveLocation + '/html outputs')
    #writeListToDestination(config.args.saveLocation +'/html outputs/typeCharacterisation__OUTPUT.html', htmlString)

    return changeInFinalSupportOfInitialTypeDict, initialSupportOfInitialTypesDict, htmlString


def writeMajMatricesToFile(coldists):
    print('Running writeMajMatricesToFile()')
    argStringCSVList = []
    #graphicHtml = []
    #graphicHtml.extend(['<table>', '<tr style=\'font-weight:bold; font-size:150%\'><td>Graphical Outputs</td></tr>'])
    for CLADE in config.args.cladeList:
        # Firstly write out a cladal Matrice before looking to do the Maj-based matrices
        listOfSamples = []
        for SAMPLE in config.abundanceList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    listOfSamples.append(CLADECOLLECTION.foundWithinSample)
        if len(listOfSamples) > 2: # There will be some clades that don't have any clade collections in them. In this case we can't make a graphic as there are no samples.
            #graphicHtml.extend(createMajMatrixFromColDists(listofsamples=listOfSamples, clade=CLADE, Fstcoldistdict=coldists, ismaj=False))
            argStringCSVList.append(createMajMatrixFromColDists(listofsamples=listOfSamples, clade=CLADE, Fstcoldistdict=coldists, ismaj=False))
            # Then write out the Maj matrices if appropriate
            # Identify the Majs
            MAJListFromSupportedNonCoDomTypes = []
            for SAMPLE in config.abundanceList:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE and not CLADECOLLECTION.initialType.coDom and CLADECOLLECTION.initialType.supportedType: # for each Maj ITS2 that isn't ONLY found in coDoms and has more than the threshold for supported types
                        MAJListFromSupportedNonCoDomTypes.append(CLADECOLLECTION.initialType.maj)
            # Create a dictionary key = Maj value = list of samples that contain a codom type that has the given Maj as one of the Majs
            # to add samples that have a coDom initial type that has the Maj in question to the sample list when doing the Majonly plots (i.e. not cladal plots)
            majToCoDomListOfSamplesDict = {}
            for SAMPLE in config.abundanceList:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE and CLADECOLLECTION.initialType.coDom:
                        for KEY_MAJ in CLADECOLLECTION.initialType.coDomDict:
                            if CLADECOLLECTION.initialType.coDomDict[KEY_MAJ]:
                                if KEY_MAJ in majToCoDomListOfSamplesDict.keys():#If the Maj is already in the dict
                                    if SAMPLE.name not in majToCoDomListOfSamplesDict[KEY_MAJ]:# If the sample not already in the sample list for that Maj
                                        tempList = majToCoDomListOfSamplesDict[KEY_MAJ]
                                        tempList.append(SAMPLE.name)
                                        majToCoDomListOfSamplesDict[KEY_MAJ] = tempList
                                else:
                                    majToCoDomListOfSamplesDict[KEY_MAJ] = [SAMPLE.name]


            for MAJ in set(MAJListFromSupportedNonCoDomTypes):# for each Maj ITS2 that isn't ONLY found in coDoms and has more than the threshold for supported types
                #supportedTypesWithinMajList = []
                listOfSamples = [] # List of samples that have the MajType in question and are not samples with a coDom type # plotting all types so that degradation with lower cutoff level can be seen
                for SAMPLE in config.abundanceList:
                    for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                        # Here we make two lists, supportedTypesWithinMajlist and listOfSamples. List of samples is what we send to make the graph
                        # supportedTypesWithinMajList is what we use to see if we have more than the maj plot threshold cuttoff i.e. if it's worth us making a plot
                        # We may need to get rid of the majPlotThreshold If we are going to put a graph in with every Maj text part
                        # if CLADECOLLECTION.initialType.maj == MAJ and not CLADECOLLECTION.initialType.coDom and CLADECOLLECTION.initialType.supportedType:
                        #     supportedTypesWithinMajList.append(CLADECOLLECTION.initialType.name)
                        if CLADECOLLECTION.initialType.maj == MAJ and not CLADECOLLECTION.initialType.coDom:
                            listOfSamples.append(CLADECOLLECTION.foundWithinSample)
                # Convert the raw maj to laJ if possible so that we have this name on the plot.
                # Also majToCoDomListOfSamplesDict's values are already converted to must use converted if in laJDict
                #Need to have a conditional here to check that the Maj in question is found in any of the coDom types of that clade
                # if it is not then we won't need to add any additional samples to the sample list
                if MAJ in config.oursToLaJDict.keys():
                    convertedMaj = config.oursToLaJDict[MAJ]
                    isConverted = True
                    if convertedMaj in majToCoDomListOfSamplesDict.keys():
                        listOfSamples.extend(majToCoDomListOfSamplesDict[convertedMaj])

                else:
                    if MAJ in majToCoDomListOfSamplesDict.keys():
                        listOfSamples.extend(majToCoDomListOfSamplesDict[MAJ])
                    isConverted = False
                    convertedMaj = MAJ
                if len(listOfSamples) >= 4: # Only make a MajPlot if there are more than or equal to the maj plot threshold cutoff number of supported types (not coDom types) within the maj
                    #graphicHtml.extend(createMajMatrixFromColDists(Fstcoldistdict=coldists, listofsamples=listOfSamples, clade=CLADE, ismaj=True, convertedmaj=convertedMaj, maj=MAJ, isconverted=isConverted))
                    argStringCSVList.append(createMajMatrixFromColDists(Fstcoldistdict=coldists, listofsamples=listOfSamples, clade=CLADE, ismaj=True, convertedmaj=convertedMaj, maj=MAJ, isconverted=isConverted))
    # Here write out the argstringcsvlist to file
    writeListToDestination(config.args.saveLocation + r'\matrices outputs\RArgList.txt', argStringCSVList)
    producePlot()
    # Then pass the directory to it to the plot function
    #graphicHtml.append('</table>')
    return #graphicHtml

def createMajMatrixFromColDists(Fstcoldistdict, listofsamples, clade, ismaj,  maj = None, convertedmaj = None, isconverted = None):

    listOfSamples = listofsamples
    # Create empty matrix that is of dimensions [rows(n samples + 1)][columns(number of samples + 1)]
    FstMatrix = [list(listOfSamples)]
    FstMatrix[0].insert(0, 'Samples')
    for sample in listOfSamples:
        t = ['N/A']*len(listOfSamples)
        t.insert(0, sample)
        FstMatrix.append(t)

    # Matrix is FstMatrix[row][column]
    # Fill matrix from the column Fsts in FstColDist
    row = 1
    while row < len(FstMatrix): # For each row
        col = 1
        while col < len(FstMatrix[row]): # For each col
            sampleOne = FstMatrix[0][col]
            sampleTwo = FstMatrix[row][0]
            if sampleOne == sampleTwo:
                FstMatrix[row][col] = 0.00
            else:
                    FstMatrix[row][col] = FstMatrix[col][row] = Fstcoldistdict[config.args.cladeList.index(clade)][frozenset({sampleOne, sampleTwo})]
            col += 1
        row += 1
    # Get rid of the column and row headers to leave just pure numerical matrix
    del FstMatrix[0]
    FstMatrix = [a[1:] for a in FstMatrix]

    if ismaj: # If it is a sub cladal Maj level Fst
        if not os.path.exists(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/IncludingMaj/' + convertedmaj):
            os.makedirs(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/IncludingMaj/' + convertedmaj)
        write2DListToDestination(config.args.saveLocation + '\\matrices outputs\\Clade' + clade + '\\' + 'IncludingMaj' + '\\' + convertedmaj + '\\Matrix.dist', FstMatrix)
        writeListToDestination(config.args.saveLocation + '\\matrices outputs\\Clade' + clade + '\\' + 'IncludingMaj' + '\\' + convertedmaj + '\\Header.head', listOfSamples)
        # Create the info list
        infoList = []
        for SAMPLE in config.abundanceList:
            if SAMPLE.name in listOfSamples:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                        if CLADECOLLECTION.initialType.coDom: # If coDom then look up maj for the given sample
                            if CLADECOLLECTION.initialType.maj[SAMPLE] in config.oursToLaJDict:
                                majITS2 = config.oursToLaJDict[CLADECOLLECTION.initialType.maj[SAMPLE]]
                            else:
                                majITS2 = CLADECOLLECTION.initialType.maj[SAMPLE]
                            infoList.append(','.join([SAMPLE.name, SAMPLE.hostTaxon, SAMPLE.reef, SAMPLE.region, majITS2, CLADECOLLECTION.initialType.name]))
                        else: # If not coDom use the Maj that is assigned to the initial type
                            if CLADECOLLECTION.initialType.maj in config.oursToLaJDict:
                                majITS2 = config.oursToLaJDict[CLADECOLLECTION.initialType.maj]
                            else:
                                majITS2 = CLADECOLLECTION.initialType.maj
                            infoList.append(','.join([SAMPLE.name, SAMPLE.hostTaxon, SAMPLE.reef, SAMPLE.region, majITS2, CLADECOLLECTION.initialType.name]))
        writeListToDestination(config.args.saveLocation + '\\matrices outputs\\Clade' + clade + '\\' + 'IncludingMaj' + '\\' + convertedmaj + '\\SampleInfo.info', infoList)
        if isconverted:
            #htmlSnippet = producePlot(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + convertedmaj, convertedmaj, str(ismaj).upper())
            #producePlot(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + convertedmaj, convertedmaj, str(ismaj).upper())
            argString = ','.join([(config.args.saveLocation).replace('\\','/'), (config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + convertedmaj).replace('\\','/'), convertedmaj, str(ismaj).upper()])
        else:
            #htmlSnippet = producePlot(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + convertedmaj, 'Clade ' + clade + ' ' + maj, ismaj)
            #producePlot(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + convertedmaj, 'Clade ' + clade + ' ' + maj, ismaj)
            argString =','.join([(config.args.saveLocation).replace('\\','/'), (config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + convertedmaj).replace('\\','/'), 'Clade ' + clade + ' ' + maj, str(ismaj).upper()])
    else: # If a cladal level Fst
        # Check to see if the save directory already exists, if not, create it.
        if not os.path.exists(config.args.saveLocation + '/matrices outputs/Clade' + clade):
            os.makedirs(config.args.saveLocation + '/matrices outputs/Clade' + clade)

        write2DListToDestination(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/Matrix.dist', FstMatrix)
        writeListToDestination(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/Header.head', listOfSamples)
        infoList = []
        for SAMPLE in config.abundanceList:
            if SAMPLE.name in listOfSamples:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                        # Check to see if the Maj name can be LaJconverted
                        # If so then convert it and add this to the info list in place of the raw Maj name
                        if CLADECOLLECTION.initialType.coDom: # If coDom then look up maj for the given sample
                            if CLADECOLLECTION.initialType.maj[SAMPLE] in config.oursToLaJDict:
                                majITS2 = config.oursToLaJDict[CLADECOLLECTION.initialType.maj[SAMPLE]]
                            else:
                                majITS2 = CLADECOLLECTION.initialType.maj[SAMPLE]
                            infoList.append(','.join([SAMPLE.name, SAMPLE.hostTaxon, SAMPLE.reef, SAMPLE.region, majITS2, CLADECOLLECTION.initialType.name]))
                        else: # If not coDom use the Maj that is assigned to the initial type
                            if CLADECOLLECTION.initialType.maj in config.oursToLaJDict:
                                majITS2 = config.oursToLaJDict[CLADECOLLECTION.initialType.maj]
                            else:
                                majITS2 = CLADECOLLECTION.initialType.maj
                            infoList.append(','.join([SAMPLE.name, SAMPLE.hostTaxon, SAMPLE.reef, SAMPLE.region, majITS2, CLADECOLLECTION.initialType.name]))
        writeListToDestination(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/SampleInfo.info', infoList)
        # htmlSnippet = producePlot(config.args.saveLocation + '/matrices outputs/Clade' + clade , 'Clade' + clade, ismaj)
        argString =','.join([(config.args.saveLocation).replace('\\','/'), (config.args.saveLocation + '/matrices outputs/Clade' + clade ).replace('\\','/'), 'Clade' + clade, str(ismaj).upper()])
        #producePlot(config.args.saveLocation + '/matrices outputs/Clade' + clade , 'Clade' + clade, ismaj)

    print('Completed createMajMatrixFromColDists()')
    return argString


def producePlot():
#def producePlot(filelocation, svgName, ismaj):
    # Create the SVG #
    # plotDir = config.args.saveLocation + r'\plot outputs'
    # try:
    #     os.makedirs(plotDir)
    # except FileExistsError:
    #     pass
    batFile = config.args.rootLocation + r"\.bat scripts\rPlotOne.R"
    argDir = (config.args.saveLocation + r'\matrices outputs\RArgList.txt').replace('\\','/')
    # cmd = [config.args.rscriptLocation, batFile, config.args.saveLocation, filelocation.replace('\\','/'), svgName, str(ismaj).upper()]
    cmd = [config.args.rscriptLocation, batFile, config.args.saveLocation, argDir]
    print('Running ' + str(cmd))
    call(cmd)

    # Create the SVG html snippet
    # htmlSnippet = []
    # if not ismaj: # If cladal then need to add both the fstplot and the pareto plot to the html snippet
    #     # This adds two rows of two cells
    #     # First row is the headings for the images
    #     # Second row is the images themselves, side by side
    #     htmlSnippet.extend(['<tr style=\'font-weight:bold; font-size:120%\'><td>' + svgName + '</td></tr>',
    #                         '<tr><td><img src=\'' + filelocation + '\\' + svgName  + '_FstPlot.svg\' alt=\'Error in displaying plot\'></td><td><img src=\'' + filelocation + '\\' + svgName + '_pareto.svg\' alt=\'Error in displaying plot\'></td></tr>'])
    # else: # If this is a Maj plot then only need the fstplot in the html
    #     htmlSnippet.extend(['<tr style=\'font-weight:bold; font-size:120%\'><td>' + svgName + '_FstPlot</td></tr>', '<tr><td><img src=\'' + filelocation + '\\' + svgName  + '_FstPlot.svg\' alt=\'Error in displaying plot\'></td></tr>'])

    return #htmlSnippet

def createAbundanceListMultiProcess(miniAbundanceList):
    col = 1
    outPut = []
    while col < len(miniAbundanceList[0][:-1]):# For each sample i.e. col not including the final clade column
        row = 6
        tempList = []
        while row < len(miniAbundanceList):# For each ITS2 seq
            if int(miniAbundanceList[row][col]) != 0:#Then this sample contains this seq
                tempList.append(miniAbundanceList[row][0] + miniAbundanceList[row][-1] + '\t' + str(miniAbundanceList[row][col]))
            row += 1

        sortedList = sorted(tempList, key=lambda x: int(x.split('\t')[1]), reverse=True)
        sortedList.insert(0, miniAbundanceList[0][col])
        outPut.append(sortedList)
        col += 1
        print(multiprocessing.current_process().name + ' ' + str(col))
    return outPut

def createLogTransedFstColDists(masterSeqDistancesDict):
    print('Running createLogTransedFstColDists()')
    FstColDist = []
    for CLADE in config.args.cladeList: # For each clade
        tempColDist = []
        listOfCladeCollections = []
        for SAMPLE in config.abundanceList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    listOfCladeCollections.append(CLADECOLLECTION)
        for a,b in itertools.combinations(listOfCladeCollections, 2): # For each paired combination of sample
            tempColDist.append([{a.foundWithinSample, b.foundWithinSample}, str(CalculateSinglePairwiseFst(a.listOfits2SequenceOccurances, b.listOfits2SequenceOccurances, masterSeqDistancesDict))])
        # Now convert the list structure to a dict structure
        tempDict = {frozenset(item[0]): item[1] for item in tempColDist}

        print('Clade ' + CLADE + ' Fst colDist complete')
        FstColDist.append(tempDict)
    print('createLogTransedFstColDists() complete')
    # write2DListToFile('columnFormattedHumeFsts', '.coldist', FstColDist)
    return FstColDist

def assignInitialTypes():

    # Once we have worked out the initial types we will add them to each of the samples that can be identified. If they can't be identified at this point then we will add unsupported as the type
    for CLADE in config.args.cladeList:
        listOfFootprints  = []
        for SAMPLE in config.abundanceList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.footPrint not in listOfFootprints:
                        listOfFootprints.append(CLADECOLLECTION.footPrint)
        for FOOTPRINT in listOfFootprints:
            coDom = False
            supportedType = False # At the given level of support

            listOfSamplesThatContainFootprint = [] # List of the sample names of samples that have a cladeCollection that match the footprint (exactly)
            listOfMajsForSamplesWithFootPrint = []
            for SAMPLE in config.abundanceList:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE and CLADECOLLECTION.footPrint == FOOTPRINT:
                        listOfSamplesThatContainFootprint.append(SAMPLE.name)
                        listOfMajsForSamplesWithFootPrint.append(CLADECOLLECTION.maj)

            if len(listOfSamplesThatContainFootprint) >= config.args.typeSupport: # There are enough samples with this footprint to define it as a type, continue checking it for being a coDom
                supportedType = True
                if len(listOfMajsForSamplesWithFootPrint) >= 2*config.args.coDomSupport and len(set(listOfMajsForSamplesWithFootPrint)) > 1:
                # This is a coDom strictly speaking and there is enough support for it to be considered a coDom according to our extra conservaive coDomSupportvalues
                # The *2 is there because for a codom we must have at least two different MajITS2s and at least the config.args.coDomSupport number of samples containing each of them.
                    coDomSupportedMajDict = {Maj: True for Maj in set(listOfMajsForSamplesWithFootPrint)} # Create a dictionary of all of the MajITS2s found in this type/footprint and evaluate whether they are found in more than 2 samples. If the count of Trues is less than 2 at the end then this is not a codom [under our extra conservative ruling]. If more than one are above the 2 then name accordingly.
                    for Maj in coDomSupportedMajDict.keys(): # Check each one to see if found in more than 2 samples
                        if listOfMajsForSamplesWithFootPrint.count(Maj) < 2: # If we don't find the support then don't call it a coDom
                            coDomSupportedMajDict[Maj] = False
                    #Now do the count of Falses vs. Positive, if Trues greater than 1 then this is a CoDom and we need to work out the name
                    if list(coDomSupportedMajDict.values()).count(True) > 1:
                        # This is a supported CoDom now name it:
                        coDom = True
                        newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL', coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE, coDomDict=coDomSupportedMajDict)
                        #newSymbiodiniumType = symbiodiniumType(coDom = coDom, supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE, abundanceList=abundanceList, oursToLaJDict=oursToLaJDict, coDomDict=coDomSupportedMajDict)
                        # Put the new Type into the samples' collectionlits that have it
                        addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
                    else:
                    # This cannot be called a coDom and we can exit out at this point
                    # This is an interesting case. So here we have a type that is not considered a coDom using our conservative criteria
                    # But this type does have multiple MajITS2s (just not enough support for it to be considred a coDom
                    # This is going to cause trouble later on when we come to write out typeCharacterisation
                    # As we go through the MajITS2s and this type will therefore come up in two MajITS2 categories.
                    # We need to find away to make sure that it is only related to
                    # We are going to solve this problem by adding self.maj to the symbiodiniumType class when it is an initial type
                    # Then when we come to identify the Majs we can use the inital type majs rather than the clade collection Majs

                        newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL',coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)#, abundanceList=abundanceList, oursToLaJDict=oursToLaJDict)
                        addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
                else: # # This is cannot be called a coDom and we can exit out at this point
                    newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL', coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)#, abundanceList=abundanceList, oursToLaJDict=oursToLaJDict)
                    addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
            else: # This is an unsupportedType
                newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL',coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)#, abundanceList=abundanceList, oursToLaJDict=oursToLaJDict)
                addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
    return

def addTypeToSamples(newSymType, listOfSamplesThatContainFootprint):
    for SAMPLE in config.abundanceList:
         if SAMPLE.name in listOfSamplesThatContainFootprint:
             for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                 if CLADECOLLECTION.clade == newSymType.clade and CLADECOLLECTION.footPrint == newSymType.footPrint: # Then this is a cladeCollection that we want to add the new Symbiodinium Type to
                    CLADECOLLECTION.addInitialType(newSymType)

def createMasterSeqDistances():
    #sequsedinFst was calculated alongside the cladeCollectionList and so only contains sequences appearing at the given cutoff percentage, i.e. those that appear in the cladeCollectionList and so will be used in the creation of matrices and identification of types
    #Then for each of the lists, create fastas for each clade
    #Then read each fasta into mothur and output the dist file
    #Read these files in and concatenate them so that there is a mast dist file
    #Turn this into a dictionary

    print('Running createMasterSeqDistances()')
    outPutFastaCollection = []
    for CLADE in config.args.cladeList:
        sequencesUsedInFst = []
        listOfNames = []
        for sample in config.abundanceList:
            for cladeCollection in sample.cladeCollectionList: # For each cladeCollection
                if cladeCollection.clade == CLADE:
                    for sequenceOccurance in cladeCollection.listOfits2SequenceOccurances:
                        if sequenceOccurance.name not in listOfNames:
                            sequencesUsedInFst.append(sequenceOccurance)
                            listOfNames.append(sequenceOccurance.name)
        outPutFasta = []
        for seqoccurance in set(sequencesUsedInFst):
            outPutFasta.extend(['>' + seqoccurance.name, seqoccurance.sequence])
        outPutFastaCollection.append(outPutFasta)

    distsList = []
    batchFile = []
    for fasta in outPutFastaCollection: # Write the fastas to the mothur directory so that they can be processed
        if len(fasta) > 4: # Check to see if there are at least two sequences in the fasta so that a pairwise comparison can be done by mothur
            clade = config.args.cladeList[outPutFastaCollection.index(fasta)]
            dest = config.args.saveLocation + r'\seq distance files\clade' + clade + r'SequencesUsedInFst.fas'
            writeListToDestination(dest, fasta) #Write out the fasta
            batchFile.append("pairwise.seqs(fasta=clade" + clade + "SequencesUsedInFst.fas, processors=" + str(config.args.numProcessors) + ", outputdir=" + config.args.saveLocation + r'\seq distance files, inputdir=' + config.args.saveLocation+ r'\seq distance files' + ')')
    batDir = config.args.saveLocation + '/.bat scripts'
    writeListToDestination(batDir + '/mothurDistCalc.bat', batchFile)
    # Here we should have written the batch file
    # Now just to attempt to fun it on the command line
    cmd = [config.args.mothurLocation, 'mothurDistCalc.bat']
    print('Running ' + str(cmd) + ' in ' + batDir)
    call(cmd, shell=True, cwd=batDir)
    for clade in config.args.cladeList:
        if os.path.exists(config.args.saveLocation + r'\seq distance files\clade' + clade + 'SequencesUsedInFst.dist'):
            distsList.extend(readDefinedFileToList(config.args.saveLocation + r'\seq distance files\clade' + clade + 'SequencesUsedInFst.dist'))# Read the output file back in

    return distsList

def assignCladeCollections():

    for SAMPLE in config.abundanceList:
        for CLADE in config.args.cladeList: # Sequences in the different clades are too divergent to be compared so we have to work on a cladlly separated basis, we are only working with the three main clades, A, C and D
            totalSeqs = SAMPLE.totalSeqs
            cladeSpecificSeqs = sum([a.abundance for a in SAMPLE.compComplement.listOfits2SequenceOccurances if a.clade == CLADE]) # Number of sequence the given sample has that are of the given clade
            if float(cladeSpecificSeqs/totalSeqs) >= config.args.cladeCollectionCutoff:  # If this sample has a proportion of clade X creater than 10% then we will add a cladeCollection to the sample
                cutOffValue = (int(totalSeqs) * float(cladeSpecificSeqs/totalSeqs)) * config.args.cutOff # CutOff that is the number of sequences that represents X% of the cladal proportion
                tempListOfits2SequenceOccurances = [its2SequenceOccurance(name=a.name, abundance=a.abundance,clade=a.clade, sequence=a.sequence) for a in SAMPLE.compComplement.listOfits2SequenceOccurances if a.abundance >= cutOffValue and a.clade == CLADE] # Need to make sure to make a copy of the sequence occurances here so that when we put them into the cladecollection we don't change the abundnaces etc. in the
                tempTotalSeqs = sum([a.abundance for a in tempListOfits2SequenceOccurances]) # Total number of seqs in the its2sequenceoccurances that were above the cutoof for the given clade
                i = 0
                while i < len(tempListOfits2SequenceOccurances):
                    tempListOfits2SequenceOccurances[i].abundance = (tempListOfits2SequenceOccurances[i].abundance/tempTotalSeqs)*1000 # Normalise the abundances to 1000
                    i += 1
                SAMPLE.addCladeCollection(cladeCollection(CLADE, config.args.cutOff, tempListOfits2SequenceOccurances, SAMPLE.name, cladeSpecificSeqs/totalSeqs))

    return



def inferFinalSymbiodiniumTypes():

    for CLADE in config.args.cladeList:

        # Get a list of all of the intialSymbiodinium types that have unique footprints
        footprintList = []
        typeList = []
        typeCountDict = {}
        for SAMPLE in config.abundanceList:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.initialType.footPrint not in footprintList and CLADECOLLECTION.initialType.supportedType:
                            footprintList.append(CLADECOLLECTION.initialType.footPrint)
                            typeList.append(CLADECOLLECTION.initialType) # We want to keep the actual type class rather than just the frozenset footprint because we want to be able to work out abundances etc. later on

        # Search for initial type footprints within the complete complement of each sample for a given clade
        # Get list of types that are found in the complete complement first.
        # We need to be careful at this stage, the unsupported types will always fall into the original sample they came from, we need a way of checking whether these insignificant types become supported once we start to use the full complent of samples
        # We also need to be careful because if we simply look for the most complicated type in a sample and we have a sample that has 5 unique intras that are all above the 10% mark but only 1 of them is found in a supported type then this sample will be identified as something very simple which is the one intra that it does have
        # So if we have a sample that still has an intra present above the 10% level (of the cladal sequences) that hasn't been found in one of the types identified in it then we should call this sample unidentified.
        # In order to work out if each of the initial unsupported types becomes supported once the completeComplements are checked and equally whether any of the supported types become unsupported, we need to do 3 phases
        # At the end of this the only Type that may go unidentified by this method is one that is defined according to an intra that is almost always found below the 10% mark - I'm guessing these types will be pretty rare in biology and although some types,
        # i.e. ST will have identifying intras below the 10% mark in at least some samples they will be above the 10% mark so we will have picked them up in our initial type anlysis
        # Phase 1 - Go through all samples pushing in all types, supported and unsupported.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        # Phase 3 - Go through all samples and remove unsupported types then do iter comparisons of types in the sample keeping only the [type with the longest footprint of any two types that share any defining intras]*** actually we need to work through
        # each of these scenarios carefully. If the shorter footprint fits into the longer footprint then yes only keep longer footprint. But if two footprints share intras but one doesn't fit into the other then call it unresolved for the time being.
        # so step one is check to see if they have any intras in common. Then step two is ask if one is subset of the other.
        # Probably if they share intras but one doesn't fit in the other then we should keep both types. We would then have to spit the final proportion of the types between the two.
        # THe above is impossible to resolve. Consider the three type foot prints [1,2,3], [1,4,5] and [7,4,6] you can see the problem. firt two share the 1 second two share the 4 but can't make a type of all three as first and third do not share any. So lump them all into one list and make a super unresolved Symbiodinium type of them. Not ideal but can't see any other way at current.
        # Phase 4 - Finally if not all intras found at >10% (the defined cutoff) of the cladal sequences in a sample then consider this sample unidentified,
        # Once all types have been identified we can add the type information (type by type to the finalTypesList) (and, list of ITS2 occurances found only in the types)

        # Phase 1 - Go through all samples pushing in all types, supported and unsupported. With the caveat that the new types must have the old initial type's footprint in it.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        # If there are no further types then we don't have a finaltypecladecollectionlist
        for SAMPLE in config.abundanceList:
            finalTypesList = []
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE: # Then this sample has a set of intras from the given clade that are above the given cladeCollectionCuttoff
                    listOfIntrasInSample = [occurance.name for occurance in SAMPLE.compComplement.listOfits2SequenceOccurances if occurance.clade == CLADE]
                    for TYPE in typeList:
                        if TYPE.footPrint.issubset(listOfIntrasInSample) and CLADECOLLECTION.initialType.footPrint.issubset(TYPE.footPrint): # Check to see if the intras in the TYPE.footPrint are found in the listOfIntrasInSample list.
                            # Update type count dictionary and add this type to the sample's finalTypesList
                            if TYPE not in typeCountDict.keys():
                                typeCountDict[TYPE] = 1
                            else:
                                typeCountDict[TYPE] = typeCountDict[TYPE] + 1
                            finalTypesList.append(TYPE)
                    if len(finalTypesList) > 0:
                        SAMPLE.finalTypeCladeCollectionList.append(finalTypeCladeCollection(foundWithinSample=SAMPLE.name, clade=CLADE, cutoff=config.args.cutOff, listOfFinalTypes=finalTypesList))
        # Phase one and two are complete here. Time for phase three


        # Phase three: first we need to go through and remove all unsupported types from the type lists

        for SAMPLE in config.abundanceList:
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == CLADE:
            # Now iter comparisons and get rid of those types that are subsets of others
                    listOfTypesToGetRidOf = []
                    for a,b in itertools.combinations(FINALTYPECLADECOLLECTION.listOfFinalTypes, 2):
                        if len([intra for intra in a.footPrint if intra in b.footPrint]) > 0: # Then the two footprints share at least one intra i.e. one may include the other
                            if a.footPrint.issubset(b.footPrint): # Then we need to get rid of a
                                if a not in listOfTypesToGetRidOf:
                                    listOfTypesToGetRidOf.append(a)
                            elif b.footPrint.issubset(a.footPrint): # Then we need to get rid of b
                                if b not in listOfTypesToGetRidOf:
                                    listOfTypesToGetRidOf.append(b)

                    for TYPE in listOfTypesToGetRidOf:
                        FINALTYPECLADECOLLECTION.listOfFinalTypes.remove(TYPE)

            # Once all of the unsupported TYPES have been removed check to see if all of the samples > cutoff (for a given clade) intras are accounted for within the list of types identified
            # Also covert each of the SymbiodiniumTypes in the FINALCLADECOLLECTION.listoffinaltypes to new types that contain the information of the abundance of
                    for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                        if CLADECOLLECTION.clade == CLADE: # Then this is the clade in questions cladeCollection and we can use the its2 occurances in this as the list of intras found at above the cutOff
                            # Then we know that all of the major intras for this clade of this sample have been accounted for by the supported final types and we can go on to write out final types list and proportions
                            if set([A.name for A in CLADECOLLECTION.listOfits2SequenceOccurances]).issubset(FINALTYPECLADECOLLECTION.typeBasedCompCollection()):
                                FINALTYPECLADECOLLECTION.identified = True
                            else:
                                FINALTYPECLADECOLLECTION.identified = False
                            # We get rid of the inital type from the finaltypecladecollection.listoffinaltypes so that it doesn't get written out later
                            i = 0
                            while i < len(FINALTYPECLADECOLLECTION.listOfFinalTypes): # Get rid of the inital type from the final identification finaltype list
                                if FINALTYPECLADECOLLECTION.listOfFinalTypes[i].name == CLADECOLLECTION.initialType.name:
                                    del FINALTYPECLADECOLLECTION.listOfFinalTypes[i]
                                    i += -1
                                else:
                                    FINALTYPECLADECOLLECTION.listOfFinalTypes[i] = symbiodiniumType(clade=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].clade, footPrint=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].footPrint, maj=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].maj, coDomDict=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].coDomDict, coDom=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].coDom, typeOfType='FINAL', totalSeqs=SAMPLE.totalSeqs, typeSupport=typeCountDict[FINALTYPECLADECOLLECTION.listOfFinalTypes[i]], listofoccurences=SAMPLE.compComplement.listOfits2SequenceOccurances)
                                i += 1
                            break # Break out of the last for loop that checks to see if finaltypecladecollection is identified
                    break # Break out of the very first for loop that identifies a finaltypecladecollection of the given clade as there is only one per clade

    return

def typePrintOutString(coDom, listOfITS2Occurence, totalSeqs,  footprint, typeSupport, clade, typeoftype, coDomDict = None ): # Outputname
    footPrint = list(footprint)
    #make intraCladalProportionDict
    intraCladalProportionDict = {}
    cladalProportion = sum([occurence.abundance for occurence in listOfITS2Occurence if occurence.clade == clade])/totalSeqs
    for OCCURENCE in listOfITS2Occurence:
        if OCCURENCE.name in footPrint:
            intraCladalProportionDict[OCCURENCE.name] = OCCURENCE.abundance/(cladalProportion*totalSeqs)
    #convert the footPrint Intras, the intracladalpropdict and the coDOmdict to LaJeunesse language
    i = 0
    while i < len(footPrint):
        if footPrint[i] in config.oursToLaJDict.keys():
            if coDomDict:
                if footPrint[i] in coDomDict.keys():
                    coDomDict[config.oursToLaJDict[footPrint[i]]] = coDomDict[footPrint[i]]
                    del coDomDict[footPrint[i]]
            if footPrint[i] in intraCladalProportionDict.keys():
                intraCladalProportionDict[config.oursToLaJDict[footPrint[i]]] = intraCladalProportionDict[footPrint[i]]
                del intraCladalProportionDict[footPrint[i]]
            footPrint[i] = config.oursToLaJDict[footPrint[i]]
        i += 1

    sortedList = [a[0] for a in sorted(intraCladalProportionDict.items(), key=lambda x: x[1], reverse=True)]
    added = []
    if coDom:
        sortedCoDomKeys = [intra for intra in sortedList if intra in coDomDict]
        namePart1 = '/'.join([codomintra + '[' + str(format(intraCladalProportionDict[codomintra], '.2f')) + ']' for codomintra in sortedCoDomKeys if coDomDict[codomintra]]) # Add any coDom intras first
        added.extend([codomintra for codomintra in coDomDict.keys() if coDomDict[codomintra]])
    namePart2 = '-'.join([noncoDomIntras + '[' + str(format(intraCladalProportionDict[noncoDomIntras], '.2f')) + ']' for noncoDomIntras in sortedList if noncoDomIntras not in added]) # If it isn't already in the name because it is a codom then add in order of abundance within the sample
    if coDom and len(namePart2) > 0:
        outputName = '-'.join([namePart1, namePart2])
    elif coDom and len(namePart2) == 0:
        outputName = namePart1
    else:
        outputName = namePart2
    typeTotalProportion = sum([intraCladalProportionDict[intrasinfootprint] for intrasinfootprint in footPrint])# The summed proportions of all of the Type's intras
    outputName = outputName + '\t[' + str(typeSupport[0]) + ':' + str(typeSupport[1]) + ']'

    if typeoftype == 'INITIAL':
        return outputName
    else:
        return outputName#, typeTotalProportion, typeSupport

def writeSampleCharacterisationOutput(changeInFinalSupportOfInitialTypeDict, initialSupportOfInitialTypesDict):

    outPut = []
    outPut.extend([(['Sample details'], 'title'),([ None, 'Initial Type', '[Initial : Final support]', 'Further types', '[Initial:Final support]'], 'notBold'),(None, 'blankLine')])
    for SAMPLE in config.abundanceList:
        cladalPropCounter = {}
        for OCCURENCE in SAMPLE.compComplement.listOfits2SequenceOccurances:
            if OCCURENCE.clade in cladalPropCounter.keys():
                cladalPropCounter[OCCURENCE.clade] = cladalPropCounter[OCCURENCE.clade] + OCCURENCE.abundance
            else:
                cladalPropCounter[OCCURENCE.clade] = OCCURENCE.abundance
        for keys in cladalPropCounter.keys():
            cladalPropCounter[keys] = cladalPropCounter[keys]/SAMPLE.totalSeqs
        # Here we have the cladl Prop Counter which is a dict of key = clade, value = proportion of total seqs

        sortedCladalAbundanceTuple = [(a[0],a[1]) for a in sorted(cladalPropCounter.items(), key=lambda x:x[1], reverse=True)]
        cladalProportionString = ':'.join([a[0] + ' ' + str(format(a[1], '.2f')) for a in sortedCladalAbundanceTuple])
        outPut.append((' / '.join([SAMPLE.name, cladalProportionString, SAMPLE.hostTaxon, SAMPLE.region, SAMPLE.reef]), 'sample'))
        # Here we have written the header line of a sample. Now time to printout the intial and final types, in order of cladal abundance in the sample
        # For the final types we have a choice of two class parameters to sort them by, typesupport or type total proportion
        # Typesupport is the number of samples that that final type was found in (i.e. a type being a given footprint)
        # Total proportion being the number of sequences of the sample that make up that type, i.e. sum of abundances of each of the occurances within that type
        # We will use total proportion
        for CLADE in [a[0] for a in sortedCladalAbundanceTuple]:# Only clades that are represented in the sample in order of abundance
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if len(SAMPLE.finalTypeCladeCollectionList) > 0:
                        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                            if FINALTYPECLADECOLLECTION.clade == CLADE:
                                orderedListOfFinalTypes = sorted(FINALTYPECLADECOLLECTION.listOfFinalTypes, key=lambda x: x.typeTotalProportion, reverse=True)
                                initialTypeOutName = typePrintOutString(coDom=CLADECOLLECTION.initialType.coDom, listOfITS2Occurence=SAMPLE.compComplement.listOfits2SequenceOccurances, totalSeqs=SAMPLE.totalSeqs,  footprint=CLADECOLLECTION.initialType.footPrint, typeSupport=(initialSupportOfInitialTypesDict[CLADECOLLECTION.initialType.footPrint], initialSupportOfInitialTypesDict[CLADECOLLECTION.initialType.footPrint] + changeInFinalSupportOfInitialTypeDict[CLADECOLLECTION.initialType.footPrint]) , clade=CLADECOLLECTION.initialType.clade, typeoftype=CLADECOLLECTION.initialType.typeOfType, coDomDict = CLADECOLLECTION.initialType.coDomDict)
                                if len(orderedListOfFinalTypes) > 0:
                                    listOfFinalTypeOutputNames = [typePrintOutString(coDom=finalType.coDom, listOfITS2Occurence=SAMPLE.compComplement.listOfits2SequenceOccurances, totalSeqs=SAMPLE.totalSeqs,  footprint=finalType.footPrint, typeSupport=(initialSupportOfInitialTypesDict[finalType.footPrint], initialSupportOfInitialTypesDict[finalType.footPrint] + changeInFinalSupportOfInitialTypeDict[finalType.footPrint]), clade=finalType.clade, typeoftype=finalType.typeOfType, coDomDict = finalType.coDomDict) for finalType in orderedListOfFinalTypes]
                                    outPut.append(([None, initialTypeOutName.split('\t')[0], initialTypeOutName.split('\t')[1], listOfFinalTypeOutputNames[0].split('\t')[0],listOfFinalTypeOutputNames[0].split('\t')[1]], 'notBold'))
                                    for finalTypeOutPutName in listOfFinalTypeOutputNames[1:]:
                                        outPut.append(([None, None, None, finalTypeOutPutName.split('\t')[0], finalTypeOutPutName.split('\t')[1]], 'notBold'))
                                else:# If there aren't any further types due to all seqs in the finaltypecladecollection being removed
                                    outPut.append(([None, initialTypeOutName.split('\t')[0], initialTypeOutName.split('\t')[1], None, None], 'notBold'))
                    else:# If there aren't any further types due to there being no finaltypecladecollections
                        outPut.append(([None, initialTypeOutName.split('\t')[0], initialTypeOutName.split('\t')[1], None, None], 'notBold'))
        outPut.append((None, 'blankLine'))
                        # For each type in the final type list and for the initial type
                        # Do the function for the write out



    jinjaEnvironment = Environment(loader=FileSystemLoader(config.args.saveLocation + r'\html templates'), trim_blocks=True)
    jinjaTemplate = jinjaEnvironment.get_template('sampleCharacterisation__TEMPLATE.html')
    stringThing = jinjaTemplate.render(outPut=outPut)
    htmlString = [a for a in stringThing.split('\n')]

    return htmlString

def createHtmlHolder():
    tempList = []
    tempList.extend(['<!DOCTYPE html>', '<html>', '<head>', '<title>SymbiodiniumType Output</title>', '</head>', '<body>'])
    return tempList

def closeAndWriteHtmlHolder(htmlfile):
    htmlfile.extend(['</body>', '</html>'])
    writeListToDestination(config.args.saveLocation + '/html outputs/symbiodiniumTyping_OUTPUT_' + time.strftime("%H%M%S") + '_' + time.strftime("%d%m%y") +'.html', htmlfile)

## MAIN FUNCTION ##
def CreateHumeFstMatrices():
    #For initial use of program
    config.__init__()

    assignCladeCollections()
    assignInitialTypes()  # This now assigns initial types to the samples' cladeCollections
    print('Assigned initial types')

    #Create masterSeqDistancesDict
    masterSeqDistancesDict = None
    if config.args.createMasterSeqDistancesFromScratch == False:
        try:
            masterSeqDistancesDict = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'\serialized objects', 'masterSeqDistancesDict')
        except:
            messagebox.showwarning('Missing Object', 'masterSeqDistancesDict object not found in specified directory\n Creating from scratch...')
    if masterSeqDistancesDict == None:
        masterSeqDistances = createMasterSeqDistances()
        masterSeqDistancesDict = {frozenset(a.split(' ')[:-1]): float(a.split(' ')[2]) for a in masterSeqDistances}
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'\serialized objects','masterSeqDistancesDict', masterSeqDistancesDict)

    # Create log transformed distance column items of Fst distances between each sample within clade
    logTransedColDists = None
    if config.args.createFstColDistsFromScratch == False:
        try:
            logTransedColDists = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'\serialized objects',
                                                                    'logTransedColDists')
        except:
            messagebox.showwarning('Missing Object',
                                   'logTransedColDists object not found in specified directory\n Creating from scratch...')
    if logTransedColDists == None:
        logTransedColDists = createLogTransedFstColDists(masterSeqDistancesDict)
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'\serialized objects', 'logTransedColDists', logTransedColDists)

    print('Writing output')
    # Create the main html List for out put
    htmlOutput = createHtmlHolder()
    # Write out the graphics so that the writeTypeBasedOutput can write them into the html
    #htmlOutput.extend(writeMajMatricesToFile(logTransedColDists)) # This produces an html table of the plots
    writeMajMatricesToFile(logTransedColDists)
    # Append information to the samples within the abundanceList refering to the finalTypeCladeCollections
    inferFinalSymbiodiniumTypes()

    # Write the type based output
    absoluteSupportOfInitialTypeDict, initialSupportOfInitialTypesDict, htmlMainString = writeTypeBasedOutput()  # Produce and write type info output also return a dictionary that has the final support for the types
    htmlOutput.extend(htmlMainString)
    # Write the sample based output
    htmlOutput.extend(writeSampleCharacterisationOutput(absoluteSupportOfInitialTypeDict, initialSupportOfInitialTypesDict))
    # Add the closing tags to the html file
    closeAndWriteHtmlHolder(htmlOutput)
    print('Program complete')

## MAIN ENTRY POINT OF PROGRAM ##
if __name__ == '__main__':
    CreateHumeFstMatrices()


