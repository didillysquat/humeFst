# Imports

import numpy as np
from tkinter.filedialog import asksaveasfilename
import re
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
import shutil
import zipfile
from collections import Counter
import timeit
from sklearn import manifold
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

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
    # two abundances together and then multiply that result with the distance for that pair.

    # Similar to the between calculations we need to make sure that the logged abundances are re-normalised to 1000 for each sample
    # First calculate average for sampleone

    # I have replaced the log transformation with the cube root transformation because the log will give a negative when numbers are small e.g. below 1

    loggedTotalOne = sum([math.pow(occurance.abundance, 1/2) for occurance in sampleoneits2occurances])
    loggedTotalTwo = sum([math.pow(occurance.abundance, 1/2) for occurance in sampletwoits2occurances])

    for occuranceOne, occuranceTwo in itertools.combinations(sampleoneits2occurances, 2): # This gives us the combos that are possible. if only one type then this may fail but we can sort that easily
        runningTotalOne += (float((math.pow(occuranceOne.abundance, 1/2)/loggedTotalOne)*1000))*(float((math.pow(occuranceTwo.abundance, 1/2)/loggedTotalTwo)*1000))*float(seqdistances[frozenset({occuranceOne.name, occuranceTwo.name})])
    # Second calculate average for sampletwo
    for occuranceOne, occuranceTwo in itertools.combinations(sampletwoits2occurances, 2): # This gives us the combos that are possible. if only one type then this may fail but we can sort that easily
        # print((float((math.pow(occuranceOne.abundance, 1/3)/loggedTotalOne)*1000))*(float((math.pow(occuranceTwo.abundance, 1/3)/loggedTotalTwo)*1000))*float(seqdistances[frozenset({occuranceOne.name, occuranceTwo.name})]))
        runningTotalTwo += (float((math.pow(occuranceOne.abundance, 1/2)/loggedTotalOne)*1000))*(float((math.pow(occuranceTwo.abundance, 1/2)/loggedTotalTwo)*1000))*float(seqdistances[frozenset({occuranceOne.name, occuranceTwo.name})])


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
    loggedTotalOne = sum([math.pow(occurance.abundance, 1/2) for occurance in sampleoneits2occurances])
    loggedTotalTwo = sum([math.pow(occurance.abundance, 1/2) for occurance in sampletwoits2occurances])

    runningTotal = 0.00
    for oneOccurs in sampleoneits2occurances:
        for twoOccurs in sampletwoits2occurances:
            if oneOccurs.name != twoOccurs.name:#we only need to calculate distances if they are not the same
                runningTotal += seqdistances[frozenset({oneOccurs.name, twoOccurs.name})]*(float((math.pow(oneOccurs.abundance, 1/2)/loggedTotalOne)*1000))*(float((math.pow(twoOccurs.abundance, 1/2)/loggedTotalTwo)*1000))

    totalSeqCombos = 1000000

    return runningTotal/totalSeqCombos

def createMatrixFromColDists(Fstcoldistdict):
    FstMatrixCollection = []
    for CLADE in config.args.cladeList:
        listOfSamples = [] # [[cladeCollection.foundWithinSample for cladeCollection in sample.cladeCollectionList if cladeCollection.clade == CLADE] for sample in abundanceList] # All samples whether they have a type defined or not, then we can display the unsupported in R
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
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
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
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

def writeByteObjectToDefinedDirectory(directory,object):
    f = open(directory , 'wb+')
    pickle.dump(object, f)

def readByteObjectFromDefinedDirectory(directory):
    f = open(directory,'rb')
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

def writeTypeBasedOutput2Old():
    print('Conducting type-based analysis')
    #23/08/16
    #TODO
    # We have replaced the need for both of the inital support and final support dicts by creating the typeDB
    # We can now use len(typeDB[type].foundininitial)

    # We will go through each of the samples' cladeCollections by clade.
    # We will count how many samples each type is found in both as an initial and a final
    initialSupportDict = {}
    finalSupportDict = {}
    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        # Count support for inital types both unsupported and supported
        # We need to have a count for the unsupported types too as these will still be written out
        # when we are writing the sample based output
        for CLADECOLLECTION in SAMPLE.cladeCollectionList:
            if CLADECOLLECTION.initialType.name not in initialSupportDict.keys():
                initialSupportDict[CLADECOLLECTION.initialType.name] = 1
            else:
                initialSupportDict[CLADECOLLECTION.initialType.name] += 1


        # Count support for final types
        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
            for FINALTYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
                if FINALTYPE.name not in finalSupportDict.keys():
                    finalSupportDict[FINALTYPE.name] = 1
                else:
                    finalSupportDict[FINALTYPE.name] += 1

    # Once we have done all of the counts, we must also make sure that inital supported
    # types that have not been carried over into final types have a final type support of zero
    # E.g. if C1 as an itial type was surpassed by C1-XXX in all final samples, then C1 would not be in the final support dict
    # so we must add it in with a value of 0
    for initalType in initialSupportDict.keys():
        if initalType not in finalSupportDict.keys():
            finalSupportDict[initalType] = 0


    # Now collect data per clade to print out in the HTML, including:
    # Total samples
    # Number of types
    # Number of codom Types
    # Number of non-CoDom Majs


    #Doc for the output
    outPutDoc = []

    # Go through config.abundance list and populate the dicts and lists
    for CLADE in config.args.cladeList:
        totalSamples = 0
        collectionOfTypes = {}
        collectionOfMajs = {}
        collectionOfCoDomTypes = {}
        collectionOfNonCoDomMajs = {}
        collectionOfIdentifiedSamples = []
        cladeSpecificInitialSupportDict = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE and CLADECOLLECTION.initialType.supportedType:
                    if CLADECOLLECTION.initialType.name not in cladeSpecificInitialSupportDict.keys():
                        cladeSpecificInitialSupportDict[CLADECOLLECTION.initialType.name] = 1
                    else:
                        cladeSpecificInitialSupportDict[CLADECOLLECTION.initialType.name] += 1
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == CLADE:
                    totalSamples += 1
                    if FINALTYPECLADECOLLECTION.identified:
                        collectionOfIdentifiedSamples.append(SAMPLE.name)
                        for FINALTYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
                            if FINALTYPE.maj not in collectionOfMajs.keys():
                                collectionOfMajs[FINALTYPE.maj] = [FINALTYPE.name]
                            else: # Update the dictionary associating the new sample with the maj
                                tempList = collectionOfMajs[FINALTYPE.maj]
                                tempList.append(FINALTYPE.name)
                                collectionOfMajs[FINALTYPE.maj] = tempList

                            if FINALTYPE.name not in collectionOfTypes.keys():
                                collectionOfTypes[FINALTYPE.name] = 1
                            else:
                                collectionOfTypes[FINALTYPE.name] = collectionOfTypes[FINALTYPE.name] + 1
                            if FINALTYPE.coDom == True:
                                if FINALTYPE.name not in collectionOfCoDomTypes.keys():
                                    collectionOfCoDomTypes[FINALTYPE.name] = 1
                                else:
                                    collectionOfCoDomTypes[FINALTYPE.name] = collectionOfCoDomTypes[FINALTYPE.name] + 1
                            else:
                                if FINALTYPE.maj not in collectionOfNonCoDomMajs:
                                    collectionOfNonCoDomMajs[FINALTYPE.maj] = 1
                                else:
                                    collectionOfNonCoDomMajs[FINALTYPE.maj] = collectionOfNonCoDomMajs[FINALTYPE.maj] + 1


        # At this point we have all of the counts and collections completed.
        # Now use them, still within the CLADE loop to write out the information
        # in the correct order to the outPutDoc
        # This will then be translated into html by Jinja at the end
        # The output doc contains lines that are tuples with a list of information
        # to put in the html and a string that tells jinja how to treat the info
        orderedListOfNonCoDomMajs = [a[0] for a in sorted(collectionOfNonCoDomMajs.items(), key=lambda x: x[1], reverse=True)]
        orderedListOfCoDomTypes = [a[0] for a in sorted(collectionOfCoDomTypes.items(), key=lambda x: x[1], reverse=True)]
        if len(orderedListOfNonCoDomMajs) > 0 or len(orderedListOfCoDomTypes) > 0:

            # ADD CLADE LEVEL INFOMATION
            outPutDoc.append(('Clade {0}: Types = {1} / non-CoDom Majs = {2} / coDom types = {3} / samples with type identified = {4}'.format(CLADE, str(len(collectionOfTypes.keys())), str(len(collectionOfNonCoDomMajs.keys())), str(len(collectionOfCoDomTypes)), str((totalSamples - len(collectionOfIdentifiedSamples)))), 'cladeHeader'))

            # INSERT CLADAL TYPE SUMMARY INFORMATION
            outPutDoc.append(('Type Summary', 'majHeader'))
            # We currently list the types in order of initial support. However we could
            # just as easily list them in final support just by using the finalSupportDict instead
            orderedListOfInitialTypes = [a[0] for a in sorted(cladeSpecificInitialSupportDict.items(), key=lambda x: x[1], reverse=True)]
            for INITIALTYPE in orderedListOfInitialTypes:
                outPutDoc.append((['{0}'.format(INITIALTYPE), '{0}:{1}'.format(str(initialSupportDict[INITIALTYPE]),
                                                                               str(finalSupportDict[INITIALTYPE]))],
                                  'typeInfo'))

            # INSERT FIGURE HERE
            # If the plots that we want to add exists
            fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/Clade' + CLADE + '_FstPlot.svg'
            paretoPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/Clade' + CLADE + '_pareto.svg'
            if os.path.isfile(fstPlotDir) and os.path.isfile(paretoPlotDir):
                outPutDoc.append(([fstPlotDir, paretoPlotDir], 'cladalPlot'))
            else:
                outPutDoc.append((config.args.saveLocation + '/html templates/image.jpg', 'plotError'))


            # NOW GO MAJ BY MAJ WITHIN THE CLADE
            orderedListOfMajs = [a[0] for a in sorted(collectionOfMajs.items(), key=lambda x: len(x[1]), reverse=True)]
            for MAJ in orderedListOfMajs:
                # if MAJ in config.oursToLaJDict:
                #     convertedMaj = config.oursToLaJDict[MAJ]
                # else:
                #     convertedMaj = MAJ
                convertedMaj = CLJ(MAJ)

                # CALCULATE NUMBER OF CODOM TYPES IN MAJ this therefore allows us to calculate non-codoms too as we know the total numbe of types.
                numberOfCoDomTypesInMaj = 0
                for typeName in set(collectionOfMajs[MAJ]):
                    if typeName in orderedListOfCoDomTypes:
                        numberOfCoDomTypesInMaj += 1

                # INSERT MAJ SUMMARY
                outPutDoc.append(('MajITS2Seq ' + convertedMaj + ': non-CoDom types = ' + str(len(set(collectionOfMajs[MAJ]))-numberOfCoDomTypesInMaj) + ' / CoDom types = ' + str(numberOfCoDomTypesInMaj) + ' / Total types = ' + str(len(set(collectionOfMajs[MAJ]))) + ' / Samples = ' + str(len(collectionOfMajs[MAJ])), 'majHeader'))

                # INSERT MAJ PLOT
                if MAJ in config.oursToLaJDict.keys():
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + CLJ(MAJ) + '/' + CLJ(MAJ) + '_FstPlot.svg'
                else:
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + CLJ(MAJ) + '/Clade C ' + CLJ(MAJ) + '_FstPlot.svg'
                if os.path.isfile(fstPlotDir):
                    outPutDoc.append((fstPlotDir, 'majPlot'))
                else:
                    outPutDoc.append(('Insufficient samples to warrant plot', 'plotError'))

                # TABLE HEADERS FOR MAJ SUMMARY NON_CODOM
                outPutDoc.append((['non-CoDom type', 'Initial:Final Support'], 'subHeader'))

                # Make list of most common types using the Counter class
                # http://stackoverflow.com/questions/3594514/how-to-find-most-common-elements-of-a-list
                listOfNonCoDomTypesInMaj = [typeName for typeName in collectionOfMajs[MAJ] if typeName not in collectionOfCoDomTypes.keys()]
                orderedListOfAbundanceTuplesTypesInMaj = Counter(listOfNonCoDomTypesInMaj).most_common(len(set(listOfNonCoDomTypesInMaj)))
                orderedListOfTypesInMaj = [a[0] for a in orderedListOfAbundanceTuplesTypesInMaj]

                # GO TYPE BY TYPE WITHIN THE MAJ NON-CODOM FIRST THEN CODOM BOTH IN ORDER OF ABUNDACNE
                for TYPENAME in orderedListOfTypesInMaj:
                    outPutDoc.append((['Type ' + TYPENAME, str(initialSupportDict[TYPENAME]) + ':' + str(finalSupportDict[TYPENAME])], 'typeInfo'))

                #CHECK TO SEE IF THERE ARE CODOMS TO ADD
                listOfCoDomTypesInMaj = [typeName for typeName in collectionOfMajs[MAJ] if typeName in collectionOfCoDomTypes.keys()]
                orderedListOfAbundanceCoDomTuplesTypesInMaj = Counter(listOfCoDomTypesInMaj).most_common(len(set(listOfCoDomTypesInMaj)))
                orderedListOfCoDomTypesInMaj = [a[0] for a in orderedListOfAbundanceCoDomTuplesTypesInMaj]
                if len(orderedListOfCoDomTypesInMaj) > 0:
                    # CODOM HEADER
                    outPutDoc.append(('Co-Dominant Types = ' + str(len(orderedListOfCoDomTypesInMaj)), 'majHeader'))
                    # ADD CODOMS
                    for TYPENAME in orderedListOfCoDomTypesInMaj:
                        outPutDoc.append((['Type {0}'.format(TYPENAME), '{0}:{1}'.format(str(initialSupportDict[TYPENAME]), str(finalSupportDict[TYPENAME]))], 'typeInfo'))

    # TRANSLATE OUTPUT DOC INTO HTML
    jinjaEnvironment = Environment(loader=FileSystemLoader(config.args.rootLocation + r'\html templates'), trim_blocks=True)
    jinjaTemplate = jinjaEnvironment.get_template('typeCharacterisation__TEMPLATE.html')
    stringThing = jinjaTemplate.render(outPut=outPutDoc)
    htmlString = [a for a in stringThing.split('\n')]

    print('Completed writeTypeBasedOutput()')
    return initialSupportDict, finalSupportDict, htmlString

def writeTypeBasedOutput2():
    print('Conducting type-based analysis')


    #Doc for the output
    outPutDoc = []

    # Go through config.abundance list and populate the dicts and lists
    for CLADE in config.args.cladeList:
        listOfTypesInClade = [config.typeDB[a] for a in config.typeDB.keys() if config.typeDB[a].clade == CLADE]
        orderedListOfTypesInClade = sorted(listOfTypesInClade, key=lambda x: len(x.samplesFoundInAsInitial), reverse=True)
        listOfMajsInClade = list(set([val for sublist in [list(set(config.typeDB[a].majList)) for a in config.typeDB.keys() if config.typeDB[a].clade == CLADE] for val in sublist]))
        listOfCoDomTypesInClade = [config.typeDB[a] for a in config.typeDB.keys() if config.typeDB[a].coDom]
        #TODO I reckon that this dictOfMajAbundancesInClade is basically what we want to plot
        dictOfMajAbundancesInClade = {}
        for TYPE in listOfTypesInClade:
            for maj in set(TYPE.majList):
                addItemPairToDict(maj, TYPE.majList.count(maj), dictOfMajAbundancesInClade)
        orderedListOfMajAbundancesInClade = [a[0] for a in sorted(dictOfMajAbundancesInClade.items(), key=lambda x: x[1], reverse=True)]



        if len(listOfTypesInClade):

            # ADD CLADE LEVEL INFOMATION
            outPutDoc.append(('Clade {0}: Types = {1} / Majs = {2} / coDom types = {3}'.format(CLADE, str(len(listOfTypesInClade)), str(len(listOfMajsInClade)), str(len(listOfCoDomTypesInClade))), 'cladeHeader'))

            # INSERT CLADAL TYPE SUMMARY INFORMATION
            outPutDoc.append(('Type Summary', 'majHeader'))
            # We currently list the types in order of initial support. However we could
            # just as easily list them in final support just by using the finalSupportDict instead

            for INITIALTYPE in orderedListOfTypesInClade:
                if INITIALTYPE.samplesFoundInAsFinal:
                    outPutDoc.append((['{0}'.format(INITIALTYPE.name), '{0}:{1}'.format(str(len(INITIALTYPE.samplesFoundInAsInitial)),
                                                                                        str(len(INITIALTYPE.samplesFoundInAsFinal)))],
                                      'typeInfo'))

            # INSERT FIGURE HERE
            # If the plots that we want to add exists
            fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/Clade' + CLADE + '_FstPlot.svg'
            paretoPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/Clade' + CLADE + '_pareto.svg'
            if os.path.isfile(fstPlotDir) and os.path.isfile(paretoPlotDir):
                outPutDoc.append(([fstPlotDir, paretoPlotDir], 'cladalPlot'))
            else:
                outPutDoc.append((config.args.saveLocation + '/html templates/image.jpg', 'plotError'))


            # NOW GO MAJ BY MAJ WITHIN THE CLADE

            for MAJ in orderedListOfMajAbundancesInClade:

                coDomTypesInMaj = []
                nonCoDomTypesInMaj = []
                for typesInClade in listOfTypesInClade:
                    if MAJ in typesInClade.majList:
                        if typesInClade.coDom:
                            coDomTypesInMaj.append(typesInClade)
                        else:
                            nonCoDomTypesInMaj.append(typesInClade)


                # if MAJ in config.oursToLaJDict:
                #     convertedMaj = config.oursToLaJDict[MAJ]
                # else:
                #     convertedMaj = MAJ
                convertedMaj = CLJ(MAJ)
                # INSERT MAJ SUMMARY
                outPutDoc.append(('MajITS2Seq {0}: non-CoDom types = {1} / CoDom types = {2} / Total types = {3}'.format(convertedMaj, str(len(nonCoDomTypesInMaj)), str(len(coDomTypesInMaj)), str(len(coDomTypesInMaj) + len(nonCoDomTypesInMaj))), 'majHeader'))

                # INSERT MAJ PLOT
                if MAJ in config.oursToLaJDict.keys():
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + CLJ(MAJ) + '/' + CLJ(MAJ) + '_FstPlot.svg'
                else:
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + MAJ + '/Clade C ' + MAJ + '_FstPlot.svg'
                if os.path.isfile(fstPlotDir):
                    outPutDoc.append((fstPlotDir, 'majPlot'))
                else:
                    outPutDoc.append(('Insufficient samples to warrant plot', 'plotError'))

                # TABLE HEADERS FOR MAJ SUMMARY NON_CODOM
                outPutDoc.append((['non-CoDom type', 'Initial:Final Support'], 'subHeader'))

                # Make list of most common types using the Counter class
                # http://stackoverflow.com/questions/3594514/how-to-find-most-common-elements-of-a-list
                orderedListOfNonCoDomTypesInMaj = sorted(nonCoDomTypesInMaj, key=lambda x: len(x.samplesFoundInAsInitial), reverse=True)
                orderedListOfCoDomTypesInMaj = sorted(coDomTypesInMaj, key=lambda x: len(x.samplesFoundInAsInitial), reverse=True)
                # GO TYPE BY TYPE WITHIN THE MAJ NON-CODOM FIRST THEN CODOM BOTH IN ORDER OF ABUNDACNE
                for TYPE in orderedListOfNonCoDomTypesInMaj:
                    outPutDoc.append((['Type: {0}'.format(TYPE.name), '{0}:{1}'.format(str(len(TYPE.samplesFoundInAsInitial)), str(len(TYPE.samplesFoundInAsFinal)))], 'typeInfo'))

                #CHECK TO SEE IF THERE ARE CODOMS TO ADD
                if orderedListOfCoDomTypesInMaj:
                    # CODOM HEADER
                    outPutDoc.append(('Co-Dominant Types = {0}'.format(str(len(orderedListOfCoDomTypesInMaj))), 'majHeader'))
                    # ADD CODOMS
                    for TYPE in orderedListOfCoDomTypesInMaj:
                        outPutDoc.append((['Type: {0}'.format(TYPE.name), '{0}:{1}'.format(str(len(TYPE.samplesFoundInAsInitial)), str(len(TYPE.samplesFoundInAsFinal)))], 'typeInfo'))

    # TRANSLATE OUTPUT DOC INTO HTML
    jinjaEnvironment = Environment(loader=FileSystemLoader(config.args.rootLocation + r'/html templates'), trim_blocks=True)
    jinjaTemplate = jinjaEnvironment.get_template('typeCharacterisation__TEMPLATE.html')
    stringThing = jinjaTemplate.render(outPut=outPutDoc)
    htmlString = [a for a in stringThing.split('\n')]

    print('Completed writeTypeBasedOutput()')
    return htmlString


def addItemPairToDict(KEY, VALUE, DICT):
    if KEY in DICT.keys():
        DICT[KEY] += VALUE
    else:
        DICT[KEY] = VALUE
    return DICT


def writeTypeBasedOutput(): # Produce all of the info such as number of Maj ITS2's per clade, number of clades, number of codoms, number of predefined types
    print('Running writeTypeBasedOutput()')

    # THis had a problem in that it might show when there is increased support for a type but it doesn't work out when support has decreased for a type. For example, a lot of th
    # samples that had initial type A1 will have alternative final types. so the final support for A1 will drop when there is an alternative
    # In order to reflect this we will make an finalSupportOfInitialTypeDict which will keep track of the changes in support.
    # If a finaltypelist contains a type then this causes an increase of 1 in the support for the initial type.
    # If a finalType list has any types in it then the inital type of that cladecollection will get a -1 support as alterantives have been found.
    # If the finalType list is empty then that means that the support will be unchanged.
    # Unfortunately for the types that have only one defining intra, for example A1 they cannot gain support as you cannot have a finaltype of A1,
    # As your initial type would have had to have been A1 and therefore it wouldn't be in the final types list.


    #23/08/16
    # We have left the inital type in the finalTypeList as we were not sure why it was being removed.
    # Again, I think the above can be simplified.
    # TODO go back and leave out the clause that means only initially supported types are forced in to look for final types
    # We will go through each of the samples' cladeCollections by clade.
    # We will do a count for the intial types and we will do a count for the final types.
    # We will conly count final types if the finaltype is .identified==True


    # Create and populate the finalSupportOfInitialTypeDict
    # The values from this dict will then be used by adding them (some values will be negative) to the intial support to get the final support
    changeInFinalSupportOfInitialTypeDict = {}
    initialSupportOfInitialTypesDict = {} # This will be used in the next method key = footprint, value = initial support
    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
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

        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
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
                for SAMPLEKEY in config.abundanceList.keys():
                    SAMPLE = config.abundanceList[SAMPLEKEY]
                    for cladeCollection in SAMPLE.cladeCollectionList:
                        if not cladeCollection.initialType.coDom and cladeCollection.initialType.maj == MAJ and cladeCollection.initialType.supportedType:
                            listOfNonCoDomTypesWithinMaj.append(cladeCollection.initialType)

                # Will now produce list of types which are in order of their abundance within the MAJ
                #orderedlistOfNonCoDomTypesWithinMaj = [a[0] for a in sorted({unique: listOfNonCoDomTypesWithinMaj.count(unique) for unique in set(listOfNonCoDomTypesWithinMaj)}.items())]
                orderedlistOfNonCoDomTypesWithinMaj = [a[0] for a in sorted({unique: listOfNonCoDomTypesWithinMaj.count(unique) for unique in set(listOfNonCoDomTypesWithinMaj)}.items(), key=lambda x: x[1], reverse=True)]

                if MAJ in config.oursToLaJDict.keys():
                    outPutDoc.append(('MajITS2Seq ' + CLJ(MAJ) + ': Number of Types = ' + str(len(orderedlistOfNonCoDomTypesWithinMaj)) + ' / Number of Samples = ' + str(len(listOfNonCoDomTypesWithinMaj)), 'majHeader'))
                else:
                    outPutDoc.append(('MajITS2Seq ' + MAJ[3:] + ': Number of Types = ' + str(len(orderedlistOfNonCoDomTypesWithinMaj)) + ' / Number of Samples = ' + str(len(listOfNonCoDomTypesWithinMaj)), 'majHeader'))
                #INSERT PLOT HERE
                if MAJ in config.oursToLaJDict.keys():
                    fstPlotDir = config.args.saveLocation + '/matrices outputs/Clade' + CLADE + '/IncludingMaj/' + CLJ(MAJ) + '/' + CLJ(MAJ) + '_FstPlot.svg'
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

    jinjaEnvironment = Environment(loader=FileSystemLoader(config.args.rootLocation + r'\html templates'), trim_blocks=True)
    jinjaTemplate = jinjaEnvironment.get_template('typeCharacterisation__TEMPLATE.html')
    stringThing = jinjaTemplate.render(outPut=outPutDoc)
    htmlString = [a for a in stringThing.split('\n')]
    if not os.path.exists(config.args.saveLocation + '/html outputs'):
        os.makedirs(config.args.saveLocation + '/html outputs')
    #writeListToDestination(config.args.saveLocation +'/html outputs/typeCharacterisation__OUTPUT.html', htmlString)

    print('Completed writeTypeBasedOutput()')
    return changeInFinalSupportOfInitialTypeDict, initialSupportOfInitialTypesDict, htmlString


def producePlotsold(coldists):
    # This script makes a plot for each clade level and each subclade or maj level dependent on whether there are enough samples to warrant it
    # It sends info to writeDataForMakingPlots to produce the actual matrix files ect that R will turn into the plot.
    print('Producing Plots')
    argStringCSVList = []
    for CLADE in config.args.cladeList:
        #### STEP ONE
        # Firstly write out a cladal Matrice before looking to do the Maj-based matrices
        listOfSamples = [] # A list of samples with cladeCollections of the given clade. i.e. all samples that will be represented in the majPlot
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            # for CLADECOLLECTION in SAMPLE.cladeCollectionList:
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                # if CLADECOLLECTION.clade == CLADE:
                if FINALTYPECLADECOLLECTION.clade == CLADE:
                    listOfSamples.append(FINALTYPECLADECOLLECTION.foundWithinSample)
        # This crazy one liner produces a lists of lists(the interal list) which is then flattened to a single list


        if len(listOfSamples) >= config.args.cladePlotCutOff: # There will be some clades that don't have any clade collections in them. In this case we can't make a graphic as there are no samples.
            # writeDataForMakingPlots returns a csv string that has the directories where the matrix, listOfSpecies and info files are located
            # R will make plots using these
            argStringCSVList.append(writeDataForMakingPlots(listofsamples=listOfSamples, clade=CLADE, Fstcoldistdict=coldists, ismaj=False))


            #### STEP TWO
            # Then write out the Maj matrices if appropriate

            # I am going to collect the Majs for each cladeCollection in each sample irrespective of whether they are a codom or not.
            # I will do this by creating one dictionary which has the Maj as the key and a list of samples that have that as a Maj as the value
            majToSamplesDict = {}
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE:
                        if SAMPLE.name == 'OMd_028' and CLADE == 'D':
                            a = 4
                        # Here I am also going to check that it has a finaltypecladecollection and if it does that it is .identified:
                        # That is to say that if there is no final type identified, or if there is a final type identified that is .identified = false
                        # i.e. all of the intras above the cutoff are not explained by the types identified.
                        # If samples which don't have a type identified for the given cladecollection will not be plotted.
                        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList: # This will ony work if there is a FINALTYPECLADECOLLECTION
                            if FINALTYPECLADECOLLECTION.clade == CLADE and FINALTYPECLADECOLLECTION.identified:
                                    if CLADECOLLECTION.maj in majToSamplesDict.keys(): # if the maj is already in the dict add sample to current list
                                        tempList = majToSamplesDict[CLADECOLLECTION.maj]
                                        tempList.append(SAMPLE.name)
                                        majToSamplesDict[CLADECOLLECTION.maj] = tempList
                                    else: # If maj not in dict yet then create new entry with the current sample's name
                                        majToSamplesDict[CLADECOLLECTION.maj] = [SAMPLE.name]

            # For each of the MAJs identified above create a plot if above the plotting threshold
            for MAJ in majToSamplesDict.keys(): # for each Maj ITS2 that isn't ONLY found in coDoms and has more than the threshold for supported types


                # If possible convert to the LaJ name
                if MAJ in config.oursToLaJDict.keys():
                    convertedMaj = CLJ(MAJ)
                    isConverted = True
                else:
                    isConverted = False
                    convertedMaj = MAJ

                if len(majToSamplesDict[MAJ]) >= config.args.majPlotCutOff: # Only make a MajPlot if there are more than or equal to the maj plot threshold cutoff number of sampels worth plotting (coDOm and non-CoDom types)

                    argStringCSVList.append(writeDataForMakingPlots(Fstcoldistdict=coldists, listofsamples=majToSamplesDict[MAJ], clade=CLADE, ismaj=True, maj=convertedMaj, isconverted=isConverted))

    # Here write out the argstringcsvlist to file
    writeListToDestination(config.args.saveLocation + r'\matrices outputs\RArgList.txt', argStringCSVList)


    producePlotWithR()

    return


def producePlots(coldists):
    # This script makes a plot for each clade level and each subclade or maj level dependent on whether there are enough samples to warrant it
    # It sends info to writeDataForMakingPlots to produce the actual matrix files ect that R will turn into the plot.
    print('Producing Plots')
    argStringCSVList = []
    for CLADE in config.args.cladeList:
        #### STEP ONE
        # Firstly write out a cladal Matrice before looking to do the Maj-based matrices
        listOfSamples = [val for sublist in [config.typeDB[a].samplesFoundInAsFinal for a in config.typeDB.keys() if
                                             config.typeDB[a].clade == CLADE] for val in sublist]
        if len(listOfSamples) >= config.args.cladePlotCutOff:  # There will be some clades that don't have any clade collections in them. In this case we can't make a graphic as there are no samples.
            # writeDataForMakingPlots returns a csv string that has the directories where the matrix, listOfSpecies and info files are located
            # R will make plots using these
            argStringCSVList.append(
                writeDataForMakingPlots(listofsamples=listOfSamples, clade=CLADE, Fstcoldistdict=coldists, ismaj=False))

            #### STEP TWO
            # Then write out the Maj matrices if appropriate
            # I am going to collect the Majs for each cladeCollection in each sample irrespective of whether they are a codom or not.
            # I will do this by creating one dictionary which has the Maj as the key and a list of samples that have that as a Maj as the value
            majToSamplesDict = {}
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE:
                                if CLADECOLLECTION.maj in majToSamplesDict.keys():  # if the maj is already in the dict add sample to current list
                                    majToSamplesDict[CLADECOLLECTION.maj].append(SAMPLE.name)
                                else:  # If maj not in dict yet then create new entry with the current sample's name
                                    majToSamplesDict[CLADECOLLECTION.maj] = [SAMPLE.name]

            # For each of the MAJs identified above create a plot if above the plotting threshold
            for MAJ in majToSamplesDict.keys():  # for each Maj ITS2 that isn't ONLY found in coDoms and has more than the threshold for supported types
                if len(majToSamplesDict[MAJ]) >= config.args.majPlotCutOff:  # Only make a MajPlot if there are more than or equal to the maj plot threshold cutoff number of sampels worth plotting (coDOm and non-CoDom types)
                    convertedMaj, isConverted = convertToLaJ(MAJ)

                    argStringCSVList.append(
                        writeDataForMakingPlots(Fstcoldistdict=coldists, listofsamples=majToSamplesDict[MAJ],
                                                clade=CLADE, ismaj=True, maj=convertedMaj, isconverted=isConverted))

    # Here write out the argstringcsvlist to file
    writeListToDestination(config.args.saveLocation + r'\matrices outputs\RArgList.txt', argStringCSVList)

    producePlotWithR()

    return

def convertToLaJ(VAR):
    try:
        if len(config.oursToLaJDict[VAR]) > 1:
            return '#'.join(config.oursToLaJDict[VAR]), True
        else:
            return config.oursToLaJDict[VAR][0], True
    except:
        return VAR, False


def writeDataForMakingPlots(Fstcoldistdict, listofsamples, clade, ismaj,  maj = None,  isconverted = None):
    # This function will make the matrix, and information file and write them along with the list of samples to file.
    # They will be used by R to create the plots for the html output
    # It returns argString, which is csv string that will be read by R and used to locate the data to make the plots from

    #######
    #  Step ONE
    # Create the matrix file
    # Create empty matrix that is of dimensions [rows(n samples + 1)][columns(number of samples + 1)]
    listOfSamples = listofsamples
    xstr = lambda s: s or ""
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

    ##########
    # STEP TWO
    # Write the files to destination. Destinations are different on whether this is a clade matrix or a subclade (e.g. maj) matrix
    # So we do them sepearately depending on ismaj boolean variable
    if ismaj: # If it is a sub cladal Maj level Fst
        # A) Check to see if the save directory already exists, if not, create it.
        if not os.path.exists(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/IncludingMaj/' + maj):
            os.makedirs(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/IncludingMaj/' + maj)
        # B) Write the matrix to file
        write2DListToDestination(config.args.saveLocation + '\\matrices outputs\\Clade' + clade + '\\' + 'IncludingMaj' + '\\' + maj + '\\Matrix.dist', FstMatrix)
        # C) Write the listOfSamples to file
        writeListToDestination(config.args.saveLocation + '\\matrices outputs\\Clade' + clade + '\\' + 'IncludingMaj' + '\\' + maj + '\\Header.head', listOfSamples)

        # D) Create the info list
        infoList = []
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name in listOfSamples:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                            # if CLADECOLLECTION.maj in config.oursToLaJDict:
                            #     majITS2 = config.oursToLaJDict[CLADECOLLECTION.maj]
                            # else:
                            #     majITS2 = CLADECOLLECTION.maj
                            majITS2 = CLJ(CLADECOLLECTION.maj)
                            infoList.append(','.join([SAMPLE.name, xstr(SAMPLE.hostTaxon), xstr(SAMPLE.reef), xstr(SAMPLE.region), majITS2, CLADECOLLECTION.initialType.name, str(CLADECOLLECTION.initialType.coDom)]))
        # E) Write info list
        writeListToDestination(config.args.saveLocation + '\\matrices outputs\\Clade' + clade + '\\' + 'IncludingMaj' + '\\' + maj + '\\SampleInfo.info', infoList)

        # F) Add directories holding data for reading into R
        if isconverted:
            argString = ','.join([ (config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + maj).replace('\\','/'), maj, str(ismaj).upper()])
        else:
            argString = ','.join([ (config.args.saveLocation + '/matrices outputs/Clade' + clade + '/' + 'IncludingMaj' + '/' + maj).replace('\\','/'), 'Clade ' + clade + ' ' + maj, str(ismaj).upper()])




    else: # If a cladal level Fst
        # A) Check to see if the save directory already exists, if not, create it.
        if not os.path.exists(config.args.saveLocation + '/matrices outputs/Clade' + clade):
            os.makedirs(config.args.saveLocation + '/matrices outputs/Clade' + clade)
        # B) Write matrix to file
        write2DListToDestination(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/Matrix.dist', FstMatrix)
        # C) Write list of samples to file
        writeListToDestination(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/Header.head', listOfSamples)

        # D) Create the info list which will contain data on sample name, host, reef, region etc.
        infoList = []
        for SAMPLE in [config.abundanceList[a] for a in listOfSamples]:
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == clade:
                        # Check to see if the Maj name can be LaJconverteds
                        # Before adding it to the infolist
                        # if CLADECOLLECTION.maj in config.oursToLaJDict:
                        #     majITS2 = config.oursToLaJDict[CLADECOLLECTION.maj]
                        # else:
                        #     majITS2 = CLADECOLLECTION.maj
                        majITS2 = CLJ(CLADECOLLECTION.maj)
                        if CLADECOLLECTION.initialType == None:
                            badger = 'cider'
                        infoList.append(','.join([SAMPLE.name, xstr(SAMPLE.hostTaxon), xstr(SAMPLE.reef), xstr(SAMPLE.region), majITS2, CLADECOLLECTION.initialType.name, str(CLADECOLLECTION.initialType.coDom)]))

        # E) Write infolist
        writeListToDestination(config.args.saveLocation + '/matrices outputs/Clade' + clade + '/SampleInfo.info', infoList)

        # F) Add directories holding data for reading into R
        argString = ','.join([(config.args.saveLocation + '/matrices outputs/Clade' + clade ).replace('\\','/'), 'Clade' + clade, str(ismaj).upper()])

    print('Completed writeDataForMakingPlots()')

    return argString


def producePlotWithR():
    batFile = (config.args.rootLocation + r"\.bat scripts\updatedR.r").replace('\\','/')
    argDir = (config.args.saveLocation + r'\matrices outputs\RArgList.txt').replace('\\','/')
    # This R file can be run in cmd.exe by entering each string between quotation marks (double not single) and leaving a space between each argument.
    # Spaces in names are fine.
    cmd = [config.args.rscriptLocation, batFile, config.args.rootLocation, argDir]
    print(cmd)
    print('Running ' + str(cmd))
    call(cmd)

    return

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
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    listOfCladeCollections.append(CLADECOLLECTION)
        for a,b in itertools.combinations(listOfCladeCollections, 2): # For each paired combination of sample
            tempColDist.append([{a.foundWithinSample, b.foundWithinSample}, str(CalculateSinglePairwiseFst(a.listOfSeqsAboveCutOff, b.listOfSeqsAboveCutOff, masterSeqDistancesDict))])
        # Now convert the list structure to a dict structure
        tempDict = {frozenset(item[0]): item[1] for item in tempColDist}

        print('Clade ' + CLADE + ' Fst colDist complete')
        FstColDist.append(tempDict)
    print('createLogTransedFstColDists() complete')

    return FstColDist

def createLogTransedFstColDistsFromFinalTypes(masterSeqDistancesDict):

    allDefiningITS2OccurancesDict = {}
    print('Running createLogTransedFstColDistsFromFinalTypes()')
    FstColDist = []
    for CLADE in config.args.cladeList:  # For each clade
        tempColDist = []
        listOfFinalCladeCollections = []

        # Make list of finalcladecollection within clade
        # For each sample containing one, add their alldefiningits2occurances to the dict for use in the iter.
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for FINALCLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALCLADECOLLECTION.clade == CLADE:
                    allDefiningITS2OccurancesDict['{0}/{1}'.format(SAMPLE.name, CLADE)] = getAllDefiningIts2OccurancesOne(FINALCLADECOLLECTION)
                    listOfFinalCladeCollections.append(FINALCLADECOLLECTION)


        for a, b in itertools.combinations(listOfFinalCladeCollections, 2):  # For each paired combination of sample
            allDefiningIts2OccurancesA = allDefiningITS2OccurancesDict['{0}/{1}'.format(a.foundWithinSample, a.clade)]
            allDefiningIts2OccurancesB = allDefiningITS2OccurancesDict['{0}/{1}'.format(b.foundWithinSample, b.clade)]
            tempColDist.append([{a.foundWithinSample, b.foundWithinSample}, str(CalculateSinglePairwiseFst(allDefiningIts2OccurancesA, allDefiningIts2OccurancesB, masterSeqDistancesDict))])
        # Now convert the list structure to a dict structure
        tempDict = {frozenset(item[0]): item[1] for item in tempColDist}

        print('Clade ' + CLADE + ' Fst colDist complete')
        FstColDist.append(tempDict)
    print('createLogTransedFstColDistsFromFinalTypes() complete')

    return FstColDist

def getAllDefiningIts2Occurances(finalTypeCladeCollectionA, finalTypeCladeCollectionB):
    allDefiningIts2OccurancesA = []
    checkedOccurancesA = []
    allDefiningIts2OccurancesB = []
    checkedOccurancesB = []


    #First calculate the relative abundances of the defining sequences to each other
    for types in finalTypeCladeCollectionA.listOfFinalTypes:
        for sortedOccurance in types.sortedDefiningIts2Occurances:
            if sortedOccurance not in checkedOccurancesA:
                checkedOccurancesA.append(sortedOccurance)
    totalInA = sum(a[1] for a in checkedOccurancesA)

    for types in finalTypeCladeCollectionB.listOfFinalTypes:
        for sortedOccurance in types.sortedDefiningIts2Occurances:
            if sortedOccurance not in checkedOccurancesB:
                checkedOccurancesB.append(sortedOccurance)
    totalInB = sum(b[1] for b in checkedOccurancesB)

    checkedOccurancesA.clear()
    checkedOccurancesB.clear()

    for types in finalTypeCladeCollectionA.listOfFinalTypes:
        for sortedOccurance in types.sortedDefiningIts2Occurances:
            if sortedOccurance not in checkedOccurancesA:
                checkedOccurancesA.append(sortedOccurance)
                allDefiningIts2OccurancesA.append(its2SequenceOccurance(abundance = sortedOccurance[1]/(totalInA) * 1000, name = sortedOccurance[0], clade=None, sequence=None))

    for types in finalTypeCladeCollectionB.listOfFinalTypes:
        for sortedOccurance in types.sortedDefiningIts2Occurances:
            if sortedOccurance not in checkedOccurancesB:
                checkedOccurancesB.append(sortedOccurance)
                allDefiningIts2OccurancesB.append(its2SequenceOccurance(abundance=sortedOccurance[1]/(totalInB) * 1000, name=sortedOccurance[0], clade=None, sequence=None))

    return allDefiningIts2OccurancesA, allDefiningIts2OccurancesB

def getAllDefiningIts2OccurancesOne(finalTypeCladeCollectionA):
    ''' Here we are trying to identify the relative abundances (within the finalTYpeCladeCOllection)
    of all intras that are found within one of the final types'''
    allDefiningIts2OccurancesA = []

    checked = []
    totalInA = 0


    #Calculate the total reads from intras found in the final types collection
    SAMPLE = config.abundanceList[finalTypeCladeCollectionA.foundWithinSample]
    for types in finalTypeCladeCollectionA.sortedListOfFinalTypes:
        for intra in config.typeDB[types].footPrint:
            if intra not in [a[0] for a in checked]:
                checked.append((intra, SAMPLE.intraAbundanceDict[intra]))
                totalInA += SAMPLE.intraAbundanceDict[intra]



    # Calculate the relative abundaces of the intras
    for intra in checked:
        allDefiningIts2OccurancesA.append(its2SequenceOccurance(abundance=(intra[1]/totalInA)*1000, name=intra[0], clade=None, sequence=None))


    return allDefiningIts2OccurancesA

def assignInitialTypesOld(cladecollectioncountdict):
    cladeCollectionCountDict = cladecollectioncountdict
    #STEP ONE
    #Make a list of all of the unique footprints
    # Do this clade by clade
    for CLADE in config.args.cladeList:
        listOfFootprints  = []
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name == 'ADa_011' or SAMPLE.name == 'OMd_028':
                a = 1
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.footPrint not in listOfFootprints:
                        listOfFootprints.append(CLADECOLLECTION.footPrint)
        # STEP TWO
        # For each footprint, work out whether it is a coDom and whether it is a supported type
        # Reassign to the samples' cladeCollections that contain it.
        unsupportedTypeList = {}
        supportedList = []
        for FOOTPRINT in listOfFootprints:
            if FOOTPRINT == frozenset({'Otu0003'}):
                a=5
            coDom = False
            supportedType = False
            listOfSamplesThatContainFootprint = [] # List of the sample names that have a cladeCollection that match the footprint (exactly)
            listOfMajsForSamplesWithFootPrint = [] # List of the sequence names that are found as the predominant seqs in the footprint in question
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE and CLADECOLLECTION.footPrint == FOOTPRINT:
                        listOfSamplesThatContainFootprint.append(SAMPLE.name)
                        listOfMajsForSamplesWithFootPrint.append(CLADECOLLECTION.maj)

            # Ascertain whether this is a supported type and if so check to see if it is a coDom
            # There are enough samples with this footprint to define it as a type, continue checking it for being a coDom or if there is only one seq in the footprint then this is a Maj and therefore must be a type.
            if len(listOfSamplesThatContainFootprint) >= max(4, math.ceil(config.args.typeSupport*cladeCollectionCountDict[CLADE])) or len(FOOTPRINT) == 1:
                supportedType = True
                if len(set(listOfMajsForSamplesWithFootPrint)) > 1: # More than one Maj ITS2 seq: This is a coDom strictly speaking now we need to check if at least two of the Maj's (majority ITS2 sequences) are found in >= config.args.coDomSupport samples
                    # # I don't think that there is any need for this extra conservative ruling on whether something is a coDom or not. I think we can just say that if it has two Majs then it is.
                    # coDomSupportedMajDict = {Maj: True for Maj in set(listOfMajsForSamplesWithFootPrint)} # Create a dictionary of all of the MajITS2s found in this type/footprint and evaluate whether they are found in more than 2 samples. If the count of Trues is less than 2 at the end then this is not a codom [under our extra conservative ruling]. If more than one are above the 2 then name accordingly.

                    # #This loop checks each of the Majs that have been identified for the footprint to see how many samples they have been found as Majs in.
                    # #If they have been found in < config.args.coDomSupport number of samples then they are assigned a False in the coDomSupportedMajDict
                    # for Maj in coDomSupportedMajDict.keys():
                    #     if listOfMajsForSamplesWithFootPrint.count(Maj) < config.args.coDomSupport:
                    #         coDomSupportedMajDict[Maj] = False

                    # #Now do the count of False vs. True, if count of True greater than 1 then this is a CoDom and we need to work out the name
                    # if list(coDomSupportedMajDict.values()).count(True) > 1:

                    # This is a supported CoDom now name it:
                    coDom = True
                    newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL', coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE, listofcodommajs = list(set(listOfMajsForSamplesWithFootPrint)))
                    supportedList.append(newSymbiodiniumType.footPrint)
                    # Put the new Type into the samples' collectionlits that have it
                    addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)

                    # else:
                    # # This cannot be called a coDom and we can exit out at this point
                    # # This is an interesting case. So here we have a type that is not considered a coDom using our conservative criteria
                    # # But this type does have multiple MajITS2s (just not enough support for it to be considred a coDom
                    # # This is going to cause trouble later on when we come to write out typeCharacterisation
                    # # As we go through the MajITS2s and this type will therefore come up in two MajITS2 categories.
                    # # We need to find away to make sure that it is only related to
                    # # We are going to solve this problem by adding self.maj to the symbiodiniumType class when it is an initial type
                    # # Then when we come to identify the Majs we can use the inital type majs rather than the clade collection Majs
                    #
                    #     newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL',coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)
                    #     addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
                else: # # This is cannot be called a coDom and we can exit out at this point
                    newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL', coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)
                    addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
                    supportedList.append(newSymbiodiniumType.footPrint)
            else: # This is an unsupportedType
                # newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL',coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)
                unsupportedTypeList[newSymbiodiniumType.footPrint] = listOfSamplesThatContainFootprint


                # Put the new Type into the samples' collectionlits that have it
                # addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
        if len(unsupportedTypeList) > 1:
            searchForFurtherInitialsAgain(unsupportedTypeList = unsupportedTypeList, reqsupport=max(4, math.ceil(config.args.typeSupport*cladeCollectionCountDict[CLADE])))
            a = 5

    return

def assignInitialTypes(cladecollectioncountdict):
    cladeCollectionCountDict = cladecollectioncountdict
    #STEP ONE
    #Make a list of all of the unique footprints
    # Do this clade by clade

    for CLADE in config.args.cladeList:
        footPrintDict  = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for i in range(len(SAMPLE.cladeCollectionList)):
                CLADECOLLECTION = SAMPLE.cladeCollectionList[i]
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.footPrint not in footPrintDict.keys():
                        footPrintDict[CLADECOLLECTION.footPrint] = [[SAMPLE.name], [CLADECOLLECTION.maj]]
                    else:
                        footPrintDict[CLADECOLLECTION.footPrint][0].append(SAMPLE.name)
                        footPrintDict[CLADECOLLECTION.footPrint][1].append(CLADECOLLECTION.maj)
        if CLADE == 'D':
            a = 5
        if len(footPrintDict) > 0:
            collapsedFootPrintDict = searchForFurtherInitialsAgain(unsupportedTypeList = footPrintDict, reqsupport=max(4, math.ceil(config.args.typeSupport*cladeCollectionCountDict[CLADE])))




            for FOOTPRINT in collapsedFootPrintDict.keys():
                # Check that the majs have been called correctly: in one case we had a maj that wasn't part of the footprint due to the collapsing
                if not set(collapsedFootPrintDict[FOOTPRINT][1]).issubset(set(FOOTPRINT)):
                    # Then we have a problem with at least one of the majs that has been called
                    # We need to redo the maj's
                    newMajList = []
                    for SAMPLESIN in collapsedFootPrintDict[FOOTPRINT][0]:
                        for CLADECOLLECTION in config.abundanceList[SAMPLESIN].cladeCollectionList:
                            if CLADECOLLECTION.clade == CLADE:
                                for intraOrdered in CLADECOLLECTION.listOfSeqsAboveCutOff:
                                    if intraOrdered.name in FOOTPRINT:
                                        newMajList.append(intraOrdered.name)
                                        break
                    collapsedFootPrintDict[FOOTPRINT][1] = newMajList
                if len(set(collapsedFootPrintDict[FOOTPRINT][1])) > 1:
                    coDom = True
                    # TODO create the type and assign to the sample(s)
                    newSymbiodiniumType = symbiodiniumType( coDom=coDom, maj=max(set(collapsedFootPrintDict[FOOTPRINT][1]), key=collapsedFootPrintDict[FOOTPRINT][1].count), footPrint=FOOTPRINT, listofSamples=collapsedFootPrintDict[FOOTPRINT][0], clade=CLADE, majList = collapsedFootPrintDict[FOOTPRINT][1])
                    addTypeToSamples(newSymbiodiniumType, collapsedFootPrintDict[FOOTPRINT][0])
                    config.typeDB.addType(newSymbiodiniumType)

                else:
                    coDom = False
                    # TODO create the type and assign to the sample(s)
                    newSymbiodiniumType = symbiodiniumType( coDom=coDom, maj=max(set(collapsedFootPrintDict[FOOTPRINT][1]), key=collapsedFootPrintDict[FOOTPRINT][1].count), footPrint=FOOTPRINT, listofSamples=collapsedFootPrintDict[FOOTPRINT][0], majList = collapsedFootPrintDict[FOOTPRINT][1], clade=CLADE)
                    addTypeToSamples(newSymbiodiniumType, collapsedFootPrintDict[FOOTPRINT][0])
                    config.typeDB.addType(newSymbiodiniumType)

        a = 5



    return

def modelIntraProbDistribution(definingintrainfodict, typename, fithpercentile):
    #ToDo maybe here is where we look for multi nomial distributions

    # orderedListOfIntrasByAvAbundance = [a[0] for a in sorted(config.typeDB[typename].definingIntrasInfo.items(), key=lambda x: x[1][1], reverse=True)]
    orderedListOfIntrasByAvAbundance = [a[0] for a in config.typeDB[typename].sortedDefiningIts2Occurances]

    # First we try the probability distribution approach
    probDistModelSupport = False
    normAbunModelSupport = False
    probDist = [definingintrainfodict[i][1] for i in range(len(orderedListOfIntrasByAvAbundance))]
    successCounter = 0
    # permutations
    perms = 100
    successRate = 0.05
    for i in range(perms):
        arrayOfSampledIntras = np.random.choice(orderedListOfIntrasByAvAbundance, int(fithpercentile), p=probDist)
        if orderedListOfIntrasByAvAbundance[-1] in arrayOfSampledIntras:
            successCounter += 1
    if successCounter/perms >= 1 - successRate:
        probDistModelSupport = True

    # If we consider all intra abundances, in the case where there are big and small abundances then we get a large s.d. and we end up with lots of negative values in the model
    # This may be discounting good intras.
    # What we are worried about is
    # Now we try the normal distribution of the one intra approach
    intraAbundances = definingintrainfodict[-1][0]
    intraAbundancesOutlierRemoved = reject_outliers(definingintrainfodict[-1][0])
    outliersRemoved = [x for x in intraAbundances if x not in intraAbundancesOutlierRemoved]
    if outliersRemoved:
        a = 6
    mu, sigma = statistics.mean(intraAbundancesOutlierRemoved), statistics.stdev(intraAbundancesOutlierRemoved)
    modelelledAbundances = np.random.normal(mu, sigma, 10000)
    proportionOfPosVal = len([x for x in modelelledAbundances if x > 0])/10000
    if proportionOfPosVal >= 1 - successRate:
        normAbunModelSupport = True
    # Currently requring that both of the support methods are True as a sort of conservative check on the intra
    if probDistModelSupport and normAbunModelSupport:
        return True, None
    else:
        return False, orderedListOfIntrasByAvAbundance[-1]

def modelIntraProbDistributionWFinalTypes(definingintrainfodict, typename, fithpercentile):
    #ToDo maybe here is where we look for multi nomial distributions

    # orderedListOfIntrasByAvAbundance = [a[0] for a in sorted(config.typeDB[typename].definingIntrasInfo.items(), key=lambda x: x[1][1], reverse=True)]
    orderedListOfIntrasByAvAbundance = [a[0] for a in config.typeDB[typename].sortedDefiningIts2Occurances]

    # First we try the probability distribution approach
    probDistModelSupport = False
    normAbunModelSupport = False
    probDist = [definingintrainfodict[i][1] for i in range(len(orderedListOfIntrasByAvAbundance))]
    successCounter = 0
    # permutations
    perms = 100
    successRate = 0.05
    for i in range(perms):
        arrayOfSampledIntras = np.random.choice(orderedListOfIntrasByAvAbundance, int(fithpercentile), p=probDist)
        if orderedListOfIntrasByAvAbundance[-1] in arrayOfSampledIntras:
            successCounter += 1
    if successCounter/perms >= 1 - successRate:
        probDistModelSupport = True

    # If we consider all intra abundances, in the case where there are big and small abundances then we get a large s.d. and we end up with lots of negative values in the model
    # This may be discounting good intras.
    # What we are worried about is
    # Now we try the normal distribution of the one intra approach
    intraAbundances = definingintrainfodict[-1][0]
    intraAbundancesOutlierRemoved = reject_outliers(definingintrainfodict[-1][0])
    outliersRemoved = [x for x in intraAbundances if x not in intraAbundancesOutlierRemoved]
    if outliersRemoved:
        a = 6
    mu, sigma = statistics.mean(intraAbundancesOutlierRemoved), statistics.stdev(intraAbundancesOutlierRemoved)
    modelelledAbundances = np.random.normal(mu, sigma, 10000)
    proportionOfPosVal = len([x for x in modelelledAbundances if x > 0])/10000
    if proportionOfPosVal >= 1 - successRate:
        normAbunModelSupport = True
    # Currently requring that both of the support methods are True as a sort of conservative check on the intra
    if probDistModelSupport and normAbunModelSupport:
        return True, None
    else:
        return False, orderedListOfIntrasByAvAbundance[-1]

def reject_outliers(data, m=1.96):
    return [datapoint for datapoint in data if abs(datapoint - np.mean(data)) < m * np.std(data)]


def chDict(value, dict):
    if value in dict.keys():
        return dict[value]
    else:
        return 0



def searchForFurtherInitialsAgain(unsupportedTypeList, reqsupport):


    # I'm going for an end n of 3 here so that at most we try to collapse footprints that are three long into twos
    # I don't think we want to be collapsing footprints containing only two intras
    for n in range(max([len(keyName) for keyName in unsupportedTypeList.keys()]), 2, -1):
        collapseDict = {}
        # This gets a list of the largest footprints that need collapsing. It takes into account that since the
        # start of analysis they may have had other footprints collapsed into them and as such may
        # now be associated with enough samples to put the footprint above the reqsupport. In this case they will
        # be left out of the largestFootprintList and as such excluded from collapsing
        # Here we will
        largestFootprintList = [keyName for keyName in unsupportedTypeList.keys() if len(keyName) == n and len(unsupportedTypeList[keyName][0]) < reqsupport]
        if largestFootprintList:
            for bigFootprint in largestFootprintList:
                for N in range(n,2,-1):
                    nextlargestFootprintList = [keyName for keyName in unsupportedTypeList.keys() if len(keyName)== N - 1]
                    if nextlargestFootprintList:
                        topScore = 0
                        for smallerFootprint in nextlargestFootprintList:
                            if set(smallerFootprint).issubset(set(bigFootprint)):
                                score = len(unsupportedTypeList[bigFootprint][0]) + len(unsupportedTypeList[smallerFootprint][0])
                                if score > topScore:
                                    topScore = score
                                    collapseDict[bigFootprint] = smallerFootprint

                        # Once we have checked the big foot print against all of the next smallest footprint we need to
                        # assess whether we have found a seq for it to collapse into. If we have then great, we have a
                        # note of this and it is onto the next big footprint for us. If not then we need to check the next
                        # biggest footprints etc. etc. until we find one to collapse into. If we don't find one to collapse into
                        # then we remove the footprint from the list and the associated samples will not have an intial supported
                        # types associated with them.
                        if topScore != 0:
                            break
            # Here we have finished going through each of the big footprints
            for footPrintToCollapse in collapseDict.keys():
                unsupportedTypeList[collapseDict[footPrintToCollapse]] = [unsupportedTypeList[collapseDict[footPrintToCollapse]][0] + unsupportedTypeList[footPrintToCollapse][0], unsupportedTypeList[collapseDict[footPrintToCollapse]][1] + unsupportedTypeList[footPrintToCollapse][1]]
                del unsupportedTypeList[footPrintToCollapse]
    '''At this point we have collapsed all but the twos. The twos pose an interesting question.
    We have the opportunity here to assign the single intra Majs as inital types to the remaining samples with
    unidentified types.  If we go through all of the two footprints and take out the intras counting as we go we can
    create an intra count dict. Then if we go back through the two footprints and associate their samples to the intra
    that has the most support in the counts dict. I'm honestly not sure how this will effect the typing over all but
    I think it it probably a fair and sensible thing to do.
    '''

    collectionOfUnsupportedFootprints = [tuples for tuples in unsupportedTypeList.items() if len(tuples[1][0]) < reqsupport and len(tuples[0]) > 1]

    for i in range(len(collectionOfUnsupportedFootprints)):
        # most common intra is equivalent here to high intra
        # The most common maj must be found in the footprint too
        highIntra = max(set(maj for maj in collectionOfUnsupportedFootprints[i][1][1] if maj in collectionOfUnsupportedFootprints[i][0]), key=collectionOfUnsupportedFootprints[i][1][1].count)
        # highIntra = returnHigherDictValue(collectionOfUnsupportedFootprints[i], intraAbundanceDict)
        if frozenset([highIntra]) in unsupportedTypeList.keys():
            unsupportedTypeList[frozenset([highIntra])] = [unsupportedTypeList[frozenset([highIntra])][0] + unsupportedTypeList[collectionOfUnsupportedFootprints[i][0]][0], unsupportedTypeList[frozenset([highIntra])][1] + unsupportedTypeList[collectionOfUnsupportedFootprints[i][0]][1]]
        else:
            unsupportedTypeList[frozenset([highIntra])] = unsupportedTypeList[collectionOfUnsupportedFootprints[i][0]]
        del unsupportedTypeList[collectionOfUnsupportedFootprints[i][0]]
    return unsupportedTypeList

def returnHigherDictValue(setOfThings, dictionary):
    topscore = 0
    for item in list(setOfThings):
        score = dictionary[item]
        if score > topscore:
            topscore = score
    for item in list(setOfThings):
        if dictionary[item] == topscore:
            return item
    return 'ERROR'

def searchForFurtherInitials(unsupportedTypeList, reqsupport, supportedtypelist):
    # convert the frozenset footprints to plain set footprints to prevent any type confusion
    supportedTypeList = [set(a) for a in supportedtypelist]
    foundButAlreadySupported = []
    foundNew = []
    combosAlreadyTested = []
    # Starting with n as the size of the biggest footprint
    # Continue to n of three as we don't want to be looking for single intra footprints
    for n in range(max([len(a[0]) for a in unsupportedTypeList]),2, -1):
        listOfFootprintsWithLenn = [footprintTuple for footprintTuple in unsupportedTypeList if len(footprintTuple[0]) == n]
        for footprint in listOfFootprintsWithLenn:
            for combo in itertools.combinations(footprint[0], n-1):
                # This should help speed up the comparisons
                if set(combo) in combosAlreadyTested:
                    continue
                else:
                    combosAlreadyTested.append(set(combo))
                combinedSupport = footprint[1]
                # Produces a tuple that contains each of the footprints intras
                # compare this subset of intras footprint to all of the other footprints
                for allFootprints in [otherfootprints for otherfootprints in unsupportedTypeList if otherfootprints[0] != footprint[0]]:
                    one = set(combo)
                    two = set(allFootprints[0])
                    if set(combo).issubset(set(allFootprints[0])):
                        # Here we have found a footprint that has additional support
                        # We add together the supports for all such footprints and see if it comes above
                        # the required support for a supported inital type
                        combinedSupport += allFootprints[1]
                if combinedSupport >= reqsupport:
                    # Here we have found what could be a new supported initial type
                    # Check to see if it an already defined initalType
                    # If it is then it means all of the samples involved will find final types of that eventually.
                    if set(combo) in supportedTypeList:
                        foundButAlreadySupported.append(set(combo))
                        #Let's see how this looks so far.
                    else:
                        # Checks to see if the found combo is a subset of any of the other footprint sets that have been found
                        # pieceone = [setfootprints[0] for setfootprints in foundNew]
                        # piecetwo = [set(combo).issubset(checkingfootprints) for checkingfootprints in pieceone]
                        # if not any([set(combo).issubset(checkingfootprints) for checkingfootprints in [setfootprints[0] for setfootprints in foundNew]]):
                        # I have relaxed the require ment of the potential type not to be a subset of any of the other potential types
                        if set(combo) not in [setfootprints[0] for setfootprints in foundNew]:
                            foundNew.append((set(combo), combinedSupport))
    return foundNew

def addTypeToSamples(newSymType, listOfSamplesThatContainFootprint):
    for SAMPLEKEY in listOfSamplesThatContainFootprint:
        SAMPLE = config.abundanceList[SAMPLEKEY]
        if SAMPLE.name in listOfSamplesThatContainFootprint:
             for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                 if CLADECOLLECTION.clade == newSymType.clade: # Then this is a cladeCollection that we want to add the new Symbiodinium Type to
                    # We are going to add the new type to the sample. However we are going to add it as a new instance of the type so that we can have different sortedDefiningITSwOccurance for each sample
                    # if we add exactly the same occurance of the type to multiple samples then when we change one it will change them all.
                    # E.g. if we change the sortedDefining... for one of the samples it will automatically be changed in all of the samples with that type.
                    CLADECOLLECTION.addInitialType(symbiodiniumType(clade=newSymType.clade, coDom=newSymType.coDom, maj=newSymType.maj, footPrint=newSymType.footPrint, listofSamples=newSymType.listOfSamples, majList=newSymType.majList, name=newSymType.name, sorteddefiningits2occurances=newSymType.sortedDefiningIts2Occurances))
                    # Here we add the CLADECOLLECTION.initialType.sortedDefiningIts2Occurances
                    # CLADECOLLECTION.initialType.sortedDefiningIts2Occurances = CLADECOLLECTION.initialType.createSortedDefiningIts2Occurances(SAMPLE.compComplement.listOfits2SequenceOccurances, SAMPLE.totalSeqs)[0]
    return

def createMasterSeqDistancesNonMothur():
    print('Running createMasterSeqDistances()')


    distList = []
    for CLADE in config.args.cladeList:
        listOfNames = []
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for cladeCollection in SAMPLE.cladeCollectionList: # For each cladeCollection
                if cladeCollection.clade == CLADE:
                    for sequenceOccurance in cladeCollection.listOfSeqsAboveCutOff:
                        if sequenceOccurance.name not in listOfNames:
                            listOfNames.append(sequenceOccurance.name)
        # We only need to know the distances of every seq within the clades rather than all seqs total
        # so we will work out the dists clade by clade
        for seq1, seq2 in itertools.combinations(listOfNames,2):
            distList.append('{0} {1} {2}'.format(seq1, seq2, str(config.JSD(config.seqToFFPProbDistDict[seq1], config.seqToFFPProbDistDict[seq2]))))
    masterSeqDistancesDict = {frozenset(a.split(' ')[:-1]): float(a.split(' ')[2]) for a in distList}


    return masterSeqDistancesDict

def assignCladeCollections():
    # This function goes through all of the samples in the abundance list. It checks to see if the sample contains a given number of sequences (above 10% of the total sequences; or the config.args.cladeCollectionCutoff value) of each of the clades
    # If the sample does contain sequences of this clade then it creates a list of the sequences that are above a given proportion of the sequences from that clade (currently also 10% or the config.args.cutoff value)
    # This list forms the basis of the cladeCollection that is appended to the sample within the config.abundanceList.

    # I will also count the number of cladeCollection created for each clade and put them into a dictionary
    # I will use this information to create an intial type support
    cladeCollectionCountDict = {clade: 0 for clade in config.args.cladeList}

    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        for CLADE in config.args.cladeList: # Sequences in the different clades are too divergent to be compared so we have to work on a cladlly separated basis, we are only working with the three main clades, A, C and D
            totalSeqs = SAMPLE.totalSeqs
            cladeSpecificSeqs = sum([a.abundance for a in SAMPLE.compComplement.listOfits2SequenceOccurances if a.clade == CLADE]) # Number of sequence the given sample has that are of the given clade
            if float(cladeSpecificSeqs/totalSeqs) >= config.args.cladeCollectionCutoff:  # If this sample has a proportion of clade X creater than 10% then we will add a cladeCollection to the sample
                cutOffValue = cladeSpecificSeqs * config.args.cutOff # config.args.cutOff is the percentage that will create the cutoff for how abundant a sequence must be in order to be used in the footprint.
                # e.g. 0.1 means that if the sample has 300 clade A sequences, only sequences that are abundant at 30 or greater will be considered in the footprint.
                #tempListOfits2SequenceOccurances = subset of the its2 occurances in the sample that are above the cutOff (so will make up the footprint) and are of the clade in question
                #This list will be in the cladeCollection with normalised abundances out of 1000.
                tempListOfits2SequenceOccurances = [its2SequenceOccurance(name=a.name, abundance=a.abundance, clade=a.clade, sequence=a.sequence) for a in SAMPLE.compComplement.listOfits2SequenceOccurances if a.abundance >= cutOffValue and a.clade == CLADE] # Need to make sure to make a copy of the sequence occurances here so that when we put them into the cladecollection we don't change the abundnaces etc. in the
                tempTotalSeqs = sum([a.abundance for a in tempListOfits2SequenceOccurances]) # Total number of seqs in the its2sequenceoccurances that were above the cutoof for the given clade
                # This loop here nomalises the real sequenced abundances out of 1000 so that samples that had different sequencing depths or different proportions of a given clade do not have more or less sequences.
                i = 0
                while i < len(tempListOfits2SequenceOccurances):
                    tempListOfits2SequenceOccurances[i].abundance = (tempListOfits2SequenceOccurances[i].abundance/tempTotalSeqs)*1000 # Normalise the abundances to 1000
                    i += 1
                # Finally we add the new normalised list to the sample as a cladeCollection
                SAMPLE.addCladeCollection(cladeCollection(CLADE, config.args.cutOff, listofseqsabovecutoff=tempListOfits2SequenceOccurances, foundwithinsample= SAMPLE.name, cladalproportion=cladeSpecificSeqs/totalSeqs))
                cladeCollectionCountDict[CLADE] =  cladeCollectionCountDict[CLADE] + 1

    writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/cladeCollectionCountDict', cladeCollectionCountDict)
    return cladeCollectionCountDict

def inferFinalSymbiodiniumTypesOld():
    print('Running inferFinalSymbiodiniumTypes()')




        # Search for initial type footprints within the complete complement of each sample for a given clade
        # Get list of types that are found in the complete complement first.

        # We need to be careful at this stage, the unsupported types will always fall into the original sample they came from, we need a way of checking whether these insignificant types become supported once we start to use the full complement of samples
        # We also need to be careful because if we simply look for the most complicated type in a sample and we have a sample that has 5 unique intras that are all above the 10% mark but only 1 of them is found in a supported type then this sample will be identified as something very simple which is the one intra that it does have
        # So if we have a sample that still has an intra present above the 10% level (of the cladal sequences) that hasn't been found in one of the types identified in it then we should call this sample unidentified.

        # In order to work out if each of the initial unsupported types becomes supported once the completeComplements are checked and equally whether any of the supported types become unsupported, we need to do 3 phases
        # At the end of this the only Type that may go unidentified by this method is one that is defined according to an intra that is almost always found below the 10% mark
        # i.e. ST will have identifying intras below the 10% mark in at least some samples they will be above the 10% mark so we will have picked them up in our initial type anlysis

        # Phase 1 - Go through all samples pushing in all types, supported and unsupported.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        # Phase 3 - Go through all samples and remove unsupported types then do iter comparisons of types in the sample keeping only the [type with the longest footprint of any two types that share any defining intras]*** actually we need to work through
        # each of these scenarios carefully. If the shorter footprint fits into the longer footprint then yes only keep longer footprint. But if two footprints share intras but one doesn't fit into the other then call it unresolved for the time being.
        # Phase 3 step 1: Check to see if they have any intras in common. Step 2: is ask if one is subset of the other.

        # Probably if they share intras but one doesn't fit in the other then we should keep both types. We would then have to spit the final proportion of the types between the two.
        # The above is impossible to resolve. Consider the three type foot prints [1,2,3], [1,4,5] and [7,4,6] you can see the problem. firt two share the 1 second two share the 4 but can't make a type of all three as first and third do not share any. So lump them all into one list and make a super unresolved Symbiodinium type of them. Not ideal but can't see any other way at current.

        # Phase 4 - Finally if not all intras found at >10% (the defined cutoff) of the cladal sequences in a sample then consider this sample unidentified
        # Also, if we have a sample that still has multiple possible types in it then we need to work out what to do with it.
        # Once all types have been identified we can add the type information (type by type to the finalTypesList) (and, list of ITS2 occurances found only in the types)

        # Keeping all of the inital types, both supported and unsupported creates a real mess at the end with lots of types present.
        # If we end up finding that an inital type supported by say one sample that has three defining intras gets lots of suppot when we feed in all sequences in it
        # The problem is that as the third intra is probably found in most samples at a very low level it is likely that there will be many samples that do belong to this type but don't contain the intra
        # So rather we want to classify more conservaively. So I think that the inital idea of having supported types only considered is a good one. This then plays along a conservative feel
        # this also plays back to the central assumption here that a footprint found many times is less likely to be due to multiple types and rather due to a single given type.
        # I think that it is perhaps a good idea at this point to use a relative cut-off for the support that we require instead of an arbitrary number e.g. 4.
        # I think that we should do say 1% of the number of clade collections we have for a given clade. So if we have 600 samples that have a
        # clade C collection then we should have a minimum type support of 6

        ###### START
        # Get a list of all of the intial supported Symbiodinium types that have unique footprints
    for CLADE in config.args.cladeList:
        footprintList = []
        typeList = []
        typeCountDict = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.initialType.footPrint not in footprintList and CLADECOLLECTION.initialType.supportedType:
                            footprintList.append(CLADECOLLECTION.initialType.footPrint)
                            typeList.append(CLADECOLLECTION.initialType) # We want to keep the actual type class rather than just the frozenset footprint because we want to be able to work out abundances etc. later on


        # Phase 1 - Go through all samples pushing in all supported types.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        # If there are no further types then we don't have a finaltypecladecollectionlist
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name == 'OMd_028' and CLADE == 'D':
                a = 5
            finalTypesList = []
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE: # Then this sample has a set of intras from the given clade that are above the given cladeCollectionCuttoff
                    listOfIntrasInSample = [occurance.name for occurance in SAMPLE.compComplement.listOfits2SequenceOccurances if occurance.clade == CLADE]
                    for TYPE in typeList:
                        if TYPE.footPrint.issubset(listOfIntrasInSample): #  Check to see if the intras in the TYPE.footPrint are found in the listOfIntrasInSample list.
                            # I don't think we need to do this count any more
                            # # Update type count dictionary and add this type to the sample's finalTypesList
                            # if TYPE not in typeCountDict.keys():
                            #     typeCountDict[TYPE] = 1
                            # else:
                            #     typeCountDict[TYPE] = typeCountDict[TYPE] + 1
                            finalTypesList.append(TYPE)
                    # If we have identified types for the final type clade collection list then add the types here
                    # However, if we have not identified a type we will call the Maj its type e.g. sample OMd_028 who's foot print is D4-Otu4554-Otu3721
                    # but this is not supported and there is no D4 on its own supported type so in this case we will add the final type as D4 alone.
                    # beacuse this type does not take into account the Otu4554 or Otu3721 this finaltypecladecollection the FINALTYPECLADECOLLECTION.identified will be False
                    # Or maybe we just leave the final type empty and then don't plot it
                    if len(finalTypesList) > 0:
                        SAMPLE.finalTypeCladeCollectionList.append(finalTypeCladeCollection(foundWithinSample=SAMPLE.name, clade=CLADE, cutoff=config.args.cutOff, listOfFinalTypes=finalTypesList))


        # Phase one and two are complete here. Time for phase three


        # Phase three:  Then check to see if any of the final type footprints are sub sets of each other. Remove any that are (remove the shorter)

        #Go through each sample's finaltypecladecollection for given clade
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == CLADE:




                    # Now iter comparisons and get rid of those types that are subsets of others
                    # FINALTYPECLADECOLLECTION.isMixedIdentification # If we have two final footprints that share intras but one is not a subset of the other then we do not call the type identified and this is True.
                    listOfTypesToGetRidOf = []
                    for a,b in itertools.combinations(FINALTYPECLADECOLLECTION.listOfFinalTypes, 2):
                        if len([intra for intra in a.footPrint if intra in b.footPrint]) > 0: # Are any of a's intras found in b? Then the two footprints share at least one intra i.e. one may include the other
                            if a.footPrint.issubset(b.footPrint): # Then we need to get rid of a
                                if a not in listOfTypesToGetRidOf:
                                    listOfTypesToGetRidOf.append(a)
                            elif b.footPrint.issubset(a.footPrint): # Then we need to get rid of b
                                if b not in listOfTypesToGetRidOf:
                                    listOfTypesToGetRidOf.append(b)
                            else: # Here we have a scenario where the two footprints share intras but one is not a subset of the other
                                FINALTYPECLADECOLLECTION.isMixedIdentification = True
                    for TYPE in listOfTypesToGetRidOf:
                        FINALTYPECLADECOLLECTION.listOfFinalTypes.remove(TYPE)

                    # Once all of the unsupported (both by abundance and subsets) TYPES have been removed check to see if all of the samples > cutoff
                    # (for a given clade) intras are accounted for within the list of types identified

                    # e.g. if an inital type of D1-D7 that was unsupported
                    # end up as a final type of D1, then the D7 component of this symbiont has not been taken into account in the final type call
                    # in this case FINALTYPECLADECOLLECTION.identified would be FALSE
                    for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                        if CLADECOLLECTION.clade == CLADE: # Then this is the clade in questions cladeCollection and we can use the its2 occurances in this as the list of intras found at above the config.args.cutOff
                            # If True then all of the intras above the config.args.cutOff have been accounted for with the final types.
                            # FINALTYPECLADECOLLECTION.typeBasedCompCollection() returns a list of all of the seqs in all of the types in the final type list
                            if set([A.name for A in CLADECOLLECTION.listOfSeqsAboveCutOff]).issubset(FINALTYPECLADECOLLECTION.typeBasedCompCollection()): #and not FINALTYPECLADECOLLECTION.isMixedIdentification: # TODO consider the three type foot prints [1,2,3], [1,4,5] and [7,4,6] you can see the problem
                                FINALTYPECLADECOLLECTION.identified = True
                            else:
                                FINALTYPECLADECOLLECTION.identified = False


                            # Currently the types are written in the listOfFinalTypes as INITIAL types.
                            # So we go through each of the types and we write it out properly as a Final Type to the listOfFinalTypes
                            # I have replaced the typesupport value here with None rather than the typesupportDict that I had been making earlier. I don't see the point of this count.
                            i = 0
                            while i < len(FINALTYPECLADECOLLECTION.listOfFinalTypes):
                                FINALTYPECLADECOLLECTION.listOfFinalTypes[i] = symbiodiniumType(clade=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].clade, footPrint=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].footPrint, maj=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].maj, listofcodommajs=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].listofcodommajs, coDom=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].coDom, typeOfType='FINAL', totalSeqs=SAMPLE.totalSeqs, typeSupport=None, listofoccurences=SAMPLE.compComplement.listOfits2SequenceOccurances, name=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].name)
                                i += 1
                            break # Break out of the for that checks to see if finaltypecladecollection is identified
                    break # Break out of the very first if loop that identifies a finaltypecladecollection of the given clade as there is only one per clade
    print('Completed inferFinalSymbiodiniumTypes()')

def inferFinalSymbiodiniumTypes():
    print('Running inferFinalSymbiodiniumTypes()')
    for CLADE in config.args.cladeList:

        # Phase 1 - Go through all samples pushing in all supported types.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        notAcceptedcount = 0
        accepted = 0
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            finalTypesList = []
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE: # Then this sample has a set of intras from the given clade that are above the given cladeCollectionCuttoff
                    listOfIntrasInSample = set([occurance.name for occurance in SAMPLE.compComplement.listOfits2SequenceOccurances if occurance.clade == CLADE])
                    for TYPE in [config.typeDB[a] for a in config.typeDB.keys() if config.typeDB[a].clade == CLADE]:
                        if TYPE.footPrint.issubset(listOfIntrasInSample): #  Check to see if the intras in the TYPE.footPrint are found in the listOfIntrasInSample list.
                            if SAMPLE.name == 'OMc_081' and TYPE.name =='Otu15163-Otu17432-C1-Otu14865':
                                a = 5
                            if TYPE.name == 'Otu15163/C1':
                                a = 6
                            if len(TYPE.footPrint) > 1:
                                # This is the second/first level of control which makes sure that the ratios of the intras are sensible
                                # This stops us finding types purely diven the presence of the intras but actually checks to the
                                # way in which the intras appear make sense.
                                if ratioAcceptable(SAMPLE, TYPE):
                                    accepted += 1
                                else:
                                    notAcceptedcount += 1
                                    continue
                                    # Check to see if a the type being considered is a subset of the types already assigned or vice versa
                                    # Only keep the largest
                            typeToDel = []
                            isSubSet = False
                            for finaltype in finalTypesList:
                                if set(TYPE.footPrint).issubset(set(finaltype.footPrint)): # Checks to see if new footprint is subset
                                    isSubSet = True
                                if set(finaltype.footPrint).issubset(set(TYPE.footPrint)): # If the current final types are subsets of the new type, delete all such types
                                    typeToDel.append(finaltype)
                            for toDel in typeToDel:
                                finalTypesList.remove(toDel)
                            if isSubSet == False:
                                finalTypesList.append(TYPE)

                    for finaltype in finalTypesList:
                        config.typeDB[finaltype.name].samplesFoundInAsFinal.append(SAMPLE.name)
                    if len(finalTypesList) > 0:
                        SAMPLE.finalTypeCladeCollectionList.append(finalTypeCladeCollection(foundWithinSample=SAMPLE.name, clade=CLADE, cutoff=config.args.cutOff, listOfFinalTypes=[finaltype.name for finaltype in finalTypesList]))
        print('Clade {0}, accepted = {1}, denied = {2}'.format(CLADE, str(accepted), str(notAcceptedcount)))
    print('Completed inferFinalSymbiodiniumTypes()')

def inferFinalSymbiodiniumTypesIterative():
    '''We are going to incorporate iterations into this.
    Once we have been through the samples once we should update the intrainfo for each of the types
    based only on the samples in which the final types have been assigned.
    Then we can go through several iterations each time updating the intrainfo.
    We should keep track of some metric that will allow us to see when no more iterations are useful.
    Maybe we can track how many types were added to samples on each iteration'''

    print('Running inferFinalSymbiodiniumTypes()')
    avNumFinTypes = []
    dictOfTypesToCheckInIterations = {}
    for CLADE in config.args.cladeList:
        typeList = [config.typeDB[a] for a in config.typeDB.keys() if config.typeDB[a].clade == CLADE]
        if typeList:
        # Phase 1 - Go through all samples pushing in all supported types.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3

            notAcceptedcount = 0
            accepted = 0


            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE:  # Then this sample has a set of intras from the given clade that are above the given cladeCollectionCuttoff
                        finalTypesList = []
                        listOfIntrasInSample = set(
                            [occurance.name for occurance in SAMPLE.compComplement.listOfits2SequenceOccurances if
                             occurance.clade == CLADE])
                        for TYPE in typeList:
                            if TYPE.footPrint.issubset(
                                    listOfIntrasInSample):  # Check to see if the intras in the TYPE.footPrint are found in the listOfIntrasInSample list.
                                addToDictList(keyval = '{0}/{1}'.format(SAMPLE.name, CLADE), value = TYPE.name, dictionary = dictOfTypesToCheckInIterations)
                                if SAMPLE.name == 'ADa_039' :
                                    a = 5
                                    if TYPE.name == 'C3/seq178-ST-seq252':
                                        a = 6
                                if len(TYPE.footPrint) > 1:
                                    # This is the second/first level of control which makes sure that the ratios of the intras are sensible
                                    # This stops us finding types purely diven the presence of the intras but actually checks to the
                                    # way in which the intras appear make sense.
                                    if ratioAcceptable(SAMPLE, TYPE, 'INITIAL'):
                                        a = 1
                                    else:
                                        notAcceptedcount += 1
                                        continue
                                        # Check to see if a the type being considered is a subset of the types already assigned or vice versa
                                        # Only keep the largest
                                else:
                                    # Check that intra found at over 0.05 of cladalcollection
                                    if SAMPLE.intraAbundanceDict[list(TYPE.footPrint)[0]]/(SAMPLE.totalSeqs*CLADECOLLECTION.cladalProportion) < 0.1:
                                        notAcceptedcount += 1
                                        continue
                                #Here this is an accepted type.
                                typeToDel = []
                                isSubSet = False
                                for finaltype in finalTypesList:
                                    if set(TYPE.footPrint).issubset(
                                            set(finaltype.footPrint)):  # Checks to see if new footprint is subset
                                        isSubSet = True
                                    if set(finaltype.footPrint).issubset(set(
                                            TYPE.footPrint)):  # If the current final types are subsets of the new type, delete all such types
                                        typeToDel.append(finaltype)
                                for toDel in typeToDel:
                                    finalTypesList.remove(toDel)
                                if isSubSet == False:
                                    finalTypesList.append(TYPE)

                        for finaltype in finalTypesList:
                            config.typeDB[finaltype.name].samplesFoundInAsFinal.append(SAMPLE.name)
                            accepted += 1
                        if len(finalTypesList) > 0:
                            SAMPLE.finalTypeCladeCollectionList.append(
                                finalTypeCladeCollection(foundWithinSample=SAMPLE.name, clade=CLADE,
                                                         cutoff=config.args.cutOff,
                                                         listOfFinalTypes=[finaltype.name for finaltype in
                                                                           finalTypesList]))
            print('Clade {0}, accepted = {1}, denied = {2}'.format(CLADE, str(accepted), str(notAcceptedcount)))

    # Get the average number of final types associated to each finaltypecladecollection
    # Use this as a metric for the degree of additional type assignment
    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
            avNumFinTypes.append(len(FINALTYPECLADECOLLECTION.sortedListOfFinalTypes))
    a = sum(avNumFinTypes)/len(avNumFinTypes)
    print('Av={0}'.format(str(a)))

    multiModalDetection()

    # Here we do the first iteration of trying to get more types in
    typesAdded = 1
    typesAddedDict = {}
    while typesAdded > 0:
        typesAdded = 0
        config.typeDB.generateIntrasInfoFinalForAllTypes()
        for CLADE in config.args.cladeList:
            # typeList = [config.typeDB[a] for a in config.typeDB.keys() if config.typeDB[a].clade == CLADE]
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE:  # Then this sample has a set of intras from the given clade that are above the given cladeCollectionCuttoff
                        setOfIntrasInSample = set([occurance.name for occurance in SAMPLE.compComplement.listOfits2SequenceOccurances if occurance.clade == CLADE])
                        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                            if FINALTYPECLADECOLLECTION.clade == CLADE:
                                # for TYPENAME in [typename for typename in dictOfTypesToCheckInIterations['{0}/{1}'.format(SAMPLE.name, CLADE)] if typename not in FINALTYPECLADECOLLECTION.sortedListOfFinalTypes]:
                                typeNameList = [typename for typename in config.typeDB.keys() if typename not in FINALTYPECLADECOLLECTION.sortedListOfFinalTypes and config.typeDB[typename].clade == CLADE]
                                for TYPENAME in typeNameList:
                                    if SAMPLE.name == 'ADa_039':
                                        a = 'foo'
                                    # if TYPENAME == 'seq162-seq168-C1-seq184':
                                    #     a=5
                                    #     if len(config.typeDB[TYPENAME].samplesFoundInAsFinal) > 0:
                                    #         a = 'fromg'
                                    TYPE = config.typeDB[TYPENAME]
                                    if TYPE.footPrint.issubset(setOfIntrasInSample):  # Check to see if the intras in the TYPE.footPrint are found in the listOfIntrasInSample list.
                                    # If type has lost all support then no longer considered
                                        if len(TYPE.footPrint) > 1 and len(TYPE.samplesFoundInAsFinal) > 0:
                                            if ratioAcceptable(SAMPLE, TYPE, 'FINAL') == False:
                                                continue
                                        elif len(TYPE.footPrint) == 1 and len(TYPE.samplesFoundInAsFinal) > 0:
                                            # Check that single intra footprint/types found at over 0.05 of cladalcollection
                                            try:
                                                if SAMPLE.intraAbundanceDict[list(TYPE.footPrint)[0]] / (SAMPLE.totalSeqs * CLADECOLLECTION.cladalProportion) < 0.1:
                                                    continue
                                            except:
                                                a = 6
                                        else: # If TYPE.samplesFoundInAsFinal not > 0
                                            continue

                                        typeToDel = []
                                        isSubSet = False
                                        try:
                                            for FINALTYPE in [config.typeDB[finaltype] for finaltype in FINALTYPECLADECOLLECTION.sortedListOfFinalTypes]:
                                                if set(TYPE.footPrint).issubset(
                                                        set(FINALTYPE.footPrint)):  # Checks to see if new footprint is subset
                                                    isSubSet = True
                                                if set(FINALTYPE.footPrint).issubset(set(
                                                        TYPE.footPrint)):  # If the current final types are subsets of the new type, delete all such types
                                                    typeToDel.append(FINALTYPE)
                                        except:
                                            a = 5
                                        for toDel in typeToDel:
                                            FINALTYPECLADECOLLECTION.sortedListOfFinalTypes.remove(toDel.name)

                                            config.typeDB[toDel.name].samplesFoundInAsFinal.remove(SAMPLE.name)

                                            # Need to delete this sample from the found in final list for this type entry in the DB
                                        if isSubSet == False:
                                            FINALTYPECLADECOLLECTION.sortedListOfFinalTypes.append(TYPE.name)

                                            config.typeDB[TYPE.name].samplesFoundInAsFinal.append(SAMPLE.name)

                                            # Need to add this sample to the found in final list for this type entry in the DB
                                            typesAdded += 1
                                            if TYPE.name in typesAddedDict.keys():
                                                typesAddedDict[TYPE.name] += 1
                                            else:
                                                typesAddedDict[TYPE.name] = 1
        avNumFinTypes = []
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                avNumFinTypes.append(len(FINALTYPECLADECOLLECTION.sortedListOfFinalTypes))
        a = sum(avNumFinTypes) / len(avNumFinTypes)
        print('Av={0}'.format(str(a)))





    print('Completed inferFinalSymbiodiniumTypes()')

def ratioAcceptable(sample, symtype, initialorfinal):
    '''
    This method will assess the type that is currently about to be put into a sample to see if the abundances in which the type's intras are found make sense.
    For example, if the type is C1/Otu1234 and the average abundances of those two are something like .65 .45 then if we find that the abundances
    of the intras in this sample are like 0.05 and .45 then these intras are likely not due to this type and they will not be included into this sample.
    This function will only take symtypes that have a footprint that contains two or more defining intras.
    I think we will use a set of ratios to define a type. The ratio for each intra will always be the intra in question to the majority intra.
    I.e. we will not keep track of ratio between intras that are not the most majority intra.
    For coDom types we will use the intra that is the maj the most often and if this a tie then we will have to just take the first intra in the name I guess.
    In fact taking the first intra from the name is probably the easiest way of doing it anyway.
    The ratio information can take the form of a dictionary were
    :param sample: A config.abundanceList sample
    :param symtype: A config.typeDB entry
    :return: Bool representing whether the proposed type should be accepted.
    '''
    if initialorfinal == 'FINAL':
        if len(symtype.samplesFoundInAsFinal)<1:
            return False
    sampleCladalProportion = 0
    for cladecollection in sample.cladeCollectionList:
        if cladecollection.clade == symtype.clade:
            sampleCladalProportion = cladecollection.cladalProportion
    # Get the abundances of the intras for the type in the sample

    abundancesOfIntrasInSample = [sample.intraAbundanceDict[intra]/(sample.totalSeqs*sampleCladalProportion) for intra in [a[0] for a in symtype.sortedDefiningIts2Occurances]]

    # Make sure that at least one of the Majs is present above 5%
    # Make sure one of the Majs is the Maj
    majAbun = False
    majMax = 0
    for majs in set(symtype.majList):
        majsAbundance = sample.intraAbundanceDict[majs]/(sample.totalSeqs*sampleCladalProportion)
        if majsAbundance > 0.1:
            majAbun = True
        if majsAbundance > majMax:
            majMax = majsAbundance
    if majAbun == False:
        return False
    # Checks to see that one of the non-maj intras is not in higher abundance than the maj
    if max(abundancesOfIntrasInSample) > majMax:
        return False

    ratiosOfIntrasForTypeInSample = []
    for i in range(len(abundancesOfIntrasInSample)):
        ratiosOfIntrasForTypeInSample.append(abundancesOfIntrasInSample[i]/abundancesOfIntrasInSample[0])

    # Now work through each of the ratios, starting with the second and compare to the ratios of the typeDB entry
    for i in range(1, len(abundancesOfIntrasInSample)):
        # If an intra is bigger than the maj, this should be relected in the ratio
        if ratiosOfIntrasForTypeInSample[i] > 1:
            ratioToCheck = 1 + (1 - (1 / ratiosOfIntrasForTypeInSample[i]))
        else:
            ratioToCheck = ratiosOfIntrasForTypeInSample[i]
        if initialorfinal == 'INITIAL':
            listOfRatios = symtype.definingIntrasInfo[i][3]
        elif initialorfinal == 'FINAL':
            listOfRatios = symtype.intrasInfoFinal[i][3]


        maxR, minR = max(listOfRatios), min(listOfRatios)

        # This is a little tricky but we allow for a 10% increase/decrease on the current deviation from the ratio
        # separately considering above, below
        thresholdValue = 0.1
        maxRThreshold = calcThresh(maxR, thresholdValue)
        if ratioToCheck  < minR/(1+thresholdValue) or ratioToCheck > maxRThreshold:
            return  False
        else:
            # Check to see what absolute read values we are dealing with here
            # If very small then add the type name to a list
            seqsInClade = sample.totalSeqs*sampleCladalProportion
            a = min(abundancesOfIntrasInSample[i-1]*seqsInClade, abundancesOfIntrasInSample[i]*seqsInClade)
            if a < 4:
                someTin = 'G'
    return True

def calcThresh(maxr, thv):
    '''
    This is a somewhat complicated transformation of the ratio data
    to attempt to confom to a linear incease of ratios above the 1 ratio
    Essentially, the further above or below 1 the ratio goes, the smaller the deviance from that ratio is allowed
    '''
    maxR = maxr
    if maxR* (1 + thv) > 1:

        if maxR < 1:
            # Obtain new >1 maxR val
            maxR = maxR*thv
            # Transpose to linear system
            maxR = 1+(1-(1/maxR))
        # Calc new thresh val as a thv% increase on the < 1 ratio equivalent
        return maxR + (thv*(1-(maxR-1)))
    else:
        return maxr * (1 + thv)


def multiModalDetection():
    # Here we identify if there is a bionomial distribution in coDom types between the two most abundant intras
    # If we find a binomial distribution which meets our parameters for clearly being two separate distributions
    # then we separate the type into two new types
    config.typeDB.generateIntrasInfoFinalForAllTypes()
    typesToCreate = []
    typesToDelete = []
    for typekey in config.typeDB.keys():
        TYPE = config.typeDB[typekey]
        # We'll start with just the configs but evetually maybe we check the rest of them too
        if TYPE.coDom:
            if len(TYPE.samplesFoundInAsFinal) > 9:
                #Check only the first two most abundant intras of the coDom
                # for i in range(1, len(TYPE.sortedDefiningIts2Occurances)):
                listOfRatios = TYPE.intrasInfoFinal[1][3]
                x_grid = np.linspace(-2, 4, 2000)
                kde = gaussian_kde(listOfRatios)
                pdf = kde.evaluate(x_grid)


                # returns index of the local max's in pdf
                # Using these index on x_grid will give you x of maxs, use on pdf will give you y
                c = list((np.diff(np.sign(np.diff(pdf))) < 0).nonzero()[0] + 1)
                modes = len(c)

                # plotHists(pdf, x_grid, listOfRatios, TYPE.name)

                # If this appears to be a bimodal distribution

                if modes == 2:
                    # Must be sufficient separation between the peaks in x axis
                    xDiffValid = False
                    if x_grid[c[1]] - x_grid[c[0]] > 0.7:
                        xDiffValid = True
                    # Must also be sufficient diff between minima y and small peak y
                    # This represents the x spread and overlap of the two peaks
                    d = list((np.diff(np.sign(np.diff(pdf))) != 0).nonzero()[0] + 1) # max and min indices
                    if pdf[d[1]]/min([pdf[d[0]], pdf[d[2]]]) > 0.85: # Insufficient separation of peaks
                        xDiffValid = False


                    # Then here we have a candidate for separating into two types
                    # Seperate the samples either side of the minima x axis
                    # Find the minimum x value and then use this ratio as the separator for the two modes
                    # Because the ratio information was put in the same order as the samplesFoundInAsFinal
                    # we can work out which samples are which side of the ratio
                    if xDiffValid:
                        orderedListOfSamplesInType = TYPE.samplesFoundInAsFinal
                        samplesForTypeA = []
                        samplesForTypeB = []
                        # returns the index of local max and mins in pdf
                        # index 1 is the min
                        minX = x_grid[list(((np.diff(np.sign(np.diff(pdf))) != 0).nonzero()[0] + 1))[1]]
                        # for each sample assign ratio to one of the two new types
                        for i in range(len(listOfRatios)):
                            if listOfRatios[i] < minX:
                                samplesForTypeA.append(orderedListOfSamplesInType[i])
                            else:
                                samplesForTypeB.append(orderedListOfSamplesInType[i])
                        # Only create new types if each type is supported by three samples
                        if len(samplesForTypeA) >= 3 and len(samplesForTypeB) >=3:
                            newTypeA = symboidiniumDBTypeEntry(clade = TYPE.clade, samplename = samplesForTypeA, footprint=TYPE.footPrint, binomparentinitsamples=TYPE.samplesFoundInAsInitial)
                            newTypeB = symboidiniumDBTypeEntry(clade = TYPE.clade, samplename = samplesForTypeB, footprint=TYPE.footPrint, binomparentinitsamples=TYPE.samplesFoundInAsInitial)
                            typesToCreate.extend([newTypeA, newTypeB])
                            typesToDelete.append(TYPE.name)
                a = 6
    #TODO if binomials were found then it might be worth passing back through the types incase some of the
    # new types have further splits in them

    # del old types first so that there are minimal unique name conflicts
    for types in typesToDelete:
        #TODO sortedListOfFinalTypes will no longer be sorted. This should be sorted.
        TYPE = config.typeDB[types]
        for SAMPLENAME in TYPE.samplesFoundInAsFinal:
            SAMPLE = config.abundanceList[SAMPLENAME]
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == TYPE.clade:
                    FINALTYPECLADECOLLECTION.sortedListOfFinalTypes.remove(types)
        del config.typeDB[types]


    for types in typesToCreate:
        # Check to make sure we haven't ended up with two types with the same name
        if types.name not in config.typeDB.keys():
            for SAMPLENAME in types.samplesFoundInAsFinal:
                SAMPLE = config.abundanceList[SAMPLENAME]
                for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                    if FINALTYPECLADECOLLECTION.clade == types.clade:
                        FINALTYPECLADECOLLECTION.sortedListOfFinalTypes.append(types.name)
            config.typeDB[types.name] = types
        else:
            for i in range(1,999):
                types.name = types.name + '({0})'.format(str(i))
                if types.name not in config.typesDB.keys():
                    for SAMPLENAME in types.samplesFoundInAsFinal:
                        SAMPLE = config.abundanceList[SAMPLENAME]
                        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                            if FINALTYPECLADECOLLECTION.clade == types.clade:
                                FINALTYPECLADECOLLECTION.sortedListOfFinalTypes.append(types.name)
                    config.typeDB[types.name] = types
                    break


def plotHists(pdf, x_grid, newlist, typename):
    plt.interactive(False)
    fig, ax = plt.subplots(2, sharex=True)
    ax[0].hist(newlist, 100)
    ax[0].set_title(typename)
    ax[0].set_xlim([-2, 4])
    ax[1].plot(x_grid, pdf, color='blue', alpha=0.5, lw=3)
    plt.show()

def addToDictList(keyval, value, dictionary):
    if keyval in dictionary.keys():
        dictionary[keyval].append(value)
    else:
        dictionary[keyval] = [value]
    return


def writeSampleCharacterisationOutput():
    xstr = lambda s: s or ""
    # This is going to be going sample by sample
    outPut = []
    outPut.extend([(['Sample details'], 'title'),([ None, 'Initial Type', '[Initial : Final support]', 'Further types', '[Initial:Final support]'], 'notBold'),(None, 'blankLine')])
    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]

        # Produces cladalPropCounter which is a dict of Key = clade, value = proportion of total seqs
        cladalPropCounter = {}
        for OCCURENCE in SAMPLE.compComplement.listOfits2SequenceOccurances:
            if OCCURENCE.clade in cladalPropCounter.keys():
                cladalPropCounter[OCCURENCE.clade] = cladalPropCounter[OCCURENCE.clade] + OCCURENCE.abundance
            else:
                cladalPropCounter[OCCURENCE.clade] = OCCURENCE.abundance
        for keys in cladalPropCounter.keys():
            cladalPropCounter[keys] = cladalPropCounter[keys]/SAMPLE.totalSeqs

        # Produces a list of tuples in order of most abundant clade
        sortedCladalAbundanceTuple = [(a[0],a[1]) for a in sorted(cladalPropCounter.items(), key=lambda x:x[1], reverse=True)]
        cladalProportionString = ':'.join([a[0] + ' ' + str(format(a[1], '.2f')) for a in sortedCladalAbundanceTuple])
        outPut.append((' / '.join([SAMPLE.name, cladalProportionString, xstr(SAMPLE.hostTaxon), xstr(SAMPLE.region), xstr(SAMPLE.reef)]), 'sample'))


        # Here we have written the header line of a sample. Now time to printout the intial and final types, in order of cladal abundance in the sample
        for CLADE in [a[0] for a in sortedCladalAbundanceTuple]:# Only clades that are represented in the sample in order of abundance
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if len(SAMPLE.finalTypeCladeCollectionList) > 0:
                        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                            if FINALTYPECLADECOLLECTION.clade == CLADE:
                                typeInfo = typeOutputString(CLADECOLLECTION)
                                if len(typeInfo) > 1:
                                    # listOfFinalTypeOutputNames = [typeOutputString(name=finalType.name, sorteddefiningits2occurances=finalType.sortedDefiningIts2Occurances, initialsupportdict=initialSupportDict, finalsupportdict=finalSupportDict) for finalType in orderedListOfFinalTypes]
                                    outPut.append(([None, typeInfo[0].split('\t')[0], typeInfo[0].split('\t')[1], typeInfo[1].split('\t')[0], typeInfo[1].split('\t')[1]], 'notBold'))
                                    # further final type
                                    if len(typeInfo) > 2:
                                        for FFT in typeInfo[2:]:
                                            outPut.append(([None, None, None, FFT.split('\t')[0], FFT.split('\t')[1]], 'notBold'))
                                else:# If there aren't any further types due to all seqs in the finaltypecladecollection being removed
                                    outPut.append(([None, typeInfo[0].split('\t')[0], typeInfo[0].split('\t')[1], None, None], 'notBold'))


        outPut.append((None, 'blankLine'))

    jinjaEnvironment = Environment(loader=FileSystemLoader(config.args.rootLocation + r'/html templates'), trim_blocks=True)
    jinjaTemplate = jinjaEnvironment.get_template('sampleCharacterisation__TEMPLATE.html')
    stringThing = jinjaTemplate.render(outPut=outPut)
    htmlString = [a for a in stringThing.split('\n')]

    return htmlString

def typeOutputString(cladecollection):
    ''' Return the inital and any final type names and the abundance of their intras'''
    SAMPLE = config.abundanceList[cladecollection.foundWithinSample]
    nameInfoToReturn = []
    intraAbundanceDict = config.abundanceList[cladecollection.foundWithinSample].intraAbundanceDict
    initialName = cladecollection.initialType.name
    abundanceInfo = '[{0}]'.format(':'.join(["{0:.2f}".format(((intraAbundanceDict[intraName[0]])/(SAMPLE.totalSeqs*cladecollection.cladalProportion))) for intraName in cladecollection.initialType.sortedDefiningIts2Occurances]))
    nameInfoToReturn.append('\t'.join([initialName,abundanceInfo]))

    for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
        if FINALTYPECLADECOLLECTION.clade == cladecollection.clade:
            for FINALTYPE in FINALTYPECLADECOLLECTION.sortedListOfFinalTypes:
                abundanceInfo = '[{0}]'.format(':'.join(["{0:.2f}".format(((intraAbundanceDict[intraName[0]])/(SAMPLE.totalSeqs*cladecollection.cladalProportion))) for intraName in config.typeDB[FINALTYPE].sortedDefiningIts2Occurances]))
                nameInfoToReturn.append('\t'.join([FINALTYPE, abundanceInfo]))

    return nameInfoToReturn

def typeOutputStringOld(name, sorteddefiningits2occurances, initialsupportdict, finalsupportdict):

    # Convert the sortedDefiningITS2Occurances into a dictionary that we can use to add the proportions with later on.
    # By converting the sequences to laJs where possible we can directly use the name to call up the abundances of the intras from this dict
    abundanceDict = {}
    for tuple in sorteddefiningits2occurances:
        # if tuple[0] in config.oursToLaJDict.keys():
        #     abundanceDict[CLJ(tuple[0])] = tuple[1]
        # else:
        #     abundanceDict[tuple[0]] = tuple[1]
        abundanceDict[CLJ(tuple[0])] = tuple[1]
    # Deconstruct the name and get a list of dashes and slashes
    dashAndSlashList = []
    for character in name:
        if character == '/' or character == '-':
            dashAndSlashList.append(character)

    # Get each intra in the name allowing for both '/'s and '-'s.
    listOfIntrasInName = []
    tempList = name.split('-')
    for item in tempList:
        if '/' in item:
            listOfIntrasInName.extend(item.split('/'))
        else:
            listOfIntrasInName.append(item)

    # Create the string for the type e.g. C3[0.65]-Otu12987[0.17]
    newString = []
    i = 0
    while i < len(listOfIntrasInName):
        # Then this is the first and not last Intra in the name and there needs the string initalized and a dash after
        if i == 0 and len(listOfIntrasInName) != 1:
            newString = listOfIntrasInName[i] + '[' + str(format(abundanceDict[listOfIntrasInName[i]], '.2f')) + ']' + dashAndSlashList[i]
        # Then this is the first and last Intra in the name and therefore needs the string initalising but not dash after
        elif i == 0 and len(listOfIntrasInName) == 1:
            newString = listOfIntrasInName[i] + '[' + str(format(abundanceDict[listOfIntrasInName[i]], '.2f')) + ']'
        # This is the last intra in the list but not the first. The string doesn't need initializing and no dash.
        elif i == len(listOfIntrasInName) - 1:
            newString = newString + listOfIntrasInName[i] + '[' + str(format(abundanceDict[listOfIntrasInName[i]], '.2f')) + ']'
        # This is not the first and not the last intra so doesn't require initialising but does require a slash or dash
        else:
            newString = newString + listOfIntrasInName[i] + '[' + str(format(abundanceDict[listOfIntrasInName[i]], '.2f')) + ']' + dashAndSlashList[i]
        i += 1

    # Now add the inital and final support information
    # this should be tab seperable from the type string
    # If the name is not found in the finalSupportDict then this means that there were 0 final types for that type
    if name in finalsupportdict.keys():
        newString = newString + '\t' + '[' + str(initialsupportdict[name]) + ':' + str(finalsupportdict[name]) + ']'
    else:
        newString = newString + '\t' + '[' + str(initialsupportdict[name]) + ':0]'

    return newString

def createHtmlHolder():
    tempList = []
    tempList.extend(['<!DOCTYPE html>', '<html>', '<head>', '<title>SymbiodiniumType Output</title>', '</head>', '<body>'])
    return tempList

def closeAndWriteHtmlHolder(htmlfile):
    htmlfile.extend(['</body>', '</html>'])
    writeListToDestination(config.args.saveLocation + '/html outputs/symbiodiniumTyping_OUTPUT_' + time.strftime("%H%M%S") + '_' + time.strftime("%d%m%y") +'.html', htmlfile)

def MDS(listmatrix):
    # This takes a list of lists that represent pairwise distances between a list of objects and returns a set of points from a 2d MDS
    # We will send in a listmatrix with column and row headers
    # First of all we can get the columnRowNames variable from this matrix
    # We can then delete the first row as these are effectively column titles we can then remove the first
    # item of each of the lists as these will be effectively the row names
    # Distance file available from RMDS project:
    #    https://github.com/cheind/rmds/blob/master/examples/european_city_distances.csv
    #For an example you can use the below
    # reader = csv.reader(open(r"C:\Users\User\Google Drive\Voolstra OTU paper\testDists.csv", "r"), delimiter=';')
    # data = list(reader)

    # Where listMatrix is a list of lists. With each list having the distances as floats
    # of that given item to all the other items in the same order
    listMatrix = listmatrix

    # Where columnRowNames is the list of the items in the correct order
    columnRowNames =  listMatrix[0][1:] # This ignores the first item which will be 'types'

    # This removes the first list in the matrix which was the column headings
    del listMatrix[0]

    # This gets rid of the row headers
    # At this point we are left with purely the numbers of the distance matrix
    i = 0
    while i < len(listMatrix):
        del listMatrix[i][0]
        i += 1

    # Convert to numpy array
    adist = np.array(listMatrix)

    # Normalise to the largest distance # This may not be necessary as we are already working with small numbers rather than the euclidean distance that the example used
    amax = np.amax(adist)
    adist /= amax

    # Configure the MDS
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    # Run the MDS and get results
    results = mds.fit(adist)

    # The coords come out in the order of the matrix input, in the example this is the order of the cities
    # Coords can be accessed as e.g. coords[0] will giev the first set of coordinates
    coords = results.embedding_

    return coords

def createMatrixFromColDistsJSD(listoffinaltypes, JSDcoldistdict):
    # Create empty matrix that is of dimensions [rows(n samples + 1)][columns(number of samples + 1)]
    Matrix = [list(listoffinaltypes)]
    Matrix[0].insert(0, 'Types')
    for sample in listoffinaltypes:
        t = ['N/A'] * len(listoffinaltypes)
        t.insert(0, sample)
        Matrix.append(t)

    # Matrix is FstMatrix[row][column]
    # Fill matrix from the column Fsts in FstColDist
    row = 1
    while row < len(Matrix):  # For each row
        col = 1
        while col < len(Matrix[row]):  # For each col
            typeOne = Matrix[0][col]
            typeTwo = Matrix[row][0]
            if typeOne == typeTwo:
                Matrix[row][col] = 0.00
            else:
                Matrix[row][col] = Matrix[col][row] = max(JSDcoldistdict[frozenset({typeOne, typeTwo})], 0.001)
            col += 1
        row += 1

    return Matrix

def investigateFinalTypeCorrelations():
    # If we do this all within clade then it will mean that we end up with far fewer pariwise comparisons which will make the
    # computation go a lot faster.
    # First get list of supported inital types, this will essentially be the same as the list of final types.
    # the list of the types will give us the order of the probability distribution
    # Then initiate a list that contains a list for each type. each one of these lists should contain the number of sample number of 0s
    # THen go sample by sample through each of the final clade collections and each of the final clde types.
    # For every final clade type that is found, add a one to the appropriate list
    # YOu can find the index of which list this is by looking for the name of that type in the list of type that you made at the begging.
    # Once you have all of the lists populated i.e. you have been through each of the samples, then you can normalise the lists so as to make them prob distributions
    # You can then do pairwise comparisons of all of the types putting their distributions into JSD.
    # From here we can make a pairwise distance matrix, feed this into the MDS to retrieve coordinates
    # At this point we should probably visualise this so that we can make a decision on how we are going to call which pairs of types to collapse into one.
    # One option is to calculate the mean distance of the group and then have a cut off for if a pair of sequence are 'particularly small'.
    # We are hoping that when we visualise the initial MDS that that there are some clear gropings of the distances.

    for CLADE in config.args.cladeList:
        cladeSpecificListOfFinalTypes = []
        cladeSpecificListOfSamples = []
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == CLADE:
                    if len(FINALTYPECLADECOLLECTION.listOfFinalTypes) > 0:
                        cladeSpecificListOfSamples.append(SAMPLE.name)
                    for FINALTYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
                        if FINALTYPE.name not in cladeSpecificListOfFinalTypes:
                            cladeSpecificListOfFinalTypes.append(FINALTYPE.name)
        # Here we have a list of all of the final types that are found in all of the samples for the given clade
        #There is no point in doing a correlation analysis if there are one or less types

        if len(cladeSpecificListOfFinalTypes) > 1:
        # Initialize a list of tuples in which there is the typeName as the first item and the second item is a list
        #  which is the typesdistribution across the samples from the clade in question which is initially filled with 0s,
        # one for ever sample
            listOfTypeDistributions = [(typeName, [0 for a in range(len(cladeSpecificListOfSamples))]) for typeName in
                                       cladeSpecificListOfFinalTypes]

            # Now we revisit each sample and populate the type lists according to the abundance of types found in each sample
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                if SAMPLE.name in cladeSpecificListOfSamples:
                    for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                        if FINALTYPECLADECOLLECTION.clade == CLADE:
                            for FINALTYPE in FINALTYPECLADECOLLECTION.listOfFinalTypes:
                            # When we get here then we have found a Final type that we need to count in the type distributions
                            # We change the 0 count to a 1
                                listOfTypeDistributions[cladeSpecificListOfFinalTypes.index(FINALTYPE.name)][1][
                                    cladeSpecificListOfSamples.index(SAMPLE.name)] = 1
            # Here we have the listOfTypeDistributions populated with the counts
            # We will now compare distributions of types
            # We will compare type pairs. If they are found in very similar samples then they are candidates for collapse
            # To do the comparison we will only look at samples in which the type with the lower abundance is found
            # Explanation below:
            '''
            I think the concept of this works well in general but I think that we should only count the distribution
            of the less frequent type. This way, if the less frequent type is found in every sample withthe more frequent
            we will get pure correlation. If we compared all of the samples in which either of the types were found
            and one type is found in far more samples, then we end up with effectively lots of samples that the less frequenct
            type is not found and so that looks like a positive lack of correlation when in fact it may be due to simply
            not having more data types for that type. If the less frequent type is indded found in samples without the more frequent
            type then that is fair game to say that that is evidnce for a lack of correlation and them being separate types.
            '''
            # We will use a cutoff of 0.8 i.e. 80% similarity to collapse

            collapseList = []
            for finalType1, finalType2 in itertools.combinations(listOfTypeDistributions, 2):

                correlationCoefficient = 0
                # Identify which of the final types is more abundant i.e. found in more samples
                if sum(finalType1[1]) > sum(finalType2[1]):
                    # In this case then we want to only be concerned with samples in which finalType2 are found
                    newDistrOne = [finalType1[1][i] for i in range(len(finalType2[1])) if finalType2[1][i] == 1]
                    newDistrTwo = [finalType2[1][i] for i in range(len(finalType2[1])) if finalType2[1][i] == 1]
                    # Work out the correlation as a ratio of the number of samples both types are found in relative
                    # to the number of samples the least abundant type is fond in
                    correlationCoefficient = sum(newDistrOne)/sum(newDistrTwo)

                else:
                    # In this case then we want to only be concerned with samples in which finalType1 are found
                    newDistrOne = [finalType1[1][i] for i in range(len(finalType1[1])) if finalType1[1][i] == 1]
                    newDistrTwo = [finalType2[1][i] for i in range(len(finalType1[1])) if finalType1[1][i] == 1]
                    correlationCoefficient = sum(newDistrTwo) / sum(newDistrOne)

                # Add the pairs that will be collapsed to list
                if correlationCoefficient > 0.8:
                    collapseList.append(({finalType1[0], finalType2[0]}, correlationCoefficient))



            a = 5

                # Here we need to consider what will happen if we have triangles (or more) of types that need collapse
                # i.e. a collapses to b and b collapses to c


            # We should then combine the two outliers into a single name
            # Then go through the samples and look for either of the pairs and replace them with the single new type
            # We will then need to do another count of the final types to look for support
            # We will also have to possibly replace the inital types with this new combined type else the initial types will
            # Look as though they have lost all support as you will never find it as a final type

            # I am jumping the gun here a little bit. We still need to write the code that collapses the types that have been found to be correlated above.
            '''
            We have another situation that we need to work on though.
            When we have a sample that contains two final types that share at least one sequence but one is not the subset of the other and they are not correlated
            then we need to come up with a way of calling which type it is.
            I think a good way to do this is with discriminant anlysis
            The code here is perfect I hope: Where X is the actaul variables, each string here represents a sample
            with each item within the string represnting the value for each variable
            y in this case represents the classification of the samples which in out case will be final t ypes:

            import numpy as np
            >>> from sklearn.lda import LDA
            >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
            >>> y = np.array([1, 1, 1, 2, 2, 2])
            >>> clf = LDA()
            >>> clf.fit(X, y)
            LDA(n_components=None, priors=None, shrinkage=None, solver='svd',
              store_covariance=False, tol=0.0001)
            >>> print(clf.predict([[-0.8, -1]]))

            ther
            '''

def distBetweenTwoPoints(coord1, coord2):
    return math.hypot(coord2[0]-coord1[0], coord2[1]-coord1[1])

def computeCutOffPerformanceParameters(cutoff):

    cladeCollectionCountDict = assignCladeCollectionsPermute(cutoff)
    assignInitialTypesPermute(cladeCollectionCountDict)
    inferFinalSymbiodiniumTypesPermute(cutoff)
    supportedInitialTypes, numMajs, finalTypesIdentified, unresolvedFinalTypes = calculatePerformaceParametersPermute()
    return supportedInitialTypes, numMajs, finalTypesIdentified, unresolvedFinalTypes

def assignCladeCollectionsPermute(cutoff):
    # This function goes through all of the samples in the abundance list. It checks to see if the sample contains a given number of sequences (above 10% of the total sequences; or the config.args.cladeCollectionCutoff value) of each of the clades
    # If the sample does contain sequences of this clade then it creates a list of the sequences that are above a given proportion of the sequences from that clade (currently also 10% or the config.args.cutoff value)
    # This list forms the basis of the cladeCollection that is appended to the sample within the config.abundanceList.

    # I will also count the number of cladeCollection created for each clade and put them into a dictionary
    # I will use this information to create an intial type support
    cladeCollectionCountDict = {clade: 0 for clade in config.args.cladeList}

    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        for CLADE in config.args.cladeList: # Sequences in the different clades are too divergent to be compared so we have to work on a cladlly separated basis, we are only working with the three main clades, A, C and D
            totalSeqs = SAMPLE.totalSeqs
            cladeSpecificSeqs = sum([a.abundance for a in SAMPLE.compComplement.listOfits2SequenceOccurances if a.clade == CLADE]) # Number of sequence the given sample has that are of the given clade
            if float(cladeSpecificSeqs/totalSeqs) >= config.args.cladeCollectionCutoff:  # If this sample has a proportion of clade X creater than 10% then we will add a cladeCollection to the sample
                cutOffValue = cladeSpecificSeqs * cutoff # config.args.cutOff is the percentage that will create the cutoff for how abundant a sequence must be in order to be used in the footprint.
                # e.g. 0.1 means that if the sample has 300 clade A sequences, only sequences that are abundant at 30 or greater will be considered in the footprint.
                #tempListOfits2SequenceOccurances = subset of the its2 occurances in the sample that are above the cutOff (so will make up the footprint) and are of the clade in question
                #This list will be in the cladeCollection with normalised abundances out of 1000.
                tempListOfits2SequenceOccurances = [its2SequenceOccurance(name=a.name, abundance=a.abundance, clade=a.clade, sequence=a.sequence) for a in SAMPLE.compComplement.listOfits2SequenceOccurances if a.abundance >= cutOffValue and a.clade == CLADE] # Need to make sure to make a copy of the sequence occurances here so that when we put them into the cladecollection we don't change the abundnaces etc. in the
                tempTotalSeqs = sum([a.abundance for a in tempListOfits2SequenceOccurances]) # Total number of seqs in the its2sequenceoccurances that were above the cutoof for the given clade
                # This loop here nomalises the real sequenced abundances out of 1000 so that samples that had different sequencing depths or different proportions of a given clade do not have more or less sequences.
                i = 0
                while i < len(tempListOfits2SequenceOccurances):
                    tempListOfits2SequenceOccurances[i].abundance = (tempListOfits2SequenceOccurances[i].abundance/tempTotalSeqs)*1000 # Normalise the abundances to 1000
                    i += 1
                # Finally we add the new normalised list to the sample as a cladeCollection
                SAMPLE.addCladeCollection(cladeCollection(CLADE, cutoff, listofseqsabovecutoff=tempListOfits2SequenceOccurances, foundwithinsample= SAMPLE.name, cladalproportion=cladeSpecificSeqs/totalSeqs))
                cladeCollectionCountDict[CLADE] =  cladeCollectionCountDict[CLADE] + 1

    writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/cladeCollectionCountDict', cladeCollectionCountDict)
    return cladeCollectionCountDict

def assignInitialTypesPermute(cladecollectioncountdict):
    cladeCollectionCountDict = cladecollectioncountdict
    #STEP ONE
    #Make a list of all of the unique footprints
    # Do this clade by clade
    for CLADE in config.args.cladeList:
        listOfFootprints  = []
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name == 'ADa_011' or SAMPLE.name == 'OMd_028':
                a = 1
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.footPrint not in listOfFootprints:
                        listOfFootprints.append(CLADECOLLECTION.footPrint)
        # STEP TWO
        # For each footprint, work out whether it is a coDom and whether it is a supported type
        # Reassign to the samples' cladeCollections that contain it.

        for FOOTPRINT in listOfFootprints:
            coDom = False
            supportedType = False
            listOfSamplesThatContainFootprint = [] # List of the sample names that have a cladeCollection that match the footprint (exactly)
            listOfMajsForSamplesWithFootPrint = [] # List of the sequence names that are found as the predominant seqs in the footprint in question
            for SAMPLEKEY in config.abundanceList.keys():
                SAMPLE = config.abundanceList[SAMPLEKEY]
                for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                    if CLADECOLLECTION.clade == CLADE and CLADECOLLECTION.footPrint == FOOTPRINT:
                        listOfSamplesThatContainFootprint.append(SAMPLE.name)
                        listOfMajsForSamplesWithFootPrint.append(CLADECOLLECTION.maj)

            # Ascertain whether this is a supported type and if so check to see if it is a coDom
            # There are enough samples with this footprint to define it as a type, continue checking it for being a coDom or if there is only one seq in the footprint then this is a Maj and therefore must be a type.
            if len(listOfSamplesThatContainFootprint) >= max(4, math.ceil(config.args.typeSupport*cladeCollectionCountDict[CLADE])) or len(FOOTPRINT) == 1:
                supportedType = True
                if len(set(listOfMajsForSamplesWithFootPrint)) > 1: # More than one Maj ITS2 seq: This is a coDom strictly speaking now we need to check if at least two of the Maj's (majority ITS2 sequences) are found in >= config.args.coDomSupport samples
                    # # I don't think that there is any need for this extra conservative ruling on whether something is a coDom or not. I think we can just say that if it has two Majs then it is.
                    # coDomSupportedMajDict = {Maj: True for Maj in set(listOfMajsForSamplesWithFootPrint)} # Create a dictionary of all of the MajITS2s found in this type/footprint and evaluate whether they are found in more than 2 samples. If the count of Trues is less than 2 at the end then this is not a codom [under our extra conservative ruling]. If more than one are above the 2 then name accordingly.

                    # #This loop checks each of the Majs that have been identified for the footprint to see how many samples they have been found as Majs in.
                    # #If they have been found in < config.args.coDomSupport number of samples then they are assigned a False in the coDomSupportedMajDict
                    # for Maj in coDomSupportedMajDict.keys():
                    #     if listOfMajsForSamplesWithFootPrint.count(Maj) < config.args.coDomSupport:
                    #         coDomSupportedMajDict[Maj] = False

                    # #Now do the count of False vs. True, if count of True greater than 1 then this is a CoDom and we need to work out the name
                    # if list(coDomSupportedMajDict.values()).count(True) > 1:

                    # This is a supported CoDom now name it:
                    coDom = True
                    newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL', coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE, listofcodommajs = list(set(listOfMajsForSamplesWithFootPrint)))

                    # Put the new Type into the samples' collectionlits that have it
                    addTypeToSamplesPermute(newSymbiodiniumType, listOfSamplesThatContainFootprint)

                    # else:
                    # # This cannot be called a coDom and we can exit out at this point
                    # # This is an interesting case. So here we have a type that is not considered a coDom using our conservative criteria
                    # # But this type does have multiple MajITS2s (just not enough support for it to be considred a coDom
                    # # This is going to cause trouble later on when we come to write out typeCharacterisation
                    # # As we go through the MajITS2s and this type will therefore come up in two MajITS2 categories.
                    # # We need to find away to make sure that it is only related to
                    # # We are going to solve this problem by adding self.maj to the symbiodiniumType class when it is an initial type
                    # # Then when we come to identify the Majs we can use the inital type majs rather than the clade collection Majs
                    #
                    #     newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL',coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), supportedType = supportedType, footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)
                    #     addTypeToSamples(newSymbiodiniumType, listOfSamplesThatContainFootprint)
                else: # # This is cannot be called a coDom and we can exit out at this point
                    newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL', coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)
                    addTypeToSamplesPermute(newSymbiodiniumType, listOfSamplesThatContainFootprint)
            else: # This is an unsupportedType
                newSymbiodiniumType = symbiodiniumType(typeOfType='INITIAL',coDom = coDom, maj = max(set(listOfMajsForSamplesWithFootPrint), key=listOfMajsForSamplesWithFootPrint.count), footPrint = FOOTPRINT, listofSamples=listOfSamplesThatContainFootprint, clade=CLADE)

                # TODO at this point there is no way of knowing whether this type is coDom or not as it is in so few samples. We may have to work this out later on.
                # Put the new Type into the samples' collectionlits that have it
                addTypeToSamplesPermute(newSymbiodiniumType, listOfSamplesThatContainFootprint)


    return

def inferFinalSymbiodiniumTypesPermute(cutoff):
    print('Running inferFinalSymbiodiniumTypes()')




        # Search for initial type footprints within the complete complement of each sample for a given clade
        # Get list of types that are found in the complete complement first.

        # We need to be careful at this stage, the unsupported types will always fall into the original sample they came from, we need a way of checking whether these insignificant types become supported once we start to use the full complement of samples
        # We also need to be careful because if we simply look for the most complicated type in a sample and we have a sample that has 5 unique intras that are all above the 10% mark but only 1 of them is found in a supported type then this sample will be identified as something very simple which is the one intra that it does have
        # So if we have a sample that still has an intra present above the 10% level (of the cladal sequences) that hasn't been found in one of the types identified in it then we should call this sample unidentified.

        # In order to work out if each of the initial unsupported types becomes supported once the completeComplements are checked and equally whether any of the supported types become unsupported, we need to do 3 phases
        # At the end of this the only Type that may go unidentified by this method is one that is defined according to an intra that is almost always found below the 10% mark
        # i.e. ST will have identifying intras below the 10% mark in at least some samples they will be above the 10% mark so we will have picked them up in our initial type anlysis

        # Phase 1 - Go through all samples pushing in all types, supported and unsupported.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        # Phase 3 - Go through all samples and remove unsupported types then do iter comparisons of types in the sample keeping only the [type with the longest footprint of any two types that share any defining intras]*** actually we need to work through
        # each of these scenarios carefully. If the shorter footprint fits into the longer footprint then yes only keep longer footprint. But if two footprints share intras but one doesn't fit into the other then call it unresolved for the time being.
        # Phase 3 step 1: Check to see if they have any intras in common. Step 2: is ask if one is subset of the other.

        # Probably if they share intras but one doesn't fit in the other then we should keep both types. We would then have to spit the final proportion of the types between the two.
        # The above is impossible to resolve. Consider the three type foot prints [1,2,3], [1,4,5] and [7,4,6] you can see the problem. firt two share the 1 second two share the 4 but can't make a type of all three as first and third do not share any. So lump them all into one list and make a super unresolved Symbiodinium type of them. Not ideal but can't see any other way at current.

        # Phase 4 - Finally if not all intras found at >10% (the defined cutoff) of the cladal sequences in a sample then consider this sample unidentified
        # Also, if we have a sample that still has multiple possible types in it then we need to work out what to do with it.
        # Once all types have been identified we can add the type information (type by type to the finalTypesList) (and, list of ITS2 occurances found only in the types)

        # Keeping all of the inital types, both supported and unsupported creates a real mess at the end with lots of types present.
        # If we end up finding that an inital type supported by say one sample that has three defining intras gets lots of suppot when we feed in all sequences in it
        # The problem is that as the third intra is probably found in most samples at a very low level it is likely that there will be many samples that do belong to this type but don't contain the intra
        # So rather we want to classify more conservaively. So I think that the inital idea of having supported types only considered is a good one. This then plays along a conservative feel
        # this also plays back to the central assumption here that a footprint found many times is less likely to be due to multiple types and rather due to a single given type.
        # I think that it is perhaps a good idea at this point to use a relative cut-off for the support that we require instead of an arbitrary number e.g. 4.
        # I think that we should do say 1% of the number of clade collections we have for a given clade. So if we have 600 samples that have a
        # clade C collection then we should have a minimum type support of 6

        ###### START
        # Get a list of all of the intial supported Symbiodinium types that have unique footprints
    for CLADE in config.args.cladeList:
        footprintList = []
        typeList = []
        typeCountDict = {}
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE:
                    if CLADECOLLECTION.initialType.footPrint not in footprintList and CLADECOLLECTION.initialType.supportedType:
                            footprintList.append(CLADECOLLECTION.initialType.footPrint)
                            typeList.append(CLADECOLLECTION.initialType) # We want to keep the actual type class rather than just the frozenset footprint because we want to be able to work out abundances etc. later on


        # Phase 1 - Go through all samples pushing in all supported types.
        # Phase 2 - Go through all samples doing a count of types (make list) and identify supported and unsupported: This will tell us which of the previously unsupported initial types are now supported. It will not tell us which of the initial types that may now not be supported but these will be dropped out in phase 3
        # If there are no further types then we don't have a finaltypecladecollectionlist
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            if SAMPLE.name == 'OMd_028' and CLADE == 'D':
                a = 5
            finalTypesList = []
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == CLADE: # Then this sample has a set of intras from the given clade that are above the given cladeCollectionCuttoff
                    listOfIntrasInSample = [occurance.name for occurance in SAMPLE.compComplement.listOfits2SequenceOccurances if occurance.clade == CLADE]
                    for TYPE in typeList:
                        if TYPE.footPrint.issubset(listOfIntrasInSample): #  Check to see if the intras in the TYPE.footPrint are found in the listOfIntrasInSample list.
                            # I don't think we need to do this count any more
                            # # Update type count dictionary and add this type to the sample's finalTypesList
                            # if TYPE not in typeCountDict.keys():
                            #     typeCountDict[TYPE] = 1
                            # else:
                            #     typeCountDict[TYPE] = typeCountDict[TYPE] + 1
                            finalTypesList.append(TYPE)
                    # If we have identified types for the final type clade collection list then add the types here
                    # However, if we have not identified a type we will call the Maj its type e.g. sample OMd_028 who's foot print is D4-Otu4554-Otu3721
                    # but this is not supported and there is no D4 on its own supported type so in this case we will add the final type as D4 alone.
                    # beacuse this type does not take into account the Otu4554 or Otu3721 this finaltypecladecollection the FINALTYPECLADECOLLECTION.identified will be False
                    # Or maybe we just leave the final type empty and then don't plot it
                    if len(finalTypesList) > 0:
                        SAMPLE.finalTypeCladeCollectionList.append(finalTypeCladeCollection(foundWithinSample=SAMPLE.name, clade=CLADE, cutoff=cutoff, listOfFinalTypes=finalTypesList))


        # Phase one and two are complete here. Time for phase three


        # Phase three:  Then check to see if any of the final type footprints are sub sets of each other. Remove any that are (remove the shorter)

        #Go through each sample's finaltypecladecollection for given clade
        for SAMPLEKEY in config.abundanceList.keys():
            SAMPLE = config.abundanceList[SAMPLEKEY]
            for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
                if FINALTYPECLADECOLLECTION.clade == CLADE:
                    # Now iter comparisons and get rid of those types that are subsets of others
                    # FINALTYPECLADECOLLECTION.isMixedIdentification # If we have two final footprints that share intras but one is not a subset of the other then we do not call the type identified and this is True.
                    listOfTypesToGetRidOf = []
                    for a,b in itertools.combinations(FINALTYPECLADECOLLECTION.listOfFinalTypes, 2):
                        if len([intra for intra in a.footPrint if intra in b.footPrint]) > 0: # Are any of a's intras found in b? Then the two footprints share at least one intra i.e. one may include the other
                            if a.footPrint.issubset(b.footPrint): # Then we need to get rid of a
                                if a not in listOfTypesToGetRidOf:
                                    listOfTypesToGetRidOf.append(a)
                            elif b.footPrint.issubset(a.footPrint): # Then we need to get rid of b
                                if b not in listOfTypesToGetRidOf:
                                    listOfTypesToGetRidOf.append(b)
                            else: # Here we have a scenario where the two footprints share intras but one is not a subset of the other
                                FINALTYPECLADECOLLECTION.isMixedIdentification = True
                    for TYPE in listOfTypesToGetRidOf:
                        FINALTYPECLADECOLLECTION.listOfFinalTypes.remove(TYPE)

                    # Once all of the unsupported (both by abundance and subsets) TYPES have been removed check to see if all of the samples > cutoff
                    # (for a given clade) intras are accounted for within the list of types identified

                    # e.g. if an inital type of D1-D7 that was unsupported
                    # end up as a final type of D1, then the D7 component of this symbiont has not been taken into account in the final type call
                    # in this case FINALTYPECLADECOLLECTION.identified would be FALSE
                    for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                        if CLADECOLLECTION.clade == CLADE: # Then this is the clade in questions cladeCollection and we can use the its2 occurances in this as the list of intras found at above the config.args.cutOff
                            # If True then all of the intras above the config.args.cutOff have been accounted for with the final types.
                            # FINALTYPECLADECOLLECTION.typeBasedCompCollection() returns a list of all of the seqs in all of the types in the final type list
                            if set([A.name for A in CLADECOLLECTION.listOfSeqsAboveCutOff]).issubset(FINALTYPECLADECOLLECTION.typeBasedCompCollection()): #and not FINALTYPECLADECOLLECTION.isMixedIdentification: # TODO consider the three type foot prints [1,2,3], [1,4,5] and [7,4,6] you can see the problem
                                FINALTYPECLADECOLLECTION.identified = True
                            else:
                                FINALTYPECLADECOLLECTION.identified = False


                            # Currently the types are written in the listOfFinalTypes as INITIAL types.
                            # So we go through each of the types and we write it out properly as a Final Type to the listOfFinalTypes
                            # I have replaced the typesupport value here with None rather than the typesupportDict that I had been making earlier. I don't see the point of this count.
                            i = 0
                            while i < len(FINALTYPECLADECOLLECTION.listOfFinalTypes):
                                FINALTYPECLADECOLLECTION.listOfFinalTypes[i] = symbiodiniumType(clade=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].clade, footPrint=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].footPrint, maj=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].maj, listofcodommajs=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].listofcodommajs, coDom=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].coDom, typeOfType='FINAL', totalSeqs=SAMPLE.totalSeqs, typeSupport=None, listofoccurences=SAMPLE.compComplement.listOfits2SequenceOccurances, name=FINALTYPECLADECOLLECTION.listOfFinalTypes[i].name)
                                i += 1
                            break # Break out of the for that checks to see if finaltypecladecollection is identified
                    break # Break out of the very first if loop that identifies a finaltypecladecollection of the given clade as there is only one per clade
    return

def addTypeToSamplesPermute(newSymType, listOfSamplesThatContainFootprint):
    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        if SAMPLE.name in listOfSamplesThatContainFootprint:
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                 if CLADECOLLECTION.clade == newSymType.clade and CLADECOLLECTION.footPrint == newSymType.footPrint: # Then this is a cladeCollection that we want to add the new Symbiodinium Type to
                    # We are going to add the new type to the sample. However we are going to add it as a new instance of the type so that we can have different sortedDefiningITSwOccurance for each sample
                    # if we add exactly the same occurance of the type to multiple samples then when we change one it will change them all.
                    # E.g. if we change the sortedDefining... for one of the samples it will automatically be changed in all of the samples with that type.
                    CLADECOLLECTION.addInitialType(symbiodiniumType(clade=newSymType.clade, coDom=newSymType.coDom, maj=newSymType.maj, footPrint=newSymType.footPrint, typeOfType=newSymType.typeOfType, listofSamples=newSymType.listOfSamples, listofcodommajs=newSymType.listofcodommajs))
                    # Here we add the CLADECOLLECTION.initialType.sortedDefiningIts2Occurances
                    CLADECOLLECTION.initialType.sortedDefiningIts2Occurances = CLADECOLLECTION.initialType.createSortedDefiningIts2Occurances(SAMPLE.compComplement.listOfits2SequenceOccurances, SAMPLE.totalSeqs)[0]
    return

def calculatePerformaceParametersPermute():
    # Number of initial supported types identified
    collectionOfMajs = []
    collectionOfInitialTypes = []
    # A total so that we can work out proportions
    numberOfCladeCollections = 0
    # Number of final type clade collections that have final types in them
    numberOfCladeCollectionsWithFinalTypes = 0
    # I'm gonna keep this super simple and say if any of the types share defining intras then we have an unresolved
    numberOfFinalTypeCladeCollectionUnsolved = 0

    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        for CLADECOLLECTION in SAMPLE.cladeCollectionList:
            numberOfCladeCollections += 1
            if CLADECOLLECTION.maj not in collectionOfMajs:
                collectionOfMajs.append(CLADECOLLECTION.maj)
            if CLADECOLLECTION.initialType.supportedType and CLADECOLLECTION.initialType.name not in collectionOfInitialTypes:
                collectionOfInitialTypes.append(CLADECOLLECTION.initialType.name)

        for FINALTYPECLADECOLLECTION in SAMPLE.finalTypeCladeCollectionList:
            unresolved = False
            numberOfCladeCollectionsWithFinalTypes += 1
            if len(FINALTYPECLADECOLLECTION.listOfFinalTypes) > 1: # THere must be at least two types for us to find an unsolved finaltypecladecollection
                for finalType1, finalType2 in itertools.combinations(FINALTYPECLADECOLLECTION.listOfFinalTypes,2):
                    for definingIntras in [a[0] for a in finalType1.sortedDefiningIts2Occurances]:
                        if definingIntras in [a[0] for a in finalType2.sortedDefiningIts2Occurances]:
                            unresolved = True
            if unresolved == True:
                numberOfFinalTypeCladeCollectionUnsolved += 1

    return len(collectionOfInitialTypes), len(collectionOfMajs), numberOfCladeCollectionsWithFinalTypes/numberOfCladeCollections, numberOfFinalTypeCladeCollectionUnsolved/numberOfCladeCollectionsWithFinalTypes

def clearConfigAbundanceList():
    for SAMPLEKEY in config.abundanceList.keys():
        SAMPLE = config.abundanceList[SAMPLEKEY]
        del SAMPLE.cladeCollectionList[:]
        del SAMPLE.finalTypeCladeCollectionList[:]

def assessTypeCollapse():
    '''
    Update all types' sortedITS2Occurances
    Generate groups for each of the majs
    Place each of the types into their respective groups
    Within groups assess each type pair for collapse
    Make note of which types collapse
    Where either of the types in a collapsing pair also collapse with alternative types but the other inital
    collapsing pair doesn't collapse with this then collapse the type with the maximum collapsing possibilities
    to the type with the closest Fst.
    Once a collapse has occured... lets get to here first and see what it looks like
    :return:
    '''

    # groupnames are the Maj types found in types
    # but if coDom types exist the groups contain all intras in the codoms that share intras
    groupNames = []
    for TYPENAME in config.typeDB.keys():
        TYPE = config.typeDB[TYPENAME]
        tempGroupList = [TYPE.majDict.keys()]
        assignTempGroup(groupNames, tempGroupList)

    a = 'partridge'



def assignTempGroup(groupNames, tempGroupList):
    for i in range(len(tempGroupList)):
        for j in range(len(groupNames)):
            if tempGroupList[i] in groupNames[j]:
                groupNames[j].update(tempGroupList)
                return
    groupNames.append(set(tempGroupList))
    return



def evaluateTypeSimilarity(listOfTypeNames):

    cwd = os.path.dirname(__file__)
    global abundanceList
    abundanceList = readByteObjectFromDefinedDirectory(
        '{0}/MEDdata/serialized objects/abundanceListWithFinalTypes'.format(cwd))
    global typeDB
    typeDB = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/typeDB'.format(cwd))
    global oursToLaJDict
    oursToLaJDict = readByteObjectFromDefinedDirectory(
        '{0}/MEDdata/serialized objects/oursToLaJDict'.format(cwd))

    listOfTypes = [typeDB[typeName] for typeName in listOfTypeNames]

    # Add the intras in the order form the largest type's footprint first
    # Get list of intras in order for each type
    orderedListOfListsOfIntras = []
    footPrintList = []
    for type in listOfTypes:
        footPrintList.append([a[0] for a in type.sortedDefiningIts2Occurances])
    orderedListOfListsOfIntras = sorted(footPrintList, key=len, reverse=True)
    a = 5

    # Add all intras in order of occurence in longest type first
    orderedListOfIntras = []
    for intraList in orderedListOfListsOfIntras:
        try:
            orderedListOfIntras.extend([intra for intra in intraList if intra not in orderedListOfIntras])
        except:
            pass

    typesInfo = []
    for type in listOfTypes:
        typesInfo.append(popIntraInfo(symType=type, orderedListOfIntras=orderedListOfIntras))

## MAIN FUNCTION ##
def CreateHumeFstMatrices():
    #For initial use of program

    # Sets all the global/argument parameters e.g. cutoffs and root directories etc.
    #Creates the master fastaDictionary a dictionary that contains the seq name and sequence for every sequence that was returned
    #Checks to see if the R libraries are installed correctly
    # Creates the oursToLaJDict which is a dict linking our sequences to those of the LaJ ITS2 seqs with names
    #Creates abundance list
        # At this point this contains a list of our samples and a compComplement class that contains a list of each of the
        # sequences found in each sample along with their frequency, name, sequence and clade
        # There is also an intraAbundanceDict which is simply a dict of sequence name to abundance
    config.__init__()


    '''
    Here I am going to implement an interation of the defining of the clade collections, inital types and final types
    with different cutoff values for which abundance of seqs are used in the definition of the inital types.
    I will calculate the number of types established, the proortion of final cladal collections which do not have final
    types and the proportion of finalcladecollections that have multiple final types that have the same majority.
    I will plot these performance parameters again the cutoff value and hopefully I will be able to work out
    a cutoff from these results which will allow me to dynamically select an appropriate cutoff value for the given
    dataset.
    '''

    # Implement a master loop for each of the cutoff values to be tested
    if config.args.conductDynamicCutoffAnalysis:
        print('Permuting defining sequence cutoffs')
        # emptyHolderCopyOfAbundnaceList = deepcopy(config.abundanceList)
        performanceParametersList = []
        for CUTOFF in np.arange(0.2, 0.00, -0.01):
            print('Permuting {0} cutoff'.format(CUTOFF))
            #TODO rather than doing this deepcopy business try writing simple loop that clears out the config.abundanceList
            # config.abundanceList = deepcopy(emptyHolderCopyOfAbundnaceList)
            p1, p2, p3, p4 = computeCutOffPerformanceParameters(CUTOFF)
            performanceParametersList.append(['{0}.2f'.format(CUTOFF) , p1, p2, p3, p4])
            clearConfigAbundanceList()
        # At this point we have computed all of the parameters we need to write the list out for plotting
        write2DListToDestination(config.args.saveLocation + '/performanceParameters/performanceParameters.txt', performanceParametersList)
        #TODO also count the number of samples that have supported initial types
        #TODO create a function that chooses the best cut off value from the information gained above
    # Here we put congif.abundanceList back to the empty copy for the last time

    # For each sample this assigns the cladeCollections. For each sample, for each clade, the program looks
    # to see if the given sample contains more than the args.cladeCollectionCutOff percentage of sequences from a given clade.
    # If it does then it will make a cladeCollection which is a subselection of the ITS2 sequences that that sample
    # contains that are from the given clade and are found at an abundance greater than the
    # args.cutOff relative to the number of sequences found in that sample for the clade in question.
    # E.g. if there are 300 clade A seqs out of a total seqs of 1000 then clade A seqs that were returned at an abundnace greater or equal to 30
    # (if the args.cutOff is set to 0.1) will be put into the clade A cladeCollection for that sample
    cladeCollectionCountDict = None
    print('Beginning clade collection assignment')
    if config.args.developingMode:
        try:
            config.abundanceList = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListCladeCollectionAssigned')
            cladeCollectionCountDict = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/cladeCollectionCountDict')
        except:
            print('Missing Object: abundanceListCladeCollectionAssigned not found in specified directory\n Creating from scratch...')
            cladeCollectionCountDict = assignCladeCollections()
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/bundanceListCladeCollectionAssigned', config.abundanceList)
    else:
        cladeCollectionCountDict = assignCladeCollections()
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListCladeCollectionAssigned', config.abundanceList)
    print('Clade collection assigned')


    # For each cladeCollection that is assigned to a sample, an itital type is assigned using the subset of sequences in that cladeCollection.
    # They are treated as a footprint. For every given footprint we look to see which other samples contain that footprint.
    # By this means we are able to define a footprint or type as supported or unsupported (more or less than the args.typeSupport).
    # We also check to see if the type is a coDom or not (a footprint or type that has several different ITS2 sequences as it majority sequence).
    # We use a slightly conservative approach for this. The type must have at least two different predominant ITS2 sequences and each
    # of those sequences must be found in at least args.coDomSupport number of samples.
    # Once these parameters have been identified the type is assigned to the cadeCollection in question and also to all cladeCollections that shared the same type.
    # If the type is a coDom we add the coDomSupportedMajDict which is a dictionary containing the name of the ITS2 seqs that are
    # MajITS2 seqs as key and true or false depending on whether they were found in > args.coDomSupport number of
    # samples as True or False to the instance of the Symbiodiniumtype class

    config.typeDB = symbiodiniumTypeDB()

    print('Assigning initial types')
    if config.args.developingMode:
        try:
            config.abundanceList = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListInitialTypesAssigned')
            config.typeDB = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/typeDB')
        except:
            print('Missing Object: abundanceListCladeCollectionAssigned not found in specified directory\n Creating from scratch...')
            assignInitialTypes(cladeCollectionCountDict)
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListInitialTypesAssigned', config.abundanceList)
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/typeDB', config.typeDB)
    else:
        assignInitialTypes(cladeCollectionCountDict)  # This now assigns initial types to the samples' cladeCollections
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListInitialTypesAssigned', config.abundanceList)
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/typeDB', config.typeDB)
    print('Assigned initial types')

    ######
    # At this point we have the initial types identified. We also have the fst values identified which are what we need for the graphical outputs.
    # From this point on we have to:
    # Recalculate the final types and calculate their support.
    # Create the matrices, header files and info files to plot data in R
    # Write all of the above to the html output.

    # Append information to the samples within the abundanceList refering to the finalTypeCladeCollections
    # This forces all detected footprints (from the initial types) into the compComplement list of sequences for each sample to see which footprints can be found.
    # I.e. we are now using all seqs, not just those above the config.args.cutOff to look for the footprints already identified
    # We create a listOfFinalTypes list for each cladeCollection that contains all of the footprints found in the compcomplement
    # We then get rid of seqeunces in this list if they are now unsupported (we re-count to see how many samples each type is in)
    # (i.e. fall below config.args.typeSupport) and we get rid of the shorter footprint (less intras) if it is a subset of another type in the list.
    # If we have a scenario where footprints that fit into a sample's comComplement share intras but neither one is a subset of the other
    # then we make isMixedIdentification true.
    # Once we have knocked out unsupported types or types that are subsets of others we finaltypecladecollection.identified to either True or False
    # depending on whether all of the intras found above the config.args.cutOff have been identified by the types in the finalTypeCladeCollection.listOfFinalTypes
    # and there is no conflict between footprints that share intras but are not subsets.

    print('Infering final symbiondinium types')
    if config.args.developingMode:
        try:
            config.abundanceList = readByteObjectFromDefinedDirectory(
                config.args.saveLocation + r'/serialized objects/abundanceListWithFinalTypes')
            config.typeDB = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/typeDB')
        except:
            print(
                'Missing Object: abundanceListWithFinalTypes  not found in specified directory\n Creating from scratch...')
            inferFinalSymbiodiniumTypesIterative()
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListWithFinalTypes', config.abundanceList)
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/typeDB', config.typeDB)
    else:
        inferFinalSymbiodiniumTypesIterative()
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/abundanceListWithFinalTypes', config.abundanceList)
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/typeDB', config.typeDB)
    print('Final type inference complete')


    assessTypeCollapse()

    #Create masterSeqDistancesDict
    #This is a dictionary that has the genetic distances between every combination of sequences from the same clade.
    # The keys are frozensets so that the order in which you look for the two sequences doesn't matter. The value is the genetic distance
    masterSeqDistancesDict = None
    print('Calculating sequence distances')
    if config.args.createMasterSeqDistancesFromScratch == False or config.args.developingMode:
        try:
            masterSeqDistancesDict = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/masterSeqDistancesDict')
        except:
            print('Missing Object: masterSeqDistancesDict  not found in specified directory\n Creating from scratch...')
            masterSeqDistancesDict = createMasterSeqDistancesNonMothur()
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/masterSeqDistancesDict', masterSeqDistancesDict)
    else:
        masterSeqDistancesDict = createMasterSeqDistancesNonMothur()
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/masterSeqDistancesDict',
                                          masterSeqDistancesDict)
    print('Sequence distance calculation complete')


    finalLogTransedColDists = None
    print('Calculating Fst distances')
    if config.args.createFinalFstColDistsFromScratch == False or config.args.developingMode:
        try:
            finalLogTransedColDists = readByteObjectFromDefinedDirectory(config.args.saveLocation + r'/serialized objects/finalLogTransedColDists')
        except:
            print('Missing Object: finalLogTransedColDists not found in specified directory\n Creating from scratch...')
            finalLogTransedColDists = createLogTransedFstColDistsFromFinalTypes(masterSeqDistancesDict)
            writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/finalLogTransedColDists',
                                              finalLogTransedColDists)
    else:
        finalLogTransedColDists = createLogTransedFstColDistsFromFinalTypes(masterSeqDistancesDict)
        writeByteObjectToDefinedDirectory(config.args.saveLocation + r'/serialized objects/finalLogTransedColDists',
                                          finalLogTransedColDists)
    print('Fst distances calculated')



    # ######
    # # At this point we have the initial and final types identified as well as their supports. We also have the fst values identified which are what we need for the graphical outputs.
    # # From this point on we have to:
    # # Create the matrices, header files and info files to plot data in R
    # # Write all of the above to the html output.
    #
    # # Write out the graphics so that the writeTypeBasedOutput can write them into the html
    # # This produces plots for cladal and subcladal/maj levels
    # # For the cladal it produces both a Pareto and a 2D PC ordination of the calculated Fst values
    # print('Producing plots')
    # if config.args.developingMode:
    #     try:
    #         config.abundanceList = readByteObjectFromDefinedDirectory(
    #             config.args.saveLocation + r'\serialized objects', 'abundanceListAfterPlotCreation')
    #     except:
    #         print(
    #             'Missing Object: abundanceListAfterPlotCreation  not found in specified directory\n Creating from scratch...')
    #         producePlots(finalLogTransedColDists)
    #         writeByteObjectToDefinedDirectory(config.args.saveLocation + r'\serialized objects',
    #                                           'abundanceListAfterPlotCreation', config.abundanceList)
    # else:
    #     producePlots(finalLogTransedColDists)
    #     writeByteObjectToDefinedDirectory(config.args.saveLocation + r'\serialized objects',
    #                                       'abundanceListAfterPlotCreation', config.abundanceList)

    #TODO incorporate at somepoint
    # multiModalDetection()


    # This simply creates a list that is the start of the html eventual outputfile
    htmlOutput = createHtmlHolder()



    # Write the type based output
    # This function creates a list that contains lines of the html output
    # It is the information relating to the type-based output, i.e. how common each of the types are and which species they are found in etc.
    # It returns this html string as well as two dictionaries which are the number of times given types are found in the inital and final types lists
    htmlMainString = writeTypeBasedOutput2()

    htmlOutput.extend(htmlMainString)
    # Write the sample based output
    htmlOutput.extend(writeSampleCharacterisationOutput())
    # Add the closing tags to the html file
    closeAndWriteHtmlHolder(htmlOutput)



    
    if config.args.deleteIntermediateFiles:
        print('Deleting intermediate files')
        for dir in ['.bat scripts', 'matrices outputs', 'seq distance files', 'serialized objects']:
            path = os.path.join(config.args.saveLocation, dir)
            print('Deleting %s' % path)
            shutil.rmtree(path)
    if config.args.archiveInputs:
        print('Zipping up input files to save space')
        zf = zipfile.ZipFile(os.path.join(config.args.inputLocation, 'inputs.zip'), mode='w')
        threeFiles = ['ITS2Abundance.txt', 'ITS2Fasta.fasta', 'LaJeunesse Types.fas']
        inputs = [os.path.join(config.args.inputLocation, x) for x in threeFiles]
        for file in inputs:
            zf.write(file, arcname=os.path.basename(file), compress_type=zipfile.ZIP_DEFLATED)
        zf.close()
        for file in inputs:
            print('Deleting %s' % file)
            os.remove(file)

    investigateFinalTypeCorrelations()

    print('Program complete')

## MAIN ENTRY POINT OF PROGRAM ##
if __name__ == '__main__':
    CreateHumeFstMatrices()




