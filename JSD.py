import csv
import os
import statistics
import numpy as np
import matplotlib.pyplot as plt
path = r'C:/Users/HUMEBC/Google Drive/EdData/screwaroundpsba/'
base = r'C:/Users/HUMEBC/Google Drive/EdData/'
import pickle

def readDefinedFileToList(filename):
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

def writeListToDestination(destination, listToWrite):
    print('Writing list to ' + destination)
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

def createDictFromFasta(fastaList):
    tempDict = {}
    i = 0
    while i < len(fastaList):
        tempDict[fastaList[i][1:]] = fastaList[i+1]
        i += 2
    return tempDict

def readDefinedFileTo2DList(filename):
    tempListPrimary = []
    tempListSecondary = []

    with open(filename, mode='r') as reader:
        for line in reader:
            tempListSecondary = line.rstrip().split(',')
            tempListPrimary.append(tempListSecondary)
    return tempListPrimary

def write2DListToDestination(destination, twoDlistToWrite):

    with open(destination, mode='w+') as writer:
        i = 0
        while i < len(twoDlistToWrite):
            if i == len(twoDlistToWrite)-1:
                writer.write(','.join([str(a) for a in twoDlistToWrite[i]]))
            else:
                writer.write(','.join([str(a) for a in twoDlistToWrite[i]]) + '\n')
            i += 1

def readByteObjectFromDefinedDirectory(object):
    f = open(object, 'rb')
    return pickle.load(f)

def readTabDelim(path):
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        d = list(reader)
    return d

def addToNewFas(seq, samp, abun, newFas, count, listID):
    for i in range(abun):
        newFas[listID].extend(['>{0}_{1}'.format(samp, str(count[listID])),
                       '{0}{1}'.format(fastaDict[seq], '-' * (largestSeqLen - len(fastaDict[seq])))])
        count[listID] += 1
    return count

def addSeqOccurToNewFas(seq, samp, abun, newFas, countM):
    if cladalDict[seq] == 'A':
        countM = addToNewFas(seq, samp, abun, newFas, countM, 0)
    elif cladalDict[seq] == 'B':
        countM = addToNewFas(seq, samp, abun, newFas, countM, 1)
    elif cladalDict[seq] == 'C':
        countM = addToNewFas(seq, samp, abun, newFas, countM, 2)
    elif cladalDict[seq] == 'D':
        countM = addToNewFas(seq, samp, abun, newFas, countM, 3)


    return countM

def prepMEDFastaForSymTyperInput(MEDFastaFileDir, outputfiledir = None):

    tempFasta = readDefinedFileToList(MEDFastaFileDir)

    i = 0
    while i < len(tempFasta):
        tempFasta[i] = tempFasta[i].split('|')[0]
        i += 2
    if outputfiledir != None:
        writeListToDestination(outputfiledir, tempFasta)
    else:
        return tempFasta

def medToCountTab(inputfile, outputfile = None):
    # read in MED count file
    edDataCount = readTabDelim(inputfile)
    # transpose count file
    edT = [list(i) for i in zip(*edDataCount)]
    if outputfile != None:
        # write file to dir
        write2DListToDestination(outputfile, edT)
    else:
        return edT

def createMEDInputFas(cladedictpicklefile, rawabundancefiledir, fastafiledir, outputdir):
    newfastaFile = [[] for i in range(4)]
    fasta = readDefinedFileToList('{0}'.format(fastafiledir))
    global fastaDict
    fastaDict = createDictFromFasta(fasta)
    global largestSeqLen
    largestSeqLen = len([a[1] for a in sorted(fastaDict.items(), key=lambda x:len(x[1]), reverse=True)][0])
    global cladalDict
    cladalDict = readByteObjectFromDefinedDirectory(cladedictpicklefile)[0]
    abundanceData = readDefinedFileTo2DList(rawabundancefiledir)

    # For each seq in each sample
    # Check clade of seq and add to one of the newfastaFile lists according to clade
    # Run a count so that for each sample the sequences are numbered as they are detected
    # Reset this count at each new sample
    for sample in range(1,len(abundanceData[0])):
        print('Processing {0}'.format(abundanceData[0][sample]))
        count = [1, 1, 1, 1]
        for Otu in range(1,len(abundanceData)):
            abun = int(abundanceData[Otu][sample])
            if abun != 0:
                count = addSeqOccurToNewFas(abundanceData[Otu][0], abundanceData[0][sample], abun, newfastaFile, count)
    cladeList = ['A', 'B', 'C', 'D']
    for i in range(len(newfastaFile)):
        writeListToDestination('{0}/fastaForMed{1}.fas'.format(outputdir, cladeList[i]), newfastaFile[i])
    return


def combineMultiCladeMED():
    cladeList = ['A', 'C', 'D']
    cd = os.getcwd()
    dirList = []
    for clade in cladeList:
        dirList.append('fastaForMed{0}'.format(clade))

    # Get MED output files, fasta and countTabs
    symTyperCounts = []
    symTyperFastas = []
    for directory in dirList:
        symTyperCounts.append(medToCountTab('{0}/rawData/testingMED/{1}/MATRIX-COUNT.txt'.format(cd, directory)))
        symTyperFastas.append(prepMEDFastaForSymTyperInput('{0}/rawData/testingMED/{1}/NODE-REPRESENTATIVES.fasta'.format(cd, directory)))
    a = 5

    # because each of the MED analyses will use the same  sequence names we need to make sure these are unique
    # Go through each of the fasta and count files translating to unique seqnumbers
    seqNum = 1
    for i in range(len(symTyperCounts)): # For each clade
        seqConvertDict = {}


        # Convert fasta to unique seq names, note conversion
        for j in range(len(symTyperFastas[i])): # For each fasta line
            if symTyperFastas[i][j][0] == '>': # This is a name line
                seqConvertDict[symTyperFastas[i][j][1:]] = seqNum # Note conversion
                symTyperFastas[i][j] = '>seq{0}'.format(str(seqNum)) # Convert
                seqNum += 1 # Update seqNum


        # Convert count table to unique seq names
        for j in range(1,len(symTyperCounts[i])): # Each seq in Count
            symTyperCounts[i][j][0] = seqConvertDict[symTyperCounts[i][j][0]]


        # Create dict from each tab where key = Sample/Seq:Abun
        # Only entries where seq found
        seqCountDict = {}
        listOfSamples = []
        listOfSeqs = []
        for i in range(len(symTyperCounts)):  # For each clade
            countTab = symTyperCounts[i]
            listOfSamples.extend([a for a in countTab[0][1:] if a not in listOfSamples]) # Update sampleList
            listOfSeqs.extend([a[0] for a in countTab[1:] if a[0] not in listOfSeqs]) # Update seqList
            for sample in range(1, len(countTab[0])):
                for sequence in range(1,len(countTab)):
                    if countTab[sequence][sample] != 0:
                        seqCountDict['{0}/{1}'.format(countTab[0][sample], countTab[sequence][0])] = int(countTab[sequence][sample])


    # Blank count tab
    countTab = [[] for sequences in listOfSeqs]
    countTab.insert(0, ['countTab'])
    countTab[0].extend([listOfSamples])
    for i in range(1, len(countTab)): # Pop Seq names
        countTab[i].append(listOfSeqs[i])
        countTab[i].extend(['x' for i in range(len(listOfSamples))])


    # Populate count tab
    for sample in range(1, len(countTab[0])):
        for sequence in range(1, len(countTab)):
            try:
                countTab[sequence][sample] = seqCountDict['{0}/{1}'.format(countTab[0][sample], countTab[sequence][0])]
            except:
                countTab[sequence][sample] = 0


    # Consolidate fastas
    completeMEDFasta = []
    for i in range(len(symTyperFastas)):
        completeMEDFasta.extend(symTyperFastas[i])


    # Write
    write2DListToDestination('{0}/rawData/testingMED/completeMedCountTable', countTab)
    write2DListToDestination('{0}/rawData/testingMED/completeMedFasta', completeMEDFasta)


def prodIntraPlotSimilarTypes():
    cwd = os.path.dirname(__file__)
    global abundanceList
    abundanceList = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/abundanceListWithFinalTypes'.format(cwd))
    global typeDB
    typeDB = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/typeDB'.format(cwd))
    oTLJD = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/oursToLaJDict'.format(cwd))

    TYPEONE = typeDB['C3-ST-seq170']
    TYPETWO = typeDB['C3-ST-seq178-seq170']

    # Add the intras in the order form the largest type's footprint first
    orderedListOfIntras = max([[a[0] for a in TYPEONE.sortedDefiningIts2Occurances], [a[0] for a in TYPETWO.sortedDefiningIts2Occurances]], key=len )
    # Then add any remaining intras if not already in list
    try:
        orderedListOfIntras.extend([intra for intra in [a[0] for a in TYPEONE.sortedDefiningIts2Occurances] if intra not in b])
    except:
        pass
    try:
        orderedListOfIntras.extend([intra for intra in [a[0] for a in TYPETWO.sortedDefiningIts2Occurances] if intra not in b])
    except:
        pass

    typeOneInfo = [[] for i in range(4)]
    typeTwoInfo = [[] for i in range(4)]

    typeOneInfo = popIntraInfo(TYPEONE, orderedListOfIntras, typeOneInfo)
    typeTwoInfo = popIntraInfo(TYPETWO, orderedListOfIntras, typeTwoInfo)

    # Here we have the info we need to make the pots
    typeOneMeans = typeOneInfo[2]
    typeOneStDevs = typeOneInfo[3]
    typeTwoMeans = typeTwoInfo[2]
    typeTwoStDevs = typeTwoInfo[3]

    # the x locations for the bar pairs
    ind = np.arrange(len(orderedListOfIntras))
    width = 0.35

    # the subplot to plot on
    fig, ax = plt.subplots()

    # draw bars
    bars1 = ax.bar(ind, typeOneMeans, width, color = 'r', yerr=typeOneStDevs)
    bars2 = ax.bar(ind, typeTwoMeans, width, color = 'y', yerr= typeTwoStDevs)

    # text and labels etc.
    ax.set_ylabel('abundance ratio to majoirty intra')
    ax.set_title('Comparison of intra ratios in similar types')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(tuple(orderedListOfIntras))

def popIntraInfo(symType, orderedListOfIntras, typeOneInfo):
    # Compile info for Type1 intra means and std
    for intra in orderedListOfIntras:
        tempPropList = []
        for SAMPLENAME in symType.samplesFoundInAsFinal:
            SAMPLE = abundanceList[SAMPLENAME]
            for CLADECOLLECTION in SAMPLE.cladeCollectionList:
                if CLADECOLLECTION.clade == symType.clade:
                    try:
                        tempPropList.append(
                            SAMPLE.intraAbundanceDict[intra] / (CLADECOLLECTION.cladalProportion * SAMPLE.totalSeqs))
                    except:
                        tempPropList.append(0)

        # Put in next set of proportions for intra and calculate ratios
        typeOneInfo[0].append(tempPropList)
        # Divide the latest set of proportions by the first intra for each
        # sample to get the ratio to append in typeOneInfo[1]
        typeOneInfo[1].append(
            [(typeOneInfo[0][-1][i] / typeOneInfo[0][0][i]) for i in range(len(symType.samplesFoundInAsFinal))])
        # Transform anyratios above 1
        for i in range(len(typeOneInfo[1][-1])):
            RATIO = typeOneInfo[1][-1][i]
            if RATIO > 1:
                typeOneInfo[1][-1][i] = 1 + (1 - (1 / RATIO))
        typeOneInfo[2].append(sum(typeOneInfo[1][-1]) / len(typeOneInfo[1][-1]))
        typeOneInfo[3].append(statistics.stdev(typeOneInfo[1][-1]))
    return typeOneInfo

prodIntraPlotSimilarTypes()