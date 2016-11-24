import csv
import os
import statistics
import numpy as np
import matplotlib.pyplot as plt
path = r'C:/Users/HUMEBC/Google Drive/EdData/screwaroundpsba/'
base = r'C:/Users/HUMEBC/Google Drive/EdData/'
import pickle
from subprocess import call
import copy

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

def writeByteObjectToDefinedDirectory(directory, object):
    os.makedirs(os.path.dirname(directory), exist_ok=True)

    f = open(directory, 'wb+')
    pickle.dump(object, f)

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


def stdInToMEDIn(cladeList):
    try:
        fastaList = readByteObjectFromDefinedDirectory('{0}/inputsForMED/cache/fastaList'.format(cd))
        return fastaList
    except:
        pass

    countTabIn = readDefinedFileTo2DList(filename='{0}/inputsForMED/countTab'.format(cd))

    fastaIn = readDefinedFileToList(filename='{0}/inputsForMED/ITS2Fasta.fas'.format(cd))

    fastaInDict = createDictFromFasta(fastaIn)

    cladeDict = readByteObjectFromDefinedDirectory('{0}/originalData/serialized objects/seqNameToCladeDict'.format(cd))[0]

    listOfSamples = [samp.replace('_', '-') for samp in countTabIn[0][1:]]
    listOfSeqs = [countTabIn[i][0] for i in range(1, len(countTabIn))]
    # Simply go through the countTab sample by sample
    # Each count you come to add that many of the seq to the cladal fasta
    fastaList = [[] for i in range(len(cladeList))]
    seqCounter = 1

    # # Create raw matrix
    # copyOf = readDefinedFileTo2DList(filename='{0}/inputsForMED/countTab'.format(cd))
    # del copyOf[0]
    # for seq in range(len(copyOf)):
    #     del copyOf[seq][0]
    # nparray = np.array(copyOf)
    # for index, val in np.ndenumerate(nparray):
    #     if int(val) != 0:
    #         try:
    #             indexOfClade = cladeList.index(cladeDict[listOfSeqs[index[0]]])
    #         except:
    #             continue
    #         fastaseq = fastaInDict[listOfSeqs[index[0]]]
    #         sample = listOfSamples[index[1]]
    #
    #         for i in range(seqCounter, seqCounter + int(val)):
    #             fastaList[indexOfClade].extend(['>{0}_seq{1}'.format(sample, i), fastaseq])
    #         seqCounter += int(val)


    for sample in range(1, len(countTabIn[0])):
        samplename = listOfSamples[sample - 1]
        for seq in range(1, len(countTabIn)):
            if int(countTabIn[seq][sample]) != 0:
                cladeOfSeq = cladeDict[countTabIn[seq][0]]
                try:
                    indexOfClade = cladeList.index(cladeOfSeq)
                except:
                    continue
                sequence = fastaInDict[countTabIn[seq][0]]
                val = int(countTabIn[seq][sample])
                for i in range(seqCounter, seqCounter + val):
                    fastaList[indexOfClade].extend(['>{0}_seq{1}'.format(samplename, i), sequence])
                seqCounter += val
        print('Processing sample {0}'.format(samplename))
    writeByteObjectToDefinedDirectory('{0}/inputsForMED/cache/fastaList'.format(cd), fastaList)
    # # First I will need to create a set of fasta file and countTab for each of the clades in question


    # # Get list of seqs in each clade
    # # Get list of samples that contain a seq for a given clade
    # # create countDict; key = sample/seq value = abundance; only when seq found, i.e. no 0s
    # seqAndSampList = []
    # countDict = {} # Keep track of which seqs samples contain
    # for CLADE in cladeList:
    #     listOfSeqs = [kys for kys in cladeDict if cladeDict[kys] == CLADE]
    #     listOfSamples = []
    #     for sample in range(1,len(countTabIn[0])):
    #         hasClade = False # If sample has any seq form current clade
    #         for its2Seq in range(1, len(countTabIn)):
    #             if int(countTabIn[its2Seq][sample]) != 0:
    #                 hasClade = True
    #                 countDict['{0}/{1}'.format(sample, its2Seq)] = countTabIn[its2Seq][sample]
    #         if hasClade:
    #             listOfSamples.append(sample)
    #     seqAndSampList.append((listOfSeqs, listOfSamples))
    #
    # # Create count tab and fasta for each clade
    # fasAndCountTabList = []
    # for i in range(len(cladeList)):
    #     infoTuple = seqAndSampList[i]
    #     # newFasta
    #     newFasta = []
    #     for seq in listOfSeqs:
    #         newFasta.append('>{0}'.format(seq))
    #         newFasta.append(fastaInDict[seq])
    #     # empty countTab
    #     emptyCountTab = initCountTab(infoTuple)
    #     # Populate countTab using countDict
    #     for sample in range(1,len(infoTuple[1]) + 1):
    #         for seq in range(1, len(infoTuple[0]) + 1):
    #             try:
    #                 emptyCountTab[seq][sample] = int(countDict['{0}/{1}'.format(sample, seq)])
    #             except:
    #                 emptyCountTab[seq][sample] = 0
    #     fasAndCountTabList.append((newFasta, emptyCountTab))

    return fastaList


def initCountTab(seqAndSampListItem):
    emptyCountTab = [[] for j in len(seqAndSampListItem[0]) + 1]
    emptyCountTab[0].extend(seqAndSampListItem[1])
    emptyCountTab[0].insert(0, 'COUNT')
    for j in range(1, len(emptyCountTab)):
        emptyCountTab[j].append(seqAndSampListItem[0][j])
        emptyCountTab[j].append('x' for k in range(len(seqAndSampListItem[1])))
    return emptyCountTab

def combineMultiCladeMED():
    cladeList = ['A', 'C', 'D']
    cd = os.path.dirname(__file__)
    dirList = []
    for clade in cladeList:
        dirList.append('fastaForMed{0}'.format(clade))

    # Get MED output files, fasta and countTabs
    symTyperCounts = []
    symTyperFastas = []
    for directory in dirList:
        symTyperCounts.append(medToCountTab('{0}/inputsForMED{1}/MATRIX-COUNT.txt'.format(cd, directory)))
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


def stdInputToMEDInput():
    global cd
    cd = os.path.dirname(__file__)

    cladeList = ['A', 'C', 'D']

    # Create a fasta for each clade with a sequence per read in each sample
    fastaList = stdInToMEDIn(cladeList)


    # Write the fastas to file
    listOfFastaNames = []
    for i in range(len(cladeList)):
        writeListToDestination('{0}/inputsForMED/MEDFasta{1}.fas'.format(cd, cladeList[i]), fastaList[i])
        listOfFastaNames.append('MEDFasta{0}.fas'.format(cladeList[i]))
    # Write shell script to run the MED on each of the MEDfasta files
    sS = []
    sS.append('#! /bin/sh')
    for i in range(len(listOfFastaNames)):
        sS.append('decompose {0}'.format(listOfFastaNames[i]))
    writeListToDestination('{0}/inputsForMED/MEDscript'.format(cd), sS)


    # Execute MED on all fasta files via console
    cmd = ['cd ~/fstProjectWorking/inputsForMED', 'sh ./MEDscript']
    call(['sh', './inputsForMED/MEDscript'])

    a = 7
def prodIntraPlotSimilarTypes(listOfTypeNames):

    cwd = os.path.dirname(__file__)
    global abundanceList
    abundanceList = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/abundanceListWithFinalTypes'.format(cwd))
    global typeDB
    typeDB = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/typeDB'.format(cwd))
    global oursToLaJDict
    oursToLaJDict = readByteObjectFromDefinedDirectory('{0}/MEDdata/serialized objects/oursToLaJDict'.format(cwd))

    listOfTypes = [typeDB[typeName] for typeName in listOfTypeNames]
    TYPEONE = typeDB['C3-ST-seq170']
    TYPETWO = typeDB['C3-ST-seq178-seq170']

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
        typesInfo.append(popIntraInfo(symType = type, orderedListOfIntras=orderedListOfIntras))


    # the x locations for the bar pairs
    adjust = len(listOfTypeNames)/2 # Change for space between groups of bars
    ind = np.arange(0, len(orderedListOfIntras)*adjust, adjust)
    width = 0.35

    # the subplot to plot on
    fig, ax = plt.subplots()

    # draw bars
    barsList = []
    colourList = ['b', 'g', 'r', 'c', 'm', 'y', 'orange']
    for i in range(len(listOfTypeNames)):
        barsList.append(ax.bar(ind + (width * i), typesInfo[i][2], width, color=colourList[i % 7], yerr=typesInfo[i][3]))


    # text and labels etc.
    ax.set_ylabel('abundance ratio to majoirty intra')
    ax.set_title('Comparison of intra ratios in similar types')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(tuple([CLJ(intra) for intra in orderedListOfIntras]))
    ax.set_ylim(bottom = 0, top=2)

    # legend
    ax.legend(tuple(bar[0] for bar in barsList), tuple(listOfTypeNames))

    plt.show()

    a = 6

def CLJ(VAR):
    try:
        if len(oursToLaJDict[VAR]) > 1:
            return '#'.join(oursToLaJDict[VAR])
        else:
            return oursToLaJDict[VAR][0]
    except:
        return VAR

def popIntraInfo(symType, orderedListOfIntras):
    # Compile info for Type1 intra means and std

    intrainfolist = [[]for i in range(4)]
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
        intrainfolist[0].append(tempPropList)
        # Divide the latest set of proportions by the first intra for each
        # sample to get the ratio to append in typeOneInfo[1]
        intrainfolist[1].append(
            [(intrainfolist[0][-1][i] / intrainfolist[0][0][i]) for i in range(len(symType.samplesFoundInAsFinal))])
        # Transform anyratios above 1
        for i in range(len(intrainfolist[1][-1])):
            RATIO = intrainfolist[1][-1][i]
            if RATIO > 1:
                intrainfolist[1][-1][i] = 1 + (1 - (1 / RATIO))
        intrainfolist[2].append(sum(intrainfolist[1][-1]) / len(intrainfolist[1][-1]))
        try:
            intrainfolist[3].append(statistics.stdev(intrainfolist[1][-1]))
        except: # If e.g. ther is only one sample
            intrainfolist[3].append(0)
    return intrainfolist

# prodIntraPlotSimilarTypes(['C3-ST-seq170', 'C3-ST-seq178-seq170', 'C3-seq220-seq189-C107', 'C3-seq170', 'C3-ST-seq178-seq242', 'C3'])
# prodIntraPlotSimilarTypes(['C39-C1-seq180', 'C39-seq186-seq180-C1', 'C39-seq186-seq180-C1-seq214', 'C39-seq180-C1-seq207', 'C39-seq186-seq213-seq180', 'C39-seq180-seq238', 'C39'])
# prodIntraPlotSimilarTypes(['C3-ST-seq170', 'C3-ST-seq178-seq170', 'C3/ST', 'C3/seq178-seq170', 'C3/seq178-ST', 'C3/seq178-ST-seq252'])
# prodIntraPlotSimilarTypes(['Cseq1-seq168-C1-seq176', 'Cseq1-seq168-C1-seq184-seq176', 'Cseq1-seq173-seq168-C1-seq198', 'C1/Cseq1-C39', 'Cseq1/C1-seq176'])
# prodIntraPlotSimilarTypes(['C15', 'C15-seq236', 'C15-seq204', 'C15n/C15', 'C15/seq196'])
# prodIntraPlotSimilarTypes(['D1-D4-seq417', 'D1-seq433-D4', 'D1-seq409-D4-D6-seq411-seq415', 'D1/D4', 'D4/D1-D6-D2-seq411-seq414', 'D1/D4-D6-seq411', 'D1/D4-seq409-D6-D2'])

#TODO so the problem appears to be that MED can't access the fasta file for somereason

# stdInputToMEDInput()
a = 'joe'
cd = os.getcwd()
# os.chmod('{0}/inputsForMED/MEDFastaA.fas'.format(cd), 0o755)
# os.system('chmod 0755 ./inputsForMED/MEDscript')
os.system('less test')
os.system('decompose ~/fstProjectWorking/inputsForMED/MEDFASTA.fas-PADDED-WITH-GAPS')
# call(['sh', './inputsForMED/script.txt'])
# a = 6
# call(['sh', './inputsForMED/MEDscript'])
# a = 5