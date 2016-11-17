import argparse

import pickle
from multiprocessing import Queue, Process
import multiprocessing
import os
from subprocess import call
import sys

import numpy as np
import itertools

import timeit

from scipy.stats import entropy

from ClassesTwo import *

def __init__():
    parser = argparse.ArgumentParser()
    
    # Computation parameters
    parser.add_argument('--conductDynamicCutoffAnalysis', type=bool, default=False, metavar='dynamCutOffAnalysis')
    parser.add_argument('--cutOff', type=float, default=0.04, metavar='CUTOFF')
    parser.add_argument('--cladePlotCutOff', type=int, default=3, metavar='cladePlotCutOff')
    parser.add_argument('--majPlotCutOff', type=int, default=4, metavar='majPlotCutOff')
    parser.add_argument('--typeSupport', type=int, default=0.01, help='The number of samples that must contain a type in order for it to be considered as a genuine type as a proportion of the number of cladeCollections for a given clade', metavar='N')
    parser.add_argument('--coDomSupport', type=int, default=2, help='The number of samples that must contain the an alternative MajITS2 sequence compared to the other samples of the type in order for it be considered a codom type', metavar='N')
    parser.add_argument('--cladeCollectionCutoff', type=float, default=0.1, metavar='CUTOFF')
    parser.add_argument('--cladeList', nargs='*', default=['A', 'B', 'C', 'D'], choices=['A', 'B', 'C', 'D'])
    parser.add_argument('--numProcessors', type=int, default=3, metavar='N')
    parser.add_argument('--majPlotThreshold', type=float, default=2, help='The number of supported types required within a Maj in order for its plotting data to be made', metavar='N')
    
    # Paths
    cwd = os.path.dirname(__file__)
    edBasePath = r'C:/Users/HUMEBC/Google Drive/EdData/screwaround'

    parser.add_argument('--rootLocation', default=cwd, help='Directory where the source code is found', metavar='PATH')
    # parser.add_argument('--inputLocation', default=cwd + '/raw data', help='Directory where the three input files are found', metavar='PATH')
    parser.add_argument('--inputLocation', default=edBasePath + '/raw data',help='Directory where the three input files are found', metavar='PATH')
    # parser.add_argument('--saveLocation', default=cwd, help='Output directory for saving matrices and output tables', metavar='PATH')
    parser.add_argument('--saveLocation', default=edBasePath, help='Output directory for saving matrices and output tables',metavar='PATH')
    # parser.add_argument('--mothurLocation', default=cwd + '/Mothur/mothur.exe', help='Full path including mothur.exe', metavar='PATH')
    parser.add_argument('--rscriptLocation', default=cwd + '/R/R-3.3.0/bin/x64/Rscript.exe', help='Full path including Rscript.exe', metavar='PATH')
    
    # Batch mode
    parser.add_argument('--logToFile', action='store_true', help='Redirect stdout/stderr to files in saveLocation')
    parser.add_argument('--deleteIntermediateFiles', action='store_true', help='When finished, delete intermediate state leaving only input/output')
    parser.add_argument('--archiveInputs', action='store_true', help='When finished, compress input files into a .zip, and delete originals')
    
    # Caching
    parser.add_argument('--developingMode', type = bool, default =False, metavar='TRUE|FALSE')
    parser.add_argument('--createAbundanceListFromScratch', type=bool, default=True, metavar='TRUE|FALSE')
    parser.add_argument('--createSeqToCladeDictFromScratch', type=bool, default=True, metavar='TRUE|FALSE')
    parser.add_argument('--createFinalFstColDistsFromScratch', type=bool, default=True, metavar='TRUE|FALSE')
    parser.add_argument('--createMasterSeqDistancesFromScratch', type=bool, default=True, metavar='TRUE|FALSE')
    # parser.add_argument('--createFstColDistsFromScratch', type=bool, default=False, metavar='TRUE|FALSE')
    parser.add_argument('--createOursToLaJDictFromScratch', type=bool, default=True, metavar='TRUE|FALSE')
    
    global args
    args = parser.parse_args()
    
    # Before we do anything else, redirect stdout/sterr if requested
    # http://stackoverflow.com/a/4110906/1688738
    if args.logToFile:
        sys.stdout = open(args.saveLocation + '/stdout.txt', 'w')
        sys.stderr = open(args.saveLocation + '/stderr.txt', 'w')

    # The masterFastaDict is simply a dict of the fasta file that contains all of the seqs that came out
    global masterFastaDict
    masterFastaDict = createMasterFastaDict()


    # Check to see if the R libraries required are installed
    print('Checking for R libraries...')
    RInstalledCorrectly = True
    RInstalledCorrectly = checkRLibraries(RInstalledCorrectly)

    # Print when sucessful
    # Go through the directory and check to see if any of the libraries are missing. If any are missing then RInstalledCorrectly = False

    if RInstalledCorrectly:
        print('Groovy, R libraries have been verified')
    else: # not RInstalledCorrectly:
        print('Ooops, R libraries missing. Installing...')
        call([args.rscriptLocation, args.rootLocation + r"\.bat scripts\installRLibs.bat", args.rootLocation.replace('\\', '/')], shell=True)

    global oursToLaJDict
    oursToLaJDict = None
    if not args.createOursToLaJDictFromScratch or args.developingMode:
        try:
            oursToLaJDict = readByteObjectFromDefinedDirectory(args.saveLocation + r'\serialized objects', 'oursToLaJDict')
        except:
            print('Missing Object', 'oursToLaJDict object not found in specified directory\n Creating from scratch...')
    if oursToLaJDict == None:
        oursToLaJDict = createOursToLaJDict()
        writeByteObjectToDefinedDirectory(args.saveLocation + r'\serialized objects', 'oursToLaJDict', oursToLaJDict)



    # Make a dict where key is the ITS2 sequence name and the value is the sequence
    global seqNameToCladeDict
    global seqToFFPProbDistDict
    seqNameToCladeDict = None
    seqToFFPProbDistDict = None
    if not args.createSeqToCladeDictFromScratch or args.developingMode:
        try:
            tempResult = readByteObjectFromDefinedDirectory(args.saveLocation + r'\serialized objects',
                                                               'seqNameToCladeDict')
            seqNameToCladeDict = tempResult[0]
            seqToFFPProbDistDict = tempResult[1]
        except:
            print('Missing Object: seqNameToCladeDict  not found in specified directory\n Creating from scratch...')
    if seqNameToCladeDict == None:
        tempResult = createSeqNameToCladeDict(masterFastaDict, 3)
        seqNameToCladeDict = tempResult[0]
        seqToFFPProbDistDict = tempResult[1]
        writeByteObjectToDefinedDirectory(args.saveLocation + r'\serialized objects', 'seqNameToCladeDict', tempResult)
    print('Completed assigning clade to sequences')

    # Create AbundanceList
    global abundanceList
    abundanceList = None
    if not args.createAbundanceListFromScratch or args.developingMode:
        try:
            abundanceList = readByteObjectFromDefinedDirectory(args.saveLocation + r'\serialized objects', 'abundanceList')
        except:
            print('Missing Object: abundanceList not found in specified directory\n Creating from scratch...')
    if abundanceList == None:
        abundanceList = {SAMPLE.name: SAMPLE for SAMPLE in createAbundnanceListMainMPRaw()}
        # abundanceList = sorted(abundanceList, key=lambda x: x.name)
        writeByteObjectToDefinedDirectory(args.saveLocation + r'\serialized objects','abundanceList', abundanceList)


    global typeDB
    typeDB = None

    print('Completed config init')

def createOursToLaJDict():
    LaJFasta = readDefinedFileToList(args.inputLocation + '/LaJeunesse Types.fas')
    its2Fasta = readDefinedFileToList(args.inputLocation + '/ITS2Fasta.fasta')

    # Remove all gaps or tail dashes
    LaJFasta = [a.replace('-','') if a[0] != '>' else a for a in LaJFasta ]
    its2Fasta = [a.replace('-','') if a[0] != '>' else a for a in its2Fasta]
    # Now find out which of our its2Fasta seqs relate to LaJ's seqs
    associatedseqs = []
    i = 0
    while i < len(LaJFasta): # For each of the LaJ seqs
        print(LaJFasta[i])
        tempList = []
        j = 0
        while j < len(its2Fasta): # Compare to our seqs
            if LaJFasta[i+1] in its2Fasta[j+1] or its2Fasta[j+1] in LaJFasta[i+1]: # If either of the sequences can be found in each other or are identical
                tempList.append(its2Fasta[j][1:])
            j += 2
        if len(tempList) > 0: # if there is something in the tempList i.e. if we have found some associates
            tempList.insert(0, LaJFasta[i].split('_')[0][1:])
            associatedseqs.append(tempList)
        i += 2

    # Create the dictionary that links our seqs to the LaJ intras
    ourstolajdict = {}
    oursToLaJList = associatedseqs
    j=0
    while j < len(oursToLaJList): # For each laJ intra cycle through the associated ours intras and put into dictionary
        k=0
        for ours in oursToLaJList[j][1:]: # For each of ours
            ourstolajdict[ours] = oursToLaJList[j][0]
            k += 1
        j += 1

    return ourstolajdict

def createMasterFastaDict():
    its2fasta = readDefinedFileToList(args.inputLocation + '/ITS2Fasta.fasta')
    return createDictFromFasta(its2fasta)


# CLADAL ASSIGNMENT
def createSeqNameToCladeDict(masterfastadict, kmerlen):
    print('Beginning cladal assignment')
    # First create the dictionary that will hold the frequency of the given features e.g. AAA : 0
    # We will count exactly the same features for each seq so rather than generating the same empty dict each time
    # we initialize a dictionary here and we will pass a copy of this on to each of the sequence characterisations


    # This will be a collection of the seq feature profiles with the seq names as the key and the feature dictionary as the value
    seqProfilesDict = {}

    # Create feature profiles of each of the sequences and stor them in the seqProfilesDict
    print('Creating Feature Frequency Profiles of submitted samples')
    for sequence in masterfastadict.keys():
        seqProfilesDict[sequence] = createFeatureProfileDict(kmerlen, masterfastadict[sequence], {''.join(list(feature)): 0 for feature in itertools.product('ACTG', repeat=kmerlen)})
    # Add in the reference sequences like any other samples.
    # We will use these ref seqs to allocate clades
    # Firstly create the reference collection of sequences
    # This is a dict where value is the ref sequence name e.g. 'refA' or 'refG' and the value is the sequence
    referenceSeqDict = createRefSeqsCollection()
    # Now create seqs for the ref seqs and add to the seqProfilesDict
    print('Creating Feature Frequency Profiles of reference samples')
    for sequence in referenceSeqDict.keys():
        seqProfilesDict[sequence] = createFeatureProfileDict(kmerlen, referenceSeqDict[sequence], {''.join(list(feature)): 0 for feature in itertools.product('ACTG', repeat=kmerlen)})
    print('Feature Frequency Profiles COMPLETE')
    # Here we create the count table
    # e.g.      AAA, AAG, AAC ...
    # SAMPLE1    0  , 49 , 2 ...
    # SAMPLE2    19 , 8,  0 ...
    # ...
    # so as to maintain order of the states when looking up the counts in the dictionary we will solidify the order of
    # the dictionary keys by making it a list and use this for our orders of the columns and rows
    listOfSeqs = list(seqProfilesDict.keys())
    listOfFeatures = [''.join(list(feature)) for feature in itertools.product('ACTG', repeat=kmerlen)]

    columnFormatedCounts2DList = convertFeatureProfileCollectionToColFormated(seqprofilesdict=seqProfilesDict, listofseqs = listOfSeqs, listoffeatures = listOfFeatures)

    # Here we normalise the counts to the total counts for each given sequence so that we end up with esentially a probability distribution

    seqToProbDistFFPDict = convertFFPCountToProbs(columnFormatedCounts2DList, listOfSeqs)

    # Here we put use the probability distributions to compare our sequences with a set of reference seqeunces
    # We assign clade according to which reference sequence we are closest to.
    print('Calculating Jensen-Shannon divergences and assigning clades')
    cladalDict = FFPProbDistToAsignedClades(seqToProbDistFFPDict)
    print('Cladal assignment COMPLETE')
    return cladalDict, seqToProbDistFFPDict

def createFeatureProfileDict(kMerLen, seq, featuredict):


    # Next work through the seq populating the dict
    # Lets try to make this dynamic so that we can be flexible with the kMerLen
    for i in range(len(seq)-(kMerLen-1)):
        # Dynamically create the feature Dict Key from the given position in the sequence we are at and the kmer length
        featureDictKey = ''

        for k in range(kMerLen):
            featureDictKey = featureDictKey + seq[i + k]

        # Then populate the featureDict with the new key and value
        featuredict[featureDictKey] = featuredict[featureDictKey] + 1


    return featuredict

def createRefSeqsCollection():
    referenceSeqDict = {}
    referenceSeqDict[
        'A'] = 'AACCAATGGCCTCTTGAACGTGCATTGCGCTCTTGGGATATGCCTGAGAGCATGTCTGCTTCAGTGCTTCTACTTTCATTTTCTGCTGCTCTTGTTATC' \
                       'AGGAGCAGTGTTGCTGCATGCTTCTGCAAGTGGCACTGGCATGCTAAATATCAAGTTTTGCTTGCTGTTGTGACTGATCAACATCTCATGTCGTTTCAGT' \
                       'TGGCGAAACAAAAGCTCATGTGTGTTCTTAACACTTCCTAGCATGAAGTCAGACAAGTGAACCCCAGACAAGTGA'
    referenceSeqDict[
        'C'] = 'AACCAATGGCCTCCTGAACGTGCGTTGCACTCTTGGGATTTCCTGAGAGTATGTCTGCTTCAGTGCTTAACTTGCCCCAACTTTGCAAGCAGGATGTGTT' \
                       'TCTGCCTTGCGTTCTTATGAGCTATTGCCCTCTGAGCCAATGGCTTGTTAATTGCTTGGTTCTTGCAAAATGCTTTGCGCGCTGTTATTCAAGTTTCTAC' \
                       'CTTCGTGGTTTTACTTGAGTGACGTTGCTCATGCTTGCAACCGCTGGGATGCAGGTGCATGCCTCTAGCATGAAGTCAGACAAGTGA'
    referenceSeqDict[
        'D'] = 'AACCAATGGCCCCGTGAACGCGCATTGCACTCTTGGGACTTCCTGAGAGTATGTTTGCTTCAGTGCTTATTTTACCTCCTTGCAAGGTTCTGTCGCAACC' \
                       'TTGTGCCCTGGCCAGCCACGGGTTAACTTGCCCATGGCTTGCTGAGTAGTGATCTTTTAGAGCAAGCTCTGGCACGCTGTTGTTTGAGGCAGCCTATATTG' \
                       'AGGCTATTTCAAATGACGTTGCTACAAGCTTGATGTGTCCTTCTGCGCCGTTGCGCATCCCATAGCATGAAGTCAAACAAGAGA'
    referenceSeqDict[
        'F'] = 'AATCAATGGCCTCCTGAACGTACGTTGCACTCTTGGGGTTTCCTGAGAGTATGTCTGCTTCAGTGCTTAGCATGCATAACCCTGCGAGCAGTTTTGTTTG' \
                       'CTTTGCGCTTTTATGAGCCATTGGTTTCCAGCCAATGGCTTGTTAAATAGTTTTTTGCAAATGAAAGCTCTGCGCGCTGTTGTTCAAGCAAGTGCCTTTC' \
                       'AGGTTTCTAGGCTTGAGTGACGCTGCTCATGCTTGCAACTGCCAGGCTGCCAGTGCACGCCTCTAGCATGAAGTCAGGCAAGTGA'
    referenceSeqDict[
        'G'] = 'AACCAATGGCCTCCTGAACGCGCATTGCACTCTTGGGCTTCCCTGAGAGTATGTTTGCTTCAGTGCTTCTTTTGCTCAACCATTGCAAGGTTTGGCAGTG' \
                       'CAATGCCTCCCTGTGCCTCGGCGTGTTGTTGGCGTGTCTGCCAATGACGTGCGACCAGCGTGGCCTATGTGCAAGCATGCACGTGCTTTGTTGTTTCACT' \
                       'GCAGACATTCTCCGGAATATGCGTGGGCGACGTGGCTGATGCTTGCGGACGCGCTAGTGTGCTGCTTGCACTTCTTCCATAGCATGAAGTCAAACAGGCA'
    return referenceSeqDict

def convertFeatureProfileCollectionToColFormated(seqprofilesdict, listofseqs, listoffeatures):
    # we are esentially making a count table here with the features across the top and seqs on the side
    # e.g.      AAA, AAG, AAC ...
    # SEQ1    0  , 49 , 2 ...
    # SEQ2    19 , 8,  0 ...
    # ...

    # This is the table in the form of a 2D list that we will populate with the feature counts
    listOfSeqProfiles = [[] for i in range(len(listofseqs))]

    # Now to iterate through the seqs' counts that are held in the seqprofsdict
    # and populate the table with them
    # I will do each of the states one at a time
    # i.e. I will cycle through states and within each state I will cycle through the seqs

    for feature in range(len(listoffeatures)): # For Each feature
        for seq in range(len(listOfSeqProfiles)): # For each sequence
            # Here we add the the correct counts to the table by looking them up in the seqprofiledict
            listOfSeqProfiles[seq].append(seqprofilesdict[listofseqs[seq]][listoffeatures[feature]])

    return listOfSeqProfiles

def convertFFPCountToProbs(colformatedcountstable, listofseqs):

    for i in range(len(colformatedcountstable)): # For every seqence or list
        # Work out the sum of all counts for each given sequence
        seqTotCount = sum(colformatedcountstable[i])

        for k in range(len(colformatedcountstable[0])): # For each item in the list i.e. each individual feature count
            colformatedcountstable[i][k] = colformatedcountstable[i][k]/seqTotCount

    FFPProbDistDict = {}
    for i in range(len(colformatedcountstable)):
        FFPProbDistDict[listofseqs[i]] = colformatedcountstable[i]


    return FFPProbDistDict

def FFPProbDistToAsignedClades(probdistdict):

    probdict = probdistdict
    cladalDict = {}
    refList = ['A', 'C', 'D', 'F','G']
    listOfKeys = probdict.keys()

    for unkSeq in listOfKeys:

        lowestCladeDivergence = 1

        for refSeq in refList:
            jsdScore = JSD(probdict[unkSeq], probdict[refSeq])
            # This was a good idea but the problem is that the 0.0165 value will change with each different data set,
            # i.e. if the sequences we are looking at are much longer than the reference sequences then this may lead
            # to smaller scores over all
            # However, it is difficult to think of a way that the minimum loosing score could become smaller.
            # This would mean that the loosing score was getting more similar to a reference seq
            # For the time being I will leave this in play but we should come back to look at a scenario in which
            # we have two clades of sequences that are quite similar e.g. clade H and clade C
            # In this scenario the smallest loosing score may get smaller.
            if jsdScore < 0.01: # I have done an empirical determinant here. I have looked at all of the jsd scores that were returned from non-correct unkseq to refseq matches.
                # The mimimum score found was 0.0165. This means that no score below this was a looser. So any score below this has to be a winner.
                cladalDict[unkSeq] = refSeq
                break
            if jsdScore < lowestCladeDivergence:
                cladalDict[unkSeq] = refSeq
                lowestCladeDivergence = jsdScore
    print('Assignment of clades through closest reference association complete')

    return cladalDict

def JSD(P, Q): # Another version of Jensen-shannon divergence making use of the entropy calculation from the scipy package
    xp = np.array(P)
    xq = np.array(Q)
    _M = 0.5 * (xp + xq)
    return 0.5 * (entropy(xp, _M) + entropy(xq, _M))

def readByteObjectFromDefinedDirectory(directory, objectname):
    f = open(directory + '\\' + objectname, 'rb')
    return pickle.load(f)

def createAbundnanceListMainMP():
    print('Starting createAbudnanceListMainMP()')
    # Create queues
    task_queue = Queue()
    done_queue = Queue()

    resultList = []

    #Put chunks for processing into task_queue
    for chunk in createChunkOfAbunanceData(args.numProcessors): #This produces n[numProcessors] subsets of the abundanceData
        task_queue.put(chunk)

    print('Starting %s subprocesses' % args.numProcessors)
    all_processes = []
    for n in range(args.numProcessors):
        all_processes.append(Process(target=worker, args=(task_queue, done_queue, args), daemon=True))
        task_queue.put('STOP')

    for p in all_processes:
        p.start()


    for i in range(args.numProcessors):
        resultList.extend(done_queue.get())

    print('Finished createAbudnanceListMainMP()')
    return list(resultList)

def createAbundnanceListMainMPRaw():
    print('Starting createAbudnanceListMainMP()')
    # Create queues
    task_queue = Queue()
    done_queue = Queue()

    resultList = []

    #Put chunks for processing into task_queue
    for chunk in createChunkOfAbunanceDataRaw(args.numProcessors): #This produces n[numProcessors] subsets of the abundanceData
        task_queue.put(chunk)

    print('Starting %s subprocesses' % args.numProcessors)
    all_processes = []
    for n in range(args.numProcessors):
        # all_processes.append(Process(target=workerRaw, args=(task_queue, done_queue, args), daemon=True))
        all_processes.append(Process(target=workerRaw, args=(task_queue, done_queue, args, seqNameToCladeDict), daemon=True))
        task_queue.put('STOP')



    for p in all_processes:
        p.start()


    for i in range(args.numProcessors):
        resultList.extend(done_queue.get())

    print('Finished createAbudnanceListMainMP()')
    return list(resultList)

def writeByteObjectToDefinedDirectory(directory, objectString, object):
    try:
        os.makedirs(directory)
    except FileExistsError:
        pass
    f = open(directory + '\\' + objectString, 'wb+')
    pickle.dump(object, f)

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

def readDefinedFileToList(filename):
    tempList = []
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

# Create dictionary from Fasta
def createDictFromFasta(fastaList):
    tempDict = {}
    i = 0
    while i < len(fastaList):
        tempDict[fastaList[i][1:]] = fastaList[i+1]
        i += 2
    return tempDict

def createChunkOfAbunanceData(n):
    chunkList = []
    abundanceData = readDefinedFileTo2DList(args.inputLocation + '/ITS2Abundance.txt')
    abundanceDataInfoOnly = [a[1:len(abundanceData[0])-5] for a in abundanceData]
    roundDown = int(len(abundanceDataInfoOnly[0])/n)
    for i in range(n):
        if i == n-1: #If this is the last one then add remainder otherwise add roundedDown(len(abundanceDataInfoOnly)/n) of abundanceDataInfoOnly to eachlist
            tempList = [[a[0]] for a in abundanceData] # This creates a list of lists with the first item already in
            j = 0
            while j < len(tempList): #For each of the lists i.e. its2 seqs and info
                tempList[j].extend(abundanceDataInfoOnly[j][:]) #Add the sample info from the abundanceDataInfoOnlylist
                if j not in [1, 2, 3, 4, 5]:
                    tempList[j].append(abundanceData[j][-4])
                j += 1

        else:
            tempList = [[a[0]] for a in abundanceData] # This creates a list of lists with the first item already in
            j = 0
            while j < len(tempList): #For each of the lists i.e. its2 seqs and info
                tempList[j].extend(abundanceDataInfoOnly[j][:roundDown]) #Add the sample info from the abundanceDataInfoOnlylist
                if j not in [1, 2, 3, 4, 5]:
                    tempList[j].append(abundanceData[j][-4])
                del abundanceDataInfoOnly[j][:roundDown] #Delete the added data from the abundanceDataInfoOnlylist
                j += 1

        print('Processed abundance data chunk ' + str(i))
        chunkList.append(tempList)

    return chunkList

def createChunkOfAbunanceDataRaw(n):
    abundanceData = readDefinedFileTo2DList(args.inputLocation + '/ITS2AbundanceRaw.txt')

    chunkList = []
    # abundanceData = readDefinedFileTo2DList(args.inputLocation + '/ITS2Abundance.txt')
    # abundanceDataInfoOnly = [a[1:len(abundanceData[0])-5] for a in abundanceData]
    roundDown = int(len(abundanceData[0])/n)
    for i in range(n):

        # If this is the last chunk to be made
        if i == n-1: # If this is the last one then add remainder otherwise add roundedDown(len(abundanceDataInfoOnly)/n) of abundanceDataInfoOnly to eachlist
            tempList = [[a[0]] for a in abundanceData] # This creates a list of lists with the first item already in
            j = 0
            while j < len(tempList): # For each of the lists i.e. its2 seqs and info
                tempList[j].extend(abundanceData[j][1:]) # Add the sample info from the abundanceDataInfoOnlylist
                j += 1

        # If this is not the last chunk to be made
        else:
            # if this is not the last chunk we are dealing with we cut out runDown number of samples from the abundance list
            # Create the new empty table that will be the new abundancedata subset table.
            # each list is an ITS2 sequence so the name of that sequence is put in as the first item of each list
            tempList = [[a[0]] for a in abundanceData]
            j = 0
            while j < len(tempList): # For each of the lists i.e. its2 seqs and info
                tempList[j].extend(abundanceData[j][1:roundDown]) # Add the sample info from the abundanceData
                del abundanceData[j][1:roundDown] # Delete the added data from abundancelist. This will shift all of the data to the left. # This way we don't have to change the roundDown number
                j += 1

        print('Processed abundance data chunk ' + str(i))
        chunkList.append(tempList)

    return chunkList

def readDefinedFileTo2DList(filename):
    tempListPrimary = []
    tempListSecondary = []

    with open(filename, mode='r') as reader:
        for line in reader:
            tempListSecondary = line.rstrip().split(',')
            tempListPrimary.append(tempListSecondary)
    return tempListPrimary

def worker(input, output, args): # Input = task_queue Output = out_Queue
    sys.stdout = open(args.saveLocation + '/subprocess-' + str(os.getpid()) + '.out', 'w')
    sys.stderr = sys.stdout
    print('worker starting with process id: ' + str(os.getpid()))
    sys.stdout.flush()
    sys.stderr.flush()
    its2fasta = readDefinedFileToList(args.inputLocation + '/ITS2Fasta.fasta')
    masterFastaDict = createDictFromFasta(its2fasta)
    for chunk in iter(input.get, 'STOP'): # For each chunk in the input queue
        partialAbundanceList = createAbundanceListMultiProcessClasses(chunk, masterFastaDict) #Pass the chunk into the createAbundanceListMultiProcess function and retrieve output
        output.put(partialAbundanceList) # Put the part of the newAbundanceList into the output queue
        # output.append(partialAbundanceList)
    print('Worker finished')
    sys.stdout.flush()
    sys.stderr.flush()

def workerRaw(input, output, args, cladedict):  # Input = task_queue Output = out_Queue
    sys.stdout = open(args.saveLocation + '/subprocess-' + str(os.getpid()) + '.out', 'w')
    sys.stderr = sys.stdout
    print('worker starting with process id: ' + str(os.getpid()))
    sys.stdout.flush()
    sys.stderr.flush()
    its2fasta = readDefinedFileToList(args.inputLocation + '/ITS2Fasta.fasta')
    masterFastaDict = createDictFromFasta(its2fasta)
    for chunk in iter(input.get, 'STOP'):  # For each chunk in the input queue
        partialAbundanceList = createAbundanceListMultiProcessClassesRaw(chunk, masterFastaDict, cladedict)  # Pass the chunk into the createAbundanceListMultiProcess function and retrieve output
        output.put(partialAbundanceList)  # Put the part of the newAbundanceList into the output queue
        # output.append(partialAbundanceList)
    print('Worker finished')
    sys.stdout.flush()
    sys.stderr.flush()

def createAbundanceListRaw():
    abundanceData = readDefinedFileTo2DList(args.inputLocation + '/ITS2AbundanceRaw.txt')
    its2fasta = readDefinedFileToList(args.inputLocation + '/ITS2Fasta.fasta')
    masterFastaDict = createDictFromFasta(its2fasta)
    SeqNameToCladeDict = createSeqNameToCladeDict(masterFastaDict, 3)
    # This esentially gives us the abundanceData
    # I think it would be a good idea to try to make this compatible with the current work flow so that we can
    # keep the multi-threading









def createAbundanceListMultiProcessClasses(chunk, fastaDict):
    col = 1
    outPut = []
    while col < len(chunk[0][:-1]):# For each sample i.e. col not including the final clade column
        its2SequenceOccuranceList = []
        row = 6

        while row < len(chunk):# For each ITS2 sequence occurance
            if int(chunk[row][col]) != 0:#Then this sample contains this seq
                # Create instance of its2SequenceOccurance
                tempIts2seqOccurance = its2SequenceOccurance(name = chunk[row][0], clade = chunk[row][-1], abundance = int(chunk[row][col]), sequence= fastaDict[chunk[row][0]] )
                # Start to create list of its2SequenceOccurances
                its2SequenceOccuranceList.append(tempIts2seqOccurance)
            row += 1
        # Sort the its2OccuranceList according to their abundance
        sortedIts2SequenceOccuranceList = sorted(its2SequenceOccuranceList, key=lambda x: x.abundance, reverse=True)
        tempCompleteComplement = completeComplement(sortedIts2SequenceOccuranceList)
        outPut.append(sample(name=chunk[0][col], compComplement=tempCompleteComplement, hostTaxon=chunk[1][col], region=chunk[2][col], reef=chunk[3][col], totalSeqs=int(chunk[5][col])))
        col += 1
        print(multiprocessing.current_process().name + ' ' + str(col))
    return outPut

def createAbundanceListMultiProcessClassesRaw(chunk, fastaDict, seqnametocladedict):

    # Here, we go coral sample by coral sample making a note of any ITS2 sequences that they contain
    # When we find an ITS2 sequence it contains we add it create an ITS2SeqOccurance and add it to a list
    # which contains other occurances
    # Once we have completed the counts for one coral we then create an instance of the sample class
    # Because we will read in other sample parameters seperate at another data, we will not add that to the
    # sample instance at this point

    outPut = []
    # col = coral sample
    # row = its2 sequence
    for col in range(1, len(chunk[0])): # For each of the coral samples; we start at one to avoid the sample title
        its2SequenceOccuranceList = []
        totSeqs = 0

        for row in range(1, len(chunk)):# For each ITS2 sequence occurance; we start at 1 to avoid the sample title
            if int(chunk[row][col]) != 0:#Then this sample contains this seq
                # Create instance of its2SequenceOccurance
                tempIts2seqOccurance = its2SequenceOccurance(name = chunk[row][0], clade = seqnametocladedict[chunk[row][0]], abundance = int(chunk[row][col]), sequence= fastaDict[chunk[row][0]] )
                # Start to create list of its2SequenceOccurances
                its2SequenceOccuranceList.append(tempIts2seqOccurance)
                # Keep track of how many reads were taken in this sequence
                totSeqs += int(chunk[row][col])

        # Sort the its2OccuranceList according to their abundance
        sortedIts2SequenceOccuranceList = sorted(its2SequenceOccuranceList, key=lambda x: x.abundance, reverse=True)
        tempCompleteComplement = completeComplement(sortedIts2SequenceOccuranceList)
        # We need to calculate the total seqs at some point
        outPut.append(sample(name=chunk[0][col], compComplement=tempCompleteComplement, hostTaxon=None, region=None, reef=None, totalSeqs=totSeqs))

        print(multiprocessing.current_process().name + ' ' + str(col))
    return outPut

def checkRLibraries(rinstalledcorrectly):
    listOfLibs = ['amap', 'bios2mds', 'cluster', 'colorspace', 'dichromat',
                  'e1071', 'labeling', 'munsell', 'plyr',
                  'RColorBrewer', 'Rcpp', 'rgl', 'scales', 'qcc']
    for module in listOfLibs:
        if not os.path.exists(args.rootLocation + '/Rlibs/' + module):
            rinstalledcorrectly = False
            return rinstalledcorrectly
    return rinstalledcorrectly
