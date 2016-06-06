import argparse
from tkinter import messagebox
import pickle
from multiprocessing import Queue, Process
import multiprocessing
import os
from subprocess import call
import sys

from ClassesTwo import *

def __init__():
    parser = argparse.ArgumentParser()
    
    # Computation parameters
    parser.add_argument('--cutOff', type=float, default=0.1, metavar='CUTOFF')
    parser.add_argument('--typeSupport', type=int, default=3, help='The number of samples that must contain a type in order for it to be considered as a genuine type', metavar='N')
    parser.add_argument('--coDomSupport', type=int, default=2, help='The number of samples that must contain the an alternative MajITS2 sequence compared to the other samples of the type in order for it be considered a codom type', metavar='N')
    parser.add_argument('--cladeCollectionCutoff', type=float, default=0.1, metavar='CUTOFF')
    parser.add_argument('--cladeList', nargs='*', default=['A', 'B', 'C', 'D'], choices=['A', 'B', 'C', 'D'])
    parser.add_argument('--numProcessors', type=int, default=3, metavar='N')
    parser.add_argument('--majPlotThreshold', type=float, default=2, help='The number of supported types required within a Maj in order for its plotting data to be made', metavar='N')
    
    # Paths
    cwd = os.path.dirname(__file__)
    parser.add_argument('--rootLocation', default=cwd, help='Directory where the source code is found', metavar='PATH')
    parser.add_argument('--inputLocation', default=cwd + '/raw data', help='Directory where the three input files are found', metavar='PATH')
    parser.add_argument('--saveLocation', default=cwd, help='Output directory for saving matrices and output tables', metavar='PATH')
    parser.add_argument('--mothurLocation', default=cwd + '/Mothur/mothur.exe', help='Full path including mothur.exe', metavar='PATH')
    parser.add_argument('--rscriptLocation', default=cwd + '/R/R-3.3.0/bin/x64/Rscript.exe', help='Full path including Rscript.exe', metavar='PATH')
    parser.add_argument('--logToFile', action='store_true', help='Redirect stdout/stderr to files in saveLocation')
    
    # Caching
    parser.add_argument('--createAbundanceListFromScratch', type=bool, default=False, metavar='TRUE|FALSE')
    parser.add_argument('--createMasterSeqDistancesFromScratch', type=bool, default=False, metavar='TRUE|FALSE')
    parser.add_argument('--createFstColDistsFromScratch', type=bool, default=False, metavar='TRUE|FALSE')
    parser.add_argument('--createOursToLaJDictFromScratch', type=bool, default=False, metavar='TRUE|FALSE')
    
    global args
    args = parser.parse_args()
    
    # Before we do anything else, redirect stdout/sterr if requested
    # http://stackoverflow.com/a/4110906/1688738
    if args.logToFile:
        sys.stdout = open(args.saveLocation + '/stdout.txt', 'w');
        sys.stderr = open(args.saveLocation + '/stderr.txt', 'w');
    
    global masterFastaDict
    masterFastaDict = createMasterFastaDict()
    global reinstated
    reinstated = True

    print('Checking for R libraries...')
    RInstalledCorrectly = True
    RInstalledCorrectly = checkRLibraries(RInstalledCorrectly)

    # Print when sucessful
    # Go through the directory and check to see if any of the libraries are missing. If any are missing then RInstalledCorrectly = False
    # call([RscriptLocation, rootLocation + r"\.bat scripts\install.bat"], shell=True)
    if RInstalledCorrectly:
        print('Groovy, R libraries have been verified')
    else: # not RInstalledCorrectly:
        print('Ooops, R libraries missing. Installing...')
        call([args.rscriptLocation, args.rootLocation + r"\.bat scripts\installRLibs.bat", args.rootLocation.replace('\\', '/')], shell=True)

    global oursToLaJDict
    oursToLaJDict = None
    if not args.createOursToLaJDictFromScratch:
        try:
            oursToLaJDict = readByteObjectFromDefinedDirectory(args.saveLocation + r'\serialized objects', 'oursToLaJDict')
        except:
            messagebox.showwarning('Missing Object', 'oursToLaJDict object not found in specified directory\n Creating from scratch...')
    if oursToLaJDict == None:
        oursToLaJDict = createOursToLaJDict()
        writeByteObjectToDefinedDirectory(args.saveLocation + r'\serialized objects', 'oursToLaJDict', oursToLaJDict)

    # Create AbundanceList
    global abundanceList
    abundanceList = None
    if not args.createAbundanceListFromScratch:
        try:
            abundanceList = readByteObjectFromDefinedDirectory(args.saveLocation + r'\serialized objects', 'abundanceList')
        except:
            messagebox.showwarning('Missing Object', 'abundanceList object not found in specified directory\n Creating from scratch...')
    if abundanceList == None:
        abundanceList = createAbudnanceListMainMP() # This is the multiProcessed creation of the abundanceList
        abundanceList = sorted(abundanceList, key=lambda x: x.name)
        writeByteObjectToDefinedDirectory(args.saveLocation + r'\serialized objects','abundanceList', abundanceList)
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

def readByteObjectFromDefinedDirectory(directory, objectname):
    f = open(directory + '\\' + objectname, 'rb')
    return pickle.load(f)

def createAbudnanceListMainMP():
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
    # for p in all_processes:
    #     p.join()

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
                # print('Chunking ' + str(i) + '\t' + str(j))
        else:
            tempList = [[a[0]] for a in abundanceData] # This creates a list of lists with the first item already in
            j = 0
            while j < len(tempList): #For each of the lists i.e. its2 seqs and info
                tempList[j].extend(abundanceDataInfoOnly[j][:roundDown]) #Add the sample info from the abundanceDataInfoOnlylist
                if j not in [1, 2, 3, 4, 5]:
                    tempList[j].append(abundanceData[j][-4])
                del abundanceDataInfoOnly[j][:roundDown] #Delete the added data from the abundanceDataInfoOnlylist
                j += 1
                # print('Chunking ' + str(i) + '\t' + str(j))
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

def checkRLibraries(rinstalledcorrectly):
    listOfLibs = ['amap', 'bios2mds', 'cluster', 'colorspace', 'dichromat',
                  'e1071', 'labeling', 'munsell', 'plyr',
                  'RColorBrewer', 'Rcpp', 'rgl', 'scales', 'qcc']
    for module in listOfLibs:
        if not os.path.exists(args.rootLocation + '/Rlibs/' + module):
            rinstalledcorrectly = False
            return rinstalledcorrectly
    return rinstalledcorrectly
