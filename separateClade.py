from scipy.stats import entropy
import numpy as np
import itertools
import os
import argparse


parser = argparse.ArgumentParser()
cwd = os.path.dirname(__file__)
parser.add_argument('--fasta', default=None ,help='Directory where the three input files are found', metavar='PATH')
parser.add_argument('--name', default=None ,help='Directory where the three input files are found', metavar='PATH')
parser.add_argument('--group', default=None ,help='Directory where the three input files are found', metavar='PATH')
parser.add_argument('--outputDir', default=None ,help='Directory where the three input files are found', metavar='PATH')
args = parser.parse_args()

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

def readDefinedFileToList(filename):
    tempList = []
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList


def createDictFromFasta(fastaList):
    tempDict = {}
    i = 0
    while i < len(fastaList):
        tempDict[fastaList[i][1:]] = fastaList[i+1]
        i += 2
    return tempDict


def createMasterFastaDict(directory):
    its2fasta = readDefinedFileToList(directory)
    # Remove all gaps or tail dashes
    its2Fasta = [a.replace('-', '') if a[0] != '>' else a.split('\t')[0] for a in its2fasta]
    return createDictFromFasta(its2Fasta)

# CLADAL ASSIGNMENT
def createSeqNameToCladeDict(masterfastadict):
    print('Beginning cladal assignment')
    # First create the dictionary that will hold the frequency of the given features e.g. AAA : 0
    # We will count exactly the same features for each seq so rather than generating the same empty dict each time
    # we initialize a dictionary here and we will pass a copy of this on to each of the sequence characterisations


    # This will be a collection of the seq feature profiles with the seq names as the key and the feature dictionary as the value
    seqProfilesDict = {}

    # Create feature profiles of each of the sequences and store them in the seqProfilesDict
    print('Creating Feature Frequency Profiles of submitted samples')
    # each feature is e.g. ('A','A',A')
    listOfFeaturesIter = itertools.product('ACTG', repeat=3)
    featureDict = {''.join(list(feature)): 0 for feature in listOfFeaturesIter}
    for sequence in masterfastadict.keys():
        seqProfilesDict[sequence] = createFeatureProfileDict(3, masterfastadict[sequence], featureDict.copy())
    # Add in the reference sequences like any other samples.
    # We will use these ref seqs to allocate clades
    # Firstly create the reference collection of sequences
    # This is a dict where value is the ref sequence name e.g. 'refA' or 'refG' and the value is the sequence
    referenceSeqDict = createRefSeqsCollection()
    # Now create seqs for the ref seqs and add to the seqProfilesDict
    print('Creating Feature Frequency Profiles of reference samples')
    for sequence in referenceSeqDict.keys():
        seqProfilesDict[sequence] = createFeatureProfileDict(3, referenceSeqDict[sequence], featureDict.copy())
    print('Feature Frequency Profiles COMPLETE')
    # Here we create the count table
    # e.g.      AAA, AAG, AAC ...
    # SAMPLE1    0  , 49 , 2 ...
    # SAMPLE2    19 , 8,  0 ...
    # ...
    # so as to maintain order of the states when looking up the counts in the dictionary we will solidify the order of
    # the dictionary keys by making it a list and use this for our orders of the columns and rows
    listOfSeqs = list(seqProfilesDict.keys())
    listOfFeatures = [''.join(list(feature)) for feature in itertools.product('ACTG', repeat=3)]
    ## TODO I think we can streamline this, I'm not sure why we are making the column format first,
    # we should be able to go straight to the probability dictionary structure
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
    lengthOfList = len(listOfKeys)
    count = 0
    next = 0.01
    print("0% complete")
    for unkSeq in listOfKeys:
        # print('Assigning clade: {0}'.format(unkSeq))
        lowestCladeDivergence = 1

        if count/lengthOfList > next:
            print("{0}%".format(round(next,2)*100))
            next += 0.01
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
            # TODO 230117 I think we just have to be safe here and go with the smallest score
            # if jsdScore < 0.01: # I have done an empirical determinant here. I have looked at all of the jsd scores that were returned from non-correct unkseq to refseq matches.
            #     # The mimimum score found was 0.0165. This means that no score below this was a looser. So any score below this has to be a winner.
            #     cladalDict[unkSeq] = refSeq
            #     break
            if jsdScore < lowestCladeDivergence:
                cladalDict[unkSeq] = refSeq
                lowestCladeDivergence = jsdScore
        count += 1
    print('Assignment of clades through closest reference association complete')

    return cladalDict

def JSD(P, Q): # Another version of Jensen-shannon divergence making use of the entropy calculation from the scipy package
    xp = np.array(P)
    xq = np.array(Q)
    _M = 0.5 * (xp + xq)
    return 0.5 * (entropy(xp, _M) + entropy(xq, _M))


def writeCladeSeparatedFilesForMothurDeunique(cladaldict, fastadict, namefilepath, directory, samplename):
    nameFile = readDefinedFileToList(namefilepath)
    nameDict = {a.split('\t')[0]: a for a in nameFile}

    cladeList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

    # For each clade a list of lists where each set of three lists is documents for fasta, name and group
    documentList = [[[], []] for i in range(len(cladeList))]

    # We sort of discard the group list concept at this point as we know that all seqs are from the same sample
    for key in nameDict.keys():
        clade = cladaldict[key]
        listNum = cladeList.index(clade)
        # Write the fasta lines
        documentList[listNum][0].extend(['>{0}'.format(key), fastadict[key]])
        #Write the names lines
        documentList[listNum][1].append('{0}'.format(nameDict[key]))


    # Here the lists should be populated and we can now write
    for i in range(len(documentList)):
        if len(documentList[i][0]) > 4:
            # Write fasta
            writeListToDestination(r'{0}clade{1}/uniqueClade{1}{2}.fasta'.format(directory, cladeList[i], samplename), documentList[i][0])
            # Write names
            writeListToDestination(r'{0}clade{1}/uniqueClade{1}{2}.names'.format(directory, cladeList[i], samplename), documentList[i][1])




def splitClades(listOfDirectories):
    for directory in listOfDirectories:
        sampleName = directory.split('/')[-2]
        fastaFilePath = r'{0}{1}stability.trim.contigs.good.unique.abund.pcr.good.fasta'.format(directory, sampleName)
        nameFilePath =  r'{0}{1}stability.trim.contigs.good.abund.pcr.good.names'.format(directory, sampleName)

        masterFastaDict = createMasterFastaDict(fastaFilePath)
        cladalDict = createSeqNameToCladeDict(masterFastaDict)[0]
        writeCladeSeparatedFilesForMothurDeunique(cladalDict, masterFastaDict,namefilepath=nameFilePath, directory=directory, samplename=sampleName)

# listOfDirs = [
#     r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/AG1199/',
#     r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/AG1198/'
# ]

# splitClades(listOfDirs)