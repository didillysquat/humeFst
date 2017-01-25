import os
import subprocess

# from separateClade import splitClades
# from MEDanalysis import MEDanlysis

# This file can take a file from one of the views. This file will be a .gz file
# For the time being we will give ourselvs the file path

def moveFile(pathfrom, pathto):
    print('moving list to ' + pathto)
    try:
        os.makedirs(os.path.dirname(pathto))
    except FileExistsError:
        pass

    os.rename(pathfrom, pathto)

def writeListToDestination(destination, listToWrite):
    #print('Writing list to ' + destination)
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



def main():
    pathToFile = r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/Archive.zip'
    unzipDestination = r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/'
    primerFwd = 'GAATTGCAGAACTCCGTG' # Written 5'-->3'
    primerRev = 'GGATCCATATGCTTAAGTTCAGCGGGT' # Written 5'-->3'
    oligoFile = [
        r'#ITSintfor2',
        'forward\t{0}'.format(primerFwd),
        r'#ITS2reverse',
        'reverse\t{0}'.format(primerRev)
    ]
    # moveFile() if necessary

    # Unzipfile inplace
    # unzip file.zip -d destination_folder

    #---------------#completedProcess = subprocess.run(["unzip", pathToFile, '-d', unzipDestination])

    # Get rid of any dashes in the .gz filenames
    #---------------#for file in os.listdir(unzipDestination):
    #---------------#    if file.endswith(".gz"):
    #---------------#        if file.find("-") > 0:
    #---------------#            os.rename(r'{0}{1}'.format(unzipDestination, file), r'{0}{1}'.format(unzipDestination, file.replace('-','')))

    a = 'thisis'

    #Make a batch file for mothur, set input and output dir and create a file file
    # mBatchFile = [
    #     r'set.dir(input=/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/)',
    #     r'set.dir(output=/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/)',
    #     r'make.file(inputdir=/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/, type=gz, numcols=3)'
    #                ]
    # writeListToDestination(r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonMSIDE/mBatchFile', mBatchFile)
    # completedProcess = subprocess.run([r'mothur', r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonMSIDE/mBatchFile' ])

    # Read in the stability.files
    # Use this to make the individual contig files
    sampleFastQPairs = readDefinedFileToList(r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/stability.files')
    # This is the main control loop that will go through each of the .gz pairs
    print('Processing fastQ files')
    #As running the mothur pipeline is going to be one of the slowest parts of dataInput
    # and running through the python wrapper doesn't seem to invoke multiple processors
    # I think it would be wise to try to run the sub.processes in parallel.
    # To do this we need to set up the environment first
    # Then call all the subprocesses in a list.
    # TODO I have not debugged or tested the multi-process component of this
    # as I have already completed a full run through and so don't want to wait to generate all of the files again
    mBatchFilePathList = []
    directoryList = []
    for contigPair in sampleFastQPairs:
        sampleName = contigPair.split('\t')[0]
        print('Sample: {0}'.format(sampleName))
        currentDir = r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/{0}/'.format(sampleName)
        os.makedirs(currentDir, exist_ok=True)
        stabilityFile = [contigPair]
        stabilityFileName = r'{0}{1}'.format(sampleName,'stability.files')
        rootName = r'{0}stability'.format(sampleName)
        stabilityFilePath = r'{0}{1}'.format(currentDir,stabilityFileName)
        writeListToDestination(stabilityFilePath, stabilityFile)
        # Write oligos file to directory
        writeListToDestination('{0}{1}'.format(currentDir, 'primers.oligos'), oligoFile)
        mBatchFile = [
            r'set.dir(input={0})'.format(currentDir),
            r'set.dir(output={0})'.format(currentDir),
            r'make.contigs(file={0}, processors=4)'.format(stabilityFileName),
            r'summary.seqs(fasta={0}.trim.contigs.fasta, processors=4)'.format(rootName),
            r'screen.seqs(fasta={0}.trim.contigs.fasta, group={0}.contigs.groups, maxambig=0, maxhomop=5, processors=4)'.format(rootName),
            r'summary.seqs(fasta={0}.trim.contigs.good.fasta, processors=4)'.format(rootName),
            r'unique.seqs(fasta={0}.trim.contigs.good.fasta)'.format(rootName),
            r'summary.seqs(fasta={0}.trim.contigs.good.unique.fasta, name={0}.trim.contigs.good.names, processors=4)'.format(rootName),
            r'split.abund(cutoff=2, fasta={0}.trim.contigs.good.unique.fasta, name={0}.trim.contigs.good.names, group={0}.contigs.good.groups)'.format(rootName),
            r'summary.seqs(fasta={0}.trim.contigs.good.unique.abund.fasta, name={0}.trim.contigs.good.abund.names, processors=4)'.format(rootName),
            r'summary.seqs(fasta={0}.trim.contigs.good.unique.rare.fasta, name={0}.trim.contigs.good.rare.names, processors=4)'.format(rootName),
            r'pcr.seqs(fasta={0}.trim.contigs.good.unique.abund.fasta, group={0}.contigs.good.abund.groups, name={0}.trim.contigs.good.abund.names, oligos=primers.oligos, pdiffs=2, rdiffs=2, processors=4)'.format(rootName),
            r'summary.seqs(fasta={0}.trim.contigs.good.unique.abund.pcr.fasta, name={0}.trim.contigs.good.abund.pcr.names, processors=4)'.format(rootName),
            r'screen.seqs(fasta={0}.trim.contigs.good.unique.abund.pcr.fasta, name={0}.trim.contigs.good.abund.pcr.names, group={0}.contigs.good.abund.pcr.groups,  minlength=239, maxlength=307, processors=4)'.format(rootName),
            r'summary.seqs(fasta={0}.trim.contigs.good.unique.abund.pcr.good.fasta, name={0}.trim.contigs.good.abund.pcr.good.names, processors=4)'.format(rootName),
        ]
        mBatchFilePath = r'{0}{1}{2}'.format(currentDir,'mBatchFile', sampleName)
        mBatchFilePathList.append(mBatchFilePath)
        directoryList.append(currentDir)
        writeListToDestination(mBatchFilePath, mBatchFile)
        # completedProcess = subprocess.run([r'mothur', r'{0}'.format(mBatchFilePath)])
    processes = [subprocess.Popen([r'mothur', r'{0}'.format(mBatchFilePath)]) for mBatchFilePath in mBatchFilePathList]
    for p in processes:
        p.wait()
    print('Complete')
        # At this point we should have the fasta, name and group files for the given sample ready to be thrown into cladal analysis.
        # It may be worth doing one final step for this which is to group all sequences from this dataset together before conducting the cladal analysis
        # To minimise the time taken to perform this.
        # Alternatively query the db with the seqs and then only send the unknown seqs onto the cladal script

    return directoryList

listOfDirs = [
    r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/AG1199/',
    r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/AG1198/'
]

def deuniqueFiles(listofdirs):
    mBatchFilePathList = []
    for givenDir in listofdirs:
        # for each cladal directory
        sampleName = givenDir.split('/')[-2]
        for directory in next(os.walk(givenDir))[1]:

            fastaFilePath = ''
            nameFilePath = ''
            groupFilePath = ''
            pathToDir = '{0}{1}'.format(givenDir, directory)
            cladeName = directory
            # For each of the files in each of the Cladal directories
            for files in next(os.walk(pathToDir))[2]:
                if '.fasta' in files:
                    fastaFilePath = '{0}/{1}'.format(pathToDir, files)
                elif '.names' in files:
                    nameFilePath = '{0}/{1}'.format(pathToDir, files)

            # Build a quick mBatchFile
            mBatchFile = [
                r'set.dir(input={0}/)'.format(pathToDir),
                r'set.dir(output={0}/)'.format(pathToDir),
                r'deunique.seqs(fasta={0}, name={1})'.format(fastaFilePath, nameFilePath)
            ]
            mBatchFilePath = '{0}/{1}'.format(pathToDir, '{0}{1}{2}'.format(sampleName, cladeName, 'mBatchFile'))
            writeListToDestination(mBatchFilePath, mBatchFile)
            mBatchFilePathList.append(mBatchFilePath)
    processes = [subprocess.Popen([r'mothur', r'{0}'.format(mBatchFilePath)]) for mBatchFilePath in mBatchFilePathList]
    for p in processes:
        p.wait()
    something = 'Apples'
#fastQDirs = main()
#splitClades(fastQDirs)
#deuniqueFiles(listOfDirs)
dataPacket = MEDanlysis(fastQDirs)
