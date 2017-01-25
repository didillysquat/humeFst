import os
import subprocess

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


# listOfDirs = [
#     r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/AG1199/',
#     r'/home/benjamin/fstProjectWorking/fastQ/pythonFastQ/pythonArchive/AG1198/'
# ]

def MEDanlysis(listOfDirs):
    #By creating a list of these paths, we can run all of the MED calls in parallel
    #This should most efficient on memory
    fastaFilePathList = []
    # This is the form that we will send the processed data back to the View in
    # FORM OF:
    # (SAMPLENAME, [(CLADE, pathToMED.fasta),(CLADE, pathToMED.fasta)]
    dataPacket = []
    for givenDir in listOfDirs:
        sampleName = givenDir.split('/')[-2]
        sampleTup = (sampleName, [])

        # list the folder names contained in the directory in question
        for directory in next(os.walk(givenDir))[1]:
            clade = directory.replace('clade','')
            pathToDir = '{0}{1}/'.format(givenDir, directory)
            # For each of the files in each of the Cladal directories
            fastaFilePath = ''
            for files in next(os.walk(pathToDir))[2]:
                if '.redundant.fasta' in files:
                    # This is the fasta file
                    fastaFilePath = '{0}{1}'.format(pathToDir, files)
                    fastaFilePathList.append(fastaFilePath)
                    sampleTup[1].append((clade,'{0}{1}'.format(fastaFilePath,'-PADDED-WITH-GAPS')))
                    dataPacket.append(sampleTup)



    processes = [subprocess.Popen(['o-pad-with-gaps', fastaFile]) for fastaFile in fastaFilePathList]
    for p in processes:
        p.wait()
    processes = [subprocess.Popen(['decompose', r'{0}{1}'.format(fastaFile, '-PADDED-WITH-GAPS'), '--skip-gen-html', '--skip-gen-figures', '--skip-gexf-files', '--skip-check-input-file', '-o', '{0}/'.format(os.path.dirname(fastaFile)) ]) for fastaFile in fastaFilePathList]
    for p in processes:
        p.wait()

    return dataPacket


# MEDanlysis(listOfDirs)