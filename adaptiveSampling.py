import numpy as np
import sys
import string
import shutil
import os
import subprocess
import json
#from pyproct.clustering.clustering import Clustering
import math
import atexit
import argparse
import time
import glob
import clustering
import pickle
import startingConformationsCalculator
import blockNames
import atomset

import multiprocessing
from functools import partial
import SimulationRunner


def copyInitialStructures(initialStructures, tmpInitialStructuresTemplate, iteration):
    for i, name in enumerate(initialStructures):
        shutil.copyfile(name, tmpInitialStructuresTemplate%(iteration,i))

def createMultipleComplexesFilenames(numberOfSnapshots, inputFileTemplate, tmpInitialStructuresTemplate, iteration):
    jsonString = "\n"
    for i in range(numberOfSnapshots-1):
        jsonString += inputFileTemplate%(tmpInitialStructuresTemplate%(iteration,i)) + ",\n"
    jsonString += inputFileTemplate%(tmpInitialStructuresTemplate%(iteration,numberOfSnapshots-1))
    return jsonString

def cleanup(tmpFolder):
    if os.path.exists(tmpFolder):
        shutil.rmtree(tmpFolder)

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

def generateFilteredListSnapshotSelectionString(currentEpoch, compressionPath):
    selectionString = ""
    for i in range(currentEpoch):
        selectionString += " \"" + compressionPath%i + "\", "
    selectionString += " \"" + compressionPath%currentEpoch + "\""
    return selectionString

def generateSnapshotSelectionString(currentEpoch, outputPathTempletized, trajectoryWildcard):
    selectionString = ""
    for i in range(currentEpoch):
        selectionString += " \"" + os.path.join(outputPathTempletized%i, trajectoryWildcard) + "\", "
    selectionString += " \"" + os.path.join(outputPathTempletized%currentEpoch, trajectoryWildcard) + "\""
    return selectionString

def generateSnapshotSelectionStringLastRound(currentEpoch, outputPathTempletized, trajectoryWildcard):
    return " \"" + os.path.join(outputPathTempletized%currentEpoch, trajectoryWildcard) + "\""

def countEpochAcceptedSteps(folder):
    acceptedSteps = os.popen("tail -qn 1 " + os.path.join(folder,"*report_*") + " | awk '{a+=$3}END{print a}'").read()
    return acceptedSteps

def runPyProCT(pyProctControlFile, PYTHON, PYPROCT):
    # example to run pyproct in life: /data2/apps/PYTHON/2.7.5/bin/python2.7 -m pyproct.main controlFile
    toRun = [PYTHON, "-m", PYPROCT, pyProctControlFile]
    toRun = " ".join(toRun)
    print toRun
    startTime = time.time()
    proc = subprocess.Popen(toRun, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print out
    if err: print err

    endTime = time.time() 
    print "Pyproct took %.2f sec" % (endTime - startTime)

def getFileAndSnapshot(model):
    splitModel = model.split("\n")
    origFile = ""
    origSnapshot = 0
    for line in splitModel:
        if "REMARK source" in line:
            origFile = line.split(":")[1].strip()
        if "REMARK original model nr" in line:
            origSnapshot = int(line.split(":")[1]) - 1
            break
    
    return origFile, origSnapshot


def getSnapshot(file, snapshot):
    inputFile = open(file, "r")
    inputFileContent = inputFile.read()
    return inputFileContent.split("ENDMDL")[snapshot]

def calculateEachClustersPopulation(resultsJSON):
    jsonFile = open(resultsJSON, 'r').read()
    parsedJSON = json.loads(jsonFile)

    bestClustering = parsedJSON["best_clustering"]
    clusters = parsedJSON["selected"][bestClustering]["clustering"]
    clustersObject = Clustering.from_dic(clusters)
    sizes = np.array([])
    for cluster in clustersObject.clusters:
        sizes = np.append(sizes, cluster.get_size())
    return sizes

def calculateDegeneracy(sizes, spawningType, numberOfProcessors):
    if spawningType == "SAME_WEIGHT":
        return [1] * len(sizes)
    else:
        if isinstance(sizes, list): sizes = np.array(sizes)

        weights = 1./sizes
        weights /= sum(weights)

        degeneracy = []
        if spawningType == "INVERSELY_PROPORTIONAL":
            trajToDistribute = numberOfProcessors-1
            for i, weight in enumerate(weights):
                degeneracy.append(int(weight*trajToDistribute)) 

        elif spawningType == "AT_LEAST_ONE":
            numberOfInitialStructures = len(weights)
            trajToDistribute = numberOfProcessors-1-numberOfInitialStructures
            for i, weight in enumerate(weights):
                degeneracy.append(1 + int(weight*trajToDistribute)) #at least one per cluster


        decimalPart = []
        for i, weight in enumerate(weights):
            decimalPart.append(math.modf(weight*trajToDistribute)[0])
        sortedDecimals = np.argsort(decimalPart)
        sortedDecimals = sortedDecimals[::-1]

        leftProcessors = numberOfProcessors-1-sum(degeneracy)
        for i in range(leftProcessors):
            degeneracy[sortedDecimals[i]] += 1

        return degeneracy

def calculateDegeneracyOfClusterRepresentatives(representativesFile, resultsJSON, spawningType, numberOfProcessors):
    sizes = calculateEachClustersPopulation(resultsJSON)
    return calculateDegeneracy(sizes, spawningType, numberOfProcessors) 

def replaceLigandFileForFullProteinFile(origFile, ligandfileBasename, trajectoryBasename):
    newFile = origFile.replace(ligandfileBasename, trajectoryBasename)
    foundFiles = glob.glob(newFile)
    if len(foundFiles) > 1:
        sys.exit("Found " + str(foundFiles) + " that match pattern for full trajectory.")
    if len(foundFiles) == 0:
        sys.exit("Found no full trajectories that match pattern: \"" + trajectoryBasename + "\"")
    return foundFiles[0]

def makeClusterRepresentativesInitialStructures(representativesFile, tmpInitialStructuresTemplate, degeneracyOfRepresentatives, ligandfileBasename, trajectoryBasename):
    models = atomset.getPDBSnapshots(representativesFile)

    counts = 0
    for i, model in enumerate(models):
        origFile, origSnapshot = getFileAndSnapshot(model)
        replacedFile = replaceLigandFileForFullProteinFile(origFile, ligandfileBasename, trajectoryBasename) 
        origModel = getSnapshot(replacedFile, origSnapshot)

        for j in range(int(degeneracyOfRepresentatives[i])):
            outputFilename = tmpInitialStructuresTemplate%counts
            print 'Writing to ', outputFilename, replacedFile, origSnapshot
            representativeFile = open(outputFilename, 'w')
            representativeFile.write(origModel)
            representativeFile.close()

            counts += 1

    print "counts & cluster centers", counts, len(models)
    return counts

def makeOwnClusteringClusterRepresentativesInitialStructures(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clustering, iteration):

    counts = 0
    for i, cluster in enumerate(clustering.clusters.clusters):
        for j in range(int(degeneracyOfRepresentatives[i])):
            outputFilename = tmpInitialStructuresTemplate%(iteration,counts)
            print 'Writing to ', outputFilename, 'cluster', i
            cluster.writePDB(outputFilename)

            counts += 1

    print "counts & cluster centers", counts, len(np.where(np.array(degeneracyOfRepresentatives) > 0))
    return counts

def findFirstRun(outputPath, CLUSTERING_OUTPUT_OBJECT):
    """ Assumes that the outputPath is XXX/%d """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    for epoch in epochFolders:
        if os.path.exists(CLUSTERING_OUTPUT_OBJECT%epoch):
            return epoch
    if epoch <= 0:
        return 0

def spawningBuilder(spawningBlock):
    spawningTypeString = spawningBlock['type']

    if spawningTypeString == blockNames.STRING_SPAWNING_TYPES.sameWeight:
        spawningType = startingConformationsCalculator.SPAWNING_TYPES.sameWeight
    elif spawningTypeString == blockNames.STRING_SPAWNING_TYPES.inverselyProportional:
        spawningType = startingConformationsCalculator.SPAWNING_TYPES.inverselyProportional
    elif spawningTypeString == blockNames.STRING_SPAWNING_TYPES.epsilon:
        spawningType = startingConformationsCalculator.SPAWNING_TYPES.epsilon
    elif spawningTypeString == blockNames.STRING_SPAWNING_TYPES.FAST:
        spawningType = startingConformationsCalculator.SPAWNING_TYPES.FAST
    else:
        sys.exit("Unknown spawning type! Choices are: " + str(blockNames.SPAWNING_TYPE_TO_STRING_DICTIONARY.values()))

    spawningCalculator = None
    if spawningType == startingConformationsCalculator.SPAWNING_TYPES.sameWeight:
        spawningCalculator = SameWeightDegeneracyCalculator()
    elif spawningType == startingConformationsCalculator.SPAWNING_TYPES.inverselyProportional:
        spawningCalculator = startingConformationsCalculator.InverselyProportionalToPopulationCalculator()
    elif spawningType == startingConformationsCalculator.SPAWNING_TYPES.epsilon:
        spawningCalculator = startingConformationsCalculator.EpsilonDegeneracyCalculator()
    elif spawningType == startingConformationsCalculator.SPAWNING_TYPES.FAST:
        spawningCalculator = startingConformationsCalculator.FASTDegeneracyCalculator()


    spawningParams = startingConformationsCalculator.SpawningParams()
    if spawningType == startingConformationsCalculator.SPAWNING_TYPES.epsilon:
        spawningParamsBlock = spawningBlock['params']
        
        spawningParams.epsilon = spawningParamsBlock[blockNames.SPAWNING_PARAMS.EPSILON]
        spawningParams.reportFilename = spawningParamsBlock[blockNames.SPAWNING_PARAMS.REPORT_FILENAME]
        spawningParams.reportCol = spawningParamsBlock[blockNames.SPAWNING_PARAMS.REPORT_COL]
        spawningParams.temperature = spawningParamsBlock[blockNames.SPAWNING_PARAMS.TEMPERATURE]

    return spawningCalculator, spawningParams
    
def loadParams(jsonParams):
    """
        TODO: change for variables in a block names file, and work it out a bit more
    """
    jsonFile = open(jsonParams, 'r').read()
    parsedJSON = json.loads(jsonFile)

    return parsedJSON["RESTART"],parsedJSON['spawning'],\
            parsedJSON['outputPath'], parsedJSON['initialStructures'],\
            parsedJSON['ligandResname'].upper(), parsedJSON['debug'],\
            parsedJSON[blockNames.SIMULATION_BLOCK.blockname]

def saveInitialControlFile(jsonParams, originalControlFile):
    file = open(originalControlFile, 'w')
    jsonFile = open(jsonParams, 'r').read()
    file.write(jsonFile)



def grepTraj(epoch, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory, runNumber):
    fullTrajFilename = os.path.join(outputPathTempletized%epoch, trajectoryWildcard%runNumber)
    ligandTrajFilename = os.path.join(outputPathTempletized%epoch, ligandTrajectory%str(runNumber))
    toRun = "grep -E \"^HETATM...........%s|ATOM.........CA|MODEL|ENDMDL\" "%ligandResname + fullTrajFilename + " | grep -v \"\-\-\|REMARK\" > " + ligandTrajFilename
    #toRun = "grep -E \"^HETATM...........%s|MODEL|ENDMDL\" "%ligandResname + fullTrajFilename + " | grep -v \"\-\-\|REMARK\" > " + ligandTrajFilename
    #toRun = "cp " + fullTrajFilename + " " + ligandTrajFilename
    proc = subprocess.Popen(toRun, shell=True)
    proc.communicate()

def grepLigand(processors, epoch, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory):

    func = partial(grepTraj, epoch, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory)

    iterable = range(1,processors)

    pool = multiprocessing.Pool(processors)
    pool.map(func, iterable)
    pool.close()
    pool.join()

    #for runNumber in range(1, processors):
    #    grepTraj(runNumber, epoch, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory)

""" ORIG SERIAL
def grepLigand(processors, epoch, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory):
    for runNumber in range(1, processors):
        fullTrajFilename = os.path.join(outputPathTempletized%epoch, trajectoryWildcard%runNumber)
        ligandTrajFilename = os.path.join(outputPathTempletized%epoch, ligandTrajectory%str(runNumber))
        toRun = "grep -E \"^HETATM...........%s|ATOM.........CA|MODEL|ENDMDL\" "%ligandResname + fullTrajFilename + " | grep -v \"\-\-\|REMARK\" > " + ligandTrajFilename
        #toRun = "grep -E \"^HETATM...........%s|MODEL|ENDMDL\" "%ligandResname + fullTrajFilename + " | grep -v \"\-\-\|REMARK\" > " + ligandTrajFilename
        #toRun = "cp " + fullTrajFilename + " " + ligandTrajFilename
        proc = subprocess.Popen(toRun, shell=True)
        proc.communicate()
"""

"""
 #MAYBE GOOD TO HAVE FOR A NEAR FUTURE
 def read_in_chunks(size_in_bytes, processors, epoch, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory):
    s = 'Lets say i have a text file of 1000 GB'
    with open('data.txt', 'r+b') as f:
        prev = ''
        count = 0
        f_read  = partial(f.read, size_in_bytes)
        for text in iter(f_read, ''):
            if not text.endswith('\n'):
                # if file contains a partial line at the end, then don't
                # use it when counting the substring count. 
                text, rest = text.rsplit('\n', 1)
                # pre-pend the previous partial line if any.
                text =  prev + text
                prev = rest
            else:
                # if the text ends with a '\n' then simple pre-pend the
                # previous partial line. 
                text =  prev + text
                prev = ''
            count += text.count(s)
        count += prev.count(s)
        print count
"""

def writeClusteringOutput(outputPath, clustering, degeneracy, outputObject):
    cleanup(outputPath) 
    makeFolder(outputPath)

    for i in enumerate(clustering.clusters.clusters):
        i = i[0]
        outputFilename = "cluster_%d.pdb"%i
        outputFilename = os.path.join(outputPath, outputFilename)
        clustering.clusters.clusters[i].writePDB(outputFilename)

    sizes = []
    thresholds = []
    contacts = []
    energies = []
    for cluster in clustering.clusters.clusters:
        sizes.append(cluster.elements)
        thresholds.append(cluster.threshold)
        contacts.append(cluster.contacts)
        energies.append(cluster.getMetric())

    summaryFilename = os.path.join(outputPath, "summary.txt")
    summaryFile = open(summaryFilename, 'w')
    summaryFile.write("#cluster size degeneracy threshold contacts metric\n")
    for i in enumerate(sizes):
        i = i[0]
        summaryFile.write("%d %d %d %.1f %d %.1f\n"%(i, sizes[i], degeneracy[i], thresholds[i], contacts[i], energies[i]))
    summaryFile.close()

    with open(outputObject, 'wb') as f:
        pickle.dump(clustering, f, pickle.HIGHEST_PROTOCOL)



def main(jsonParams=None):
    if jsonParams is None:
        jsonParams = sys.argv[1]

    RESTART, spawningBlock, outputPath, initialStructures, ligandResname, DEBUG, simulationrunnerBlock = loadParams(jsonParams)

    startingConformationsCalculator, spawningParams = spawningBuilder(spawningBlock)
    runnerbuilder = SimulationRunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    print "================================"
    print "            PARAMS              "
    print "================================"
    print "Restarting simulations", RESTART
    print "Debug:", DEBUG

    print "Iterations: %d, Mpi processors: %d, Pele steps: %d"%(simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps)

    print "SpawningType:", blockNames.SPAWNING_TYPE_TO_STRING_DICTIONARY[startingConformationsCalculator.type]

    print "Output path: ", outputPath
    print "Initial Structures: ", initialStructures
    print "================================\n\n"


    #PRIVATE CONSTANTS
    ORIGINAL_CONTROLFILE = os.path.join(outputPath, "originalControlFile.conf")

    outputPathTempletized = os.path.join(outputPath, "%d")

    tmpFolder = "tmp_" +  outputPath.replace("/", "_")

    inputFileTemplate = "{ \"files\" : [ { \"path\" : \"%s\" } ] }"
    tmpInitialStructuresTemplate = tmpFolder+"/initial_%d_%d.pdb"
    tmpControlFilename = tmpFolder+"/controlFile%d.conf"
    tmpControlFilename_pyproct = tmpFolder+"/pyproct_controlFile%d.conf"
    tmpControlFilename_pyproct_thisEpoch = tmpFolder+"/pyproct_thisEpoch_controlFile%d.conf"
    trajectoryBasename = "*traj*"
    ligandTrajectoryBasename = "ligandSnapshots"

    trajectoryWildcard = trajectoryBasename + "_%d.pdb"
    ligandTrajectory = ligandTrajectoryBasename + "_%s.pdb" #e.g. ligandSnapshots_%s.pdb



    PYPROCT = "pyproct.main"
    PYPROCT_OUTPUT_DIR = outputPathTempletized+"/pyproct"
    PYPROCT_OUTPUT_DIR_SINGLE_EPOCH = outputPathTempletized+"/pyproct_thisEpoch"
    PYPROCT_REPRESENTATIVE_OUTPUT = PYPROCT_OUTPUT_DIR + "/results/representatives.pdb"
    PYPROCT_RESULTS_OUTPUT = PYPROCT_OUTPUT_DIR + "/results/results.json"
    PYPROCT_COMPRESSION_OUTPUT = PYPROCT_OUTPUT_DIR_SINGLE_EPOCH + "/results/compressed.pdb"

    PYPROCT_COMPRESSION_BLOCK = ",\"compression\":{\n\
            \"keep_remarks\": \"ALL\",\n\
            \"add_source_details\": true,\n\
            \"final_number_of_frames\": %d,\n\
            \"file\": \"compressed.pdb\",\n\
            \"type\":\"KMEDOIDS\"\n\
        }"#note leading comma

    CLUSTERING_OUTPUT_DIR = os.path.join(outputPathTempletized, "clustering")
    CLUSTERING_OUTPUT_OBJECT = os.path.join(CLUSTERING_OUTPUT_DIR, "object.pkl")
    #END PRIVATE CONSTANTS


    #if not DEBUG: atexit.register(cleanup, tmpFolder)

    makeFolder(outputPath)
    makeFolder(tmpFolder)
    saveInitialControlFile(jsonParams, ORIGINAL_CONTROLFILE)



    if RESTART:
        firstRun = findFirstRun(outputPath, CLUSTERING_OUTPUT_OBJECT)

        if firstRun != 0:
            """
            degeneracyOfRepresentatives = calculateDegeneracyOfClusterRepresentatives(PYPROCT_REPRESENTATIVE_OUTPUT%(firstRun-1), PYPROCT_RESULTS_OUTPUT%(firstRun-1), spawning, processors)
            seedingPoints = makeClusterRepresentativesInitialStructures(PYPROCT_REPRESENTATIVE_OUTPUT%(firstRun-1), tmpInitialStructuresTemplate, degeneracyOfRepresentatives, ligandTrajectoryBasename, trajectoryBasename)
            """



            #TO BE REMOVED
            """
            snapshotsJSONSelectionString = generateSnapshotSelectionString(firstRun-1, outputPathTempletized, trajectoryBasename)
            snapshotsJSONSelectionString = "[" + snapshotsJSONSelectionString + "]"
            paths = eval(snapshotsJSONSelectionString)
            clusteringMethod = clustering.Clustering(ligandResname)
            clusteringMethod.cluster(paths, withinClusterThreshold)
            """

            with open(CLUSTERING_OUTPUT_OBJECT%(firstRun), 'rb') as f:
                clusteringMethod = pickle.load(f)


            sizes = []
            for cluster in clusteringMethod.clusters.clusters:
                sizes.append(cluster.elements)
            #degeneracyOfRepresentatives = calculateDegeneracy(sizes, spawning, processors)
            degeneracyOfRepresentatives = startingConformationsCalculator.calculate(clusteringMethod.clusters.clusters, processors-1, spawningParams, firstRun)
            startingConformationsCalculator.log()
            print "Degeneracy", degeneracyOfRepresentatives

            seedingPoints = makeOwnClusteringClusterRepresentativesInitialStructures(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clusteringMethod, firstRun)

            initialStructuresAsString = createMultipleComplexesFilenames(seedingPoints, inputFileTemplate, tmpInitialStructuresTemplate, firstRun)
            
    if not RESTART or firstRun == 0: # if RESTART and firstRun = 0, it must check into the initial structures
        #Choose initial structures
        if not DEBUG: shutil.rmtree(outputPath)
        makeFolder(outputPath)
        firstRun = 0
        saveInitialControlFile(jsonParams, ORIGINAL_CONTROLFILE)
        copyInitialStructures(initialStructures, tmpInitialStructuresTemplate, firstRun)
        initialStructuresAsString = createMultipleComplexesFilenames(len(initialStructures), inputFileTemplate, tmpInitialStructuresTemplate, firstRun)


    peleControlFileDictionary = {"COMPLEXES":initialStructuresAsString, "PELE_STEPS":simulationRunner.parameters.peleSteps}
    """
    pyproctControlFileDictionary = {"MIN_SEEDING_POINTS":minSeedingPoints, "MAX_SEEDING_POINTS":maxSeedingPoints, "NUMBER_OF_PROCESSORS":processors, "RESNAME": ligandResname}

    if DEBUG:
        pyproctControlFileDictionary["NUMBER_OF_PROCESSORS"] = 8
    """

    #Make control file
    outputDir = outputPathTempletized%firstRun
    makeFolder(outputDir)
    peleControlFileDictionary["OUTPUT_PATH"] = outputDir
    peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + firstRun*simulationRunner.parameters.processors
    simulationRunner.makeWorkingControlFile(simulationRunner.parameters.runningControlfilename, tmpControlFilename%firstRun, peleControlFileDictionary) 

    for i in range(firstRun, simulationRunner.parameters.iterations):
        print "Iteration", i
        
        print "Production run..."
        if not DEBUG:
            startTime = time.time() 
           # runPeleSimulations(processors, tmpControlFilename%i, PELE_EXECUTABLE, DATA_FOLDER, DOCUMENTS_FOLDER)
            simulationRunner.run_simulation()
            endTime = time.time() 
            print "PELE %s sec" % (endTime - startTime)

            #filter ligand
            """
            startTime = time.time() 
            grepLigand(processors, i, ligandResname, outputPathTempletized, trajectoryWildcard, ligandTrajectory)
            endTime = time.time() 
            print "Grep ligand: %s sec" % (endTime - startTime)
            """

        #allsnapshotsjsonselectionstring = generatesnapshotselectionstring(i)
        """
        lastRoundSnapshotsJSONSelectionString = generateSnapshotSelectionStringLastRound(i, outputPathTempletized, ligandTrajectory%"*")
        pyproctControlFileDictionary["COMPLEXES"] = lastRoundSnapshotsJSONSelectionString
        pyproctControlFileDictionary["OUTPUT_PATH"] = PYPROCT_OUTPUT_DIR_SINGLE_EPOCH%i
        typicalAcceptance = 0.3
        compressionFactor = 0.1
        estimatedAcceptedSteps = (processors-1)*peleSteps*typicalAcceptance*compressionFactor
        estimatedAcceptedSteps = max(1,estimatedAcceptedSteps) #at least one
        pyproctControlFileDictionary["COMPRESSION"] = PYPROCT_COMPRESSION_BLOCK%min(estimatedAcceptedSteps, countEpochAcceptedSteps(outputPathTempletized%i))
        pyproctControlFileDictionary["CLUSTERING"] = "gromos"
        #pyproctControlFileDictionary["MAX_SEEDING_POINTS"] = 32
        makeWorkingControlFile(TEMPLETIZED_PYPROCT_CONTROLFILENAME, tmpControlFilename_pyproct_thisEpoch%i, pyproctControlFileDictionary) 


        print "Filtering current epoch..."
        startTime = time.time() 
        runPyProCT(tmpControlFilename_pyproct_thisEpoch%i, PYTHON, PYPROCT) #Handle case when pyProCt fails
        endTime = time.time() 
        print "Pyproct filtering ligand: %s sec" % (endTime - startTime)


        print "Clustering with filtered snapshots"
        filteredList = generateFilteredListSnapshotSelectionString(i, PYPROCT_COMPRESSION_OUTPUT)
        pyproctControlFileDictionary["COMPLEXES"] = filteredList
        pyproctControlFileDictionary["OUTPUT_PATH"] = PYPROCT_OUTPUT_DIR%i
        pyproctControlFileDictionary["COMPRESSION"] = "" #won't be use, we'll use the single trajectory compressed/filtered list
        pyproctControlFileDictionary["CLUSTERING"] = "gromos"
        #pyproctControlFileDictionary["MAX_SEEDING_POINTS"] = maxSeedingPoints
        makeWorkingControlFile(TEMPLETIZED_PYPROCT_CONTROLFILENAME, tmpControlFilename_pyproct%i, pyproctControlFileDictionary) 

        print "Clustering..."
        startTime = time.time() 
        runPyProCT(tmpControlFilename_pyproct%i, PYTHON, PYPROCT) #Handle case when pyProCt fails
        endTime = time.time() 
        print "Pyproct clustering ligand: %s sec" % (endTime - startTime)
        """

        print "Clustering..."
        startTime = time.time() 
        #snapshotsJSONSelectionString = generateSnapshotSelectionString(i, outputPathTempletized, trajectoryBasename)
        snapshotsJSONSelectionString = generateSnapshotSelectionStringLastRound(i, outputPathTempletized, trajectoryBasename)
        snapshotsJSONSelectionString = "[" + snapshotsJSONSelectionString + "]"
        paths = eval(snapshotsJSONSelectionString)
        if len(glob.glob(paths[-1])) == 0: sys.exit("No more trajectories to cluster")
        if i == 0:
            clusteringMethod = clustering.Clustering(ligandResname, spawningParams.reportFilename, spawningParams.reportCol)
        #else:
        #    #CAN'T THIS BE REMOVED????
        #    with open(CLUSTERING_OUTPUT_OBJECT%(i-1), 'rb') as input:
        #        clusteringMethod = pickle.load(input)
        clusteringMethod.cluster(paths)
        endTime = time.time() 
        print "Clustering ligand: %s sec" % (endTime - startTime)


        """
        degeneracyOfRepresentatives = calculateDegeneracyOfClusterRepresentatives(PYPROCT_REPRESENTATIVE_OUTPUT%i, PYPROCT_RESULTS_OUTPUT%i, spawning, processors)
        numberOfSeedingPoints = makeClusterRepresentativesInitialStructures(PYPROCT_REPRESENTATIVE_OUTPUT%i, tmpInitialStructuresTemplate, degeneracyOfRepresentatives, ligandTrajectoryBasename, trajectoryBasename)
        """
        degeneracyOfRepresentatives = startingConformationsCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, i)
        startingConformationsCalculator.log()
        print "Degeneracy", degeneracyOfRepresentatives

        numberOfSeedingPoints = makeOwnClusteringClusterRepresentativesInitialStructures(tmpInitialStructuresTemplate, degeneracyOfRepresentatives, clusteringMethod, i+1)

        writeClusteringOutput(CLUSTERING_OUTPUT_DIR%i, clusteringMethod, degeneracyOfRepresentatives, CLUSTERING_OUTPUT_OBJECT%i)

        #Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            initialStructuresAsString = createMultipleComplexesFilenames(numberOfSeedingPoints, inputFileTemplate, tmpInitialStructuresTemplate, i+1)
            peleControlFileDictionary["COMPLEXES"] = initialStructuresAsString

            outputDir = outputPathTempletized%(i+1)
            makeFolder(outputDir) #PELE does not do it automatically
            peleControlFileDictionary["OUTPUT_PATH"] = outputDir
            peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + (i+1)*simulationRunner.parameters.processors
            simulationRunner.makeWorkingControlFile(simulationRunner.parameters.runningControlfilename, tmpControlFilename%(i+1), peleControlFileDictionary) 

    #cleanup
    #cleanup(tmpFolder)

if __name__ == '__main__':
    main()
