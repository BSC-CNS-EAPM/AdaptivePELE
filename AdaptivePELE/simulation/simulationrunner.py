import time
import os
import json
import subprocess
import shutil
import string
import sys
from AdaptivePELE.constants import constants, blockNames
from AdaptivePELE.simulation import simulationTypes
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities


class SimulationParameters:
    def __init__(self):
        self.processors = 0
        self.executable = ""
        self.templetizedControlFile = ""
        self.dataFolder = ""
        self.documentsFolder = ""
        self.iterations = 0
        self.peleSteps = 0
        self.seed = 0
        self.exitCondition = None
        self.boxCenter = None
        self.boxRadius = 20
        self.modeMovingBox = None
        self.runEquilibration = False
        self.destination = None
        self.origin = None


class SimulationRunner:
    def __init__(self, parameters):
        self.parameters = parameters
        self.processorsToClusterMapping = []

    def runSimulation(self, runningControlFile=""):
        pass

    def hasExitCondition(self):
        """
            Check if an exit condition has been set

            :returns: bool -- True if an exit condition is set
        """
        return self.parameters.exitCondition is not None

    def checkExitCondition(self, clustering):
        """
            Check if the exit condition has been met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- True if the exit condition is met
        """
        if self.parameters.exitCondition:
            return self.parameters.exitCondition.checkExitCondition(clustering)
        return False

    def makeWorkingControlFile(self, workingControlFilename, dictionary, inputTemplate=None):
        """
            Substitute the values in the templetized control file

            :param workingControlFilename: Name of the template control file
            :type workingControlFilename: str
            :param dictionary: Dictonary containing the parameters to substitute
                in the control file
            :param inputFileTemplate: Template control file
            :type inputFileTemplate: str
            :type dictionary: dict
        """
        if inputTemplate is None:
            with open(self.parameters.templetizedControlFile, "r") as inputFile:
                inputFileContent = inputFile.read()
        else:
            inputFileContent = inputTemplate

        inputFileTemplate = string.Template(inputFileContent)
        outputFileContent = inputFileTemplate.substitute(dictionary)

        with open(workingControlFilename, "w") as outputFile:
            outputFileContent = outputFileContent.replace("'", '"')
            outputFile.write(outputFileContent)

    def updateMappingProcessors(self, mapping):
        """
            Update the value of the processorsToClusterMapping, a list with the
            snapshot from which the trajectories will start in the next iteration

            :param mapping: List with the snapshot from which the trajectories
                will start in the next iteration
            :type mapping: list
        """
        self.processorsToClusterMapping = mapping[1:]+[mapping[0]]

    def writeMappingToDisk(self, epochDir):
        """
            Write the processorsToClusterMapping to disk

            :param epochDir: Name of the folder where to write the
                processorsToClusterMapping
            :type epochDir: str
        """
        if len(self.processorsToClusterMapping) == 0:
            return
        with open(epochDir+"/processorMapping.txt", "w") as f:
            f.write(':'.join(map(str, self.processorsToClusterMapping)))

    def readMappingFromDisk(self, epochDir):
        """
            Read the processorsToClusterMapping from disk

            :param epochDir: Name of the folder where to write the
                processorsToClusterMapping
            :type epochDir: str
        """
        try:
            with open(epochDir+"/processorMapping.txt") as f:
                self.processorsToClusterMapping = map(int, f.read().rstrip().split(':'))
        except IOError:
            sys.stderr.write("WARNING: processorMapping.txt not found, you might not be able to recronstruct fine-grained pathways\n")

    def setZeroMapping(self):
        """
            Set the processorsToClusterMapping to zero
        """
        self.processorsToClusterMapping = [0 for _ in xrange(1, self.parameters.processors)]


class PeleSimulation(SimulationRunner):
    def __init__(self, parameters):
        SimulationRunner.__init__(self, parameters)
        self.type = simulationTypes.SIMULATION_TYPE.PELE

    def createSymbolicLinks(self):
        """
            Create symbolic links to Data and Documents folder if they don't exist
        """
        if not os.path.islink("Data"):
            os.system("ln -s " + self.parameters.dataFolder + " Data")
        if not os.path.islink("Documents"):
            os.system("ln -s " + self.parameters.documentsFolder + " Documents")

    def selectInitialBoxCenter(self, initialStructuresAsString, resname):
        """
            Select the coordinates of the first box, currently as the center of
            mass of the first initial structure provided

            :param initialStructuresAsString: String containing the files of the initial structures
            :type initialStructuresAsString: str
            :param resname: Residue name of the ligand in the system pdb
            :type resname: str

            :returns str: -- string to be substitued in PELE control file
        """
        # This is a dictionary because it's prepared to be subtitued in the PELE
        # control files (i.e. JSON format)
        try:
            initialStructuresDict = json.loads(initialStructuresAsString.split(",")[0])
            initialStruct = str(initialStructuresDict['files'][0]['path'])
        except ValueError:
            # If a json valid string is not passed, interpret the input string
            # as a path to a file
            initialStruct = initialStructuresAsString

        PDBinitial = atomset.PDB()
        PDBinitial.initialise(initialStruct, resname=resname)
        return repr(PDBinitial.getCOM())

    def runSimulation(self, runningControlFile=""):
        """
            Run a short PELE simulation

            :param runningControlFile: PELE control file to run
            :type runningControlFile: str
        """
        self.createSymbolicLinks()

        toRun = ["mpirun -np " + str(self.parameters.processors), self.parameters.executable, runningControlFile]
        toRun = " ".join(toRun)
        print toRun
        startTime = time.time()
        proc = subprocess.Popen(toRun, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        print out
        if err:
            print err

        endTime = time.time()
        print "PELE took %.2f sec" % (endTime - startTime)

    def getEquilibrationControlFile(self, peleControlFileDict):
        """
            Filter unnecessary parameters and return a minimal PELE control
            file for equilibration runs

            :param peleControlFileDict: Dictionary with pele control file options
            :type peleControlFileDict: dict

            :returns: dict -- Dictionary with pele control file options
        """
        # Set small rotations and translations
        peleControlFileDict["commands"][0]["Perturbation"]["translationRange"] = 0.5
        peleControlFileDict["commands"][0]["Perturbation"]["rotationScalingFactor"] = 0.01
        # Remove dynamical changes in control file
        peleControlFileDict["commands"][0]["PeleTasks"][0].pop("exitConditions", None)
        peleControlFileDict["commands"][0]["PeleTasks"][0].pop("parametersChanges", None)

        return peleControlFileDict

    def equilibrate(self, intialStructures, outputPathConstants, reportFilename, outputPath, resname):
        """
            Run short simulation to equilibrate the system. It will run one
            such simulation for every initial structure and select appropiate
            structures to start the simulation

            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: list of str
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param reportBaseFilename: Name of the file that contains the metrics of the snapshots to cluster
            :type reportBaseFilename: str
            :param outputPath: Path where trajectories are found
            :type outputPath: str
            :param resname: Residue name of the ligand in the system pdb
            :type resname: str

            :returns: list --  List with initial structures
        """
        newInitialStructures = []
        equilibrationPeleDict = {"PELE_STEPS": 50, "SEED": self.parameters.seed}
        with open(self.parameters.templetizedControlFile) as fc:
            peleControlFile = fc.read()

        templateNames = {ele[1]: '"$%s"' % ele[1] for ele in string.Template.pattern.findall(peleControlFile)}
        templateNames.pop("OUTPUT_PATH", None)
        peleControlFileDict = json.loads(string.Template(peleControlFile).safe_substitute(templateNames))
        peleControlFileDict = self.getEquilibrationControlFile(peleControlFileDict)

        for i, structure in enumerate(intialStructures):
            equilibrationOutput = os.path.join(outputPath, "equilibration_%d" % (i+1))
            equilibrationControlFile = outputPathConstants.tmpControlFilename % (i+1)
            utilities.makeFolder(equilibrationOutput)
            shutil.copyfile(structure, outputPathConstants.tmpInitialStructuresEquilibrationTemplate % (i+1))
            initialStructureString = self.createMultipleComplexesFilenames(1, outputPathConstants.tmpInitialStructuresEquilibrationTemplate, i+1, equilibration=True)
            equilibrationPeleDict["OUTPUT_PATH"] = equilibrationOutput
            equilibrationPeleDict["COMPLEXES"] = initialStructureString
            equilibrationPeleDict["BOX_CENTER"] = self.selectInitialBoxCenter(structure, resname)
            equilibrationPeleDict["BOX_RADIUS"] = self.parameters.boxRadius
            print "Running equilibration for initial structure number %d" % (i+1)
            peleControlString = json.dumps(peleControlFileDict, indent=4)
            for key, value in templateNames.iteritems():
                # Remove double quote around template keys, so that PELE
                # understands the options
                peleControlString.replace(value, '$%s' % key)
            self.makeWorkingControlFile(equilibrationControlFile, equilibrationPeleDict, )
            self.runSimulation(equilibrationControlFile)
            # TODO: Implement a selectEquilibratedStructure procedure
            # newInitialStructures.append(self.selectEquilibratedStructure())
            newInitialStructures.append(structure)
        return newInitialStructures

    def createMultipleComplexesFilenames(self, numberOfSnapshots, tmpInitialStructuresTemplate, iteration, equilibration=False):
        """
            Creates the string to substitute the complexes in the PELE control file

            :param numberOfSnapshots: Number of complexes to write
            :type numberOfSnapshots: int
            :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
            :type tmpInitialStructuresTemplate: str
            :param iteration: Epoch number
            :type iteration: int
            :param equilibration: Flag to mark wether the complexes are part of an
            equilibration run
            :type equilibration: bool

            :returns: str -- jsonString to be substituted in PELE control file
        """
        jsonString = ["\n"]
        for i in range(numberOfSnapshots-1):
            if equilibration:
                jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration)) + ",\n")
            else:
                jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, i)) + ",\n")
        if equilibration:
            jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration)))
        else:
            jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, numberOfSnapshots-1)))
        return "".join(jsonString)


class TestSimulation(SimulationRunner):
    """
        Class used for testing
    """
    def __init__(self, parameters):
        SimulationRunner.__init__(self, parameters)
        self.type = simulationTypes.SIMULATION_TYPE.TEST
        self.copied = False
        self.parameters = parameters

    def runSimulation(self, runningControlFile=""):
        """
            Copy file to test the rest of the AdaptivePELE procedure
        """
        if not self.copied:
            if os.path.exists(self.parameters.destination):
                shutil.rmtree(self.parameters.destination)
            shutil.copytree(self.parameters.origin, self.parameters.destination)
            self.copied = True

    def makeWorkingControlFile(self, workingControlFilename, dictionary):
        pass


class ExitConditionBuilder:
    def build(self, exitConditionBlock):
        """
            Build the selected exit condition object

            :param exitConditionBlock: Block of the control file
                corresponding to the exit condition
            :type exitConditionBlock: dict

            :returns: :py:class:`.MetricExitCondition` -- MetricExitCondition object
                selected
        """
        exitConditionType = exitConditionBlock[blockNames.ExitConditionType.type]
        exitConditionParams = exitConditionBlock[blockNames.SimulationParams.params]
        if exitConditionType == blockNames.ExitConditionType.metric:
            metricCol = exitConditionParams[blockNames.SimulationParams.metricCol]
            metricValue = exitConditionParams[blockNames.SimulationParams.exitValue]
            return MetricExitCondition(metricCol, metricValue)
        elif exitConditionType == blockNames.ExitConditionType.clustering:
            ntrajs = exitConditionParams[blockNames.SimulationParams.trajectories]
            return ClusteringExitCondition(ntrajs)
        else:
            sys.exit("Unknown exit condition type! Choices are: " + str(simulationTypes.EXITCONDITION_TYPE_TO_STRING_DICTIONARY.values()))


class ClusteringExitCondition:
    def __init__(self, ntrajs):
        self.clusterNum = 0
        self.ntrajs = ntrajs
        self.type = simulationTypes.EXITCONDITION_TYPE.CLUSTERING

    def checkExitCondition(self, clustering):
        """
            Iterate over all unchecked cluster and check if the exit condtion
            is met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- Returns True if the exit condition has been met
        """
        newClusterNum = clustering.getNumberClusters()
        clusterDiff = newClusterNum - self.clusterNum
        self.clusterNum = newClusterNum
        return clusterDiff < 0.1*self.ntrajs


class MetricExitCondition:
    def __init__(self, metricCol, metricValue):
        self.metricCol = metricCol
        self.metricValue = metricValue
        self.lastCheckedCluster = 0
        self.type = simulationTypes.EXITCONDITION_TYPE.METRIC

    def checkExitCondition(self, clustering):
        """
            Iterate over all unchecked cluster and check if the exit condtion
            is met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- Returns True if the exit condition has been met
        """
        for cluster in clustering.clusters.clusters:
            metric = cluster.getMetricFromColumn(self.metricCol)
            if metric is not None and metric < self.metricValue:
                return True
        return False


class RunnerBuilder:
    def build(self, simulationRunnerBlock):
        """
            Build the selected  SimulationRunner object

            :param simulationRunnerBlock: Block of the control file
                corresponding to the simulation step
            :type simulationRunnerBlock: dict

            :returns: :py:class:`.SimulationRunner` -- SimulationRunner object
                selected
        """
        simulationType = simulationRunnerBlock[blockNames.SimulationType.type]
        paramsBlock = simulationRunnerBlock[blockNames.SimulationParams.params]
        params = SimulationParameters()
        if simulationType == blockNames.SimulationType.pele:
            params.processors = paramsBlock[blockNames.SimulationParams.processors]
            params.dataFolder = paramsBlock.get(blockNames.SimulationParams.dataFolder, constants.DATA_FOLDER)
            params.documentsFolder = paramsBlock.get(blockNames.SimulationParams.documentsFolder, constants.DOCUMENTS_FOLDER)
            params.executable = paramsBlock.get(blockNames.SimulationParams.executable, constants.PELE_EXECUTABLE)
            params.templetizedControlFile = paramsBlock[blockNames.SimulationParams.templetizedControlFile]
            params.iterations = paramsBlock[blockNames.SimulationParams.iterations]
            params.peleSteps = paramsBlock[blockNames.SimulationParams.peleSteps]
            params.seed = paramsBlock[blockNames.SimulationParams.seed]
            params.modeMovingBox = paramsBlock.get(blockNames.SimulationParams.modeMovingBox)
            params.boxCenter = paramsBlock.get(blockNames.SimulationParams.boxCenter)
            params.boxRadius = paramsBlock.get(blockNames.SimulationParams.boxRadius, 20)
            params.runEquilibration = paramsBlock.get(blockNames.SimulationParams.runEquilibration, False)
            exitConditionBlock = paramsBlock.get(blockNames.SimulationParams.exitCondition, None)
            if exitConditionBlock:
                exitConditionBuilder = ExitConditionBuilder()
                params.exitCondition = exitConditionBuilder.build(exitConditionBlock)

            return PeleSimulation(params)
        elif simulationType == blockNames.SimulationType.md:
            pass
        elif simulationType == blockNames.SimulationType.test:
            params.processors = paramsBlock[blockNames.SimulationParams.processors]
            params.destination = paramsBlock[blockNames.SimulationParams.destination]
            params.origin = paramsBlock[blockNames.SimulationParams.origin]
            params.iterations = paramsBlock[blockNames.SimulationParams.iterations]
            params.peleSteps = paramsBlock[blockNames.SimulationParams.peleSteps]
            params.seed = paramsBlock[blockNames.SimulationParams.seed]
            return TestSimulation(params)
        else:
            sys.exit("Unknown simulation type! Choices are: " + str(simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY.values()))
