import time
import constants
import subprocess
import blockNames

class SIMULATION_TYPE:
    PELE, MD, TEST = range(3)

class SimulationParameters:
    def __init__(self):
        self.processors = 0
        self.executable = ""
        self.runningControlfilename = ""
        self.Datafolder = ""
        self.Documentsfolder = ""

class SimulationRunner:
    def __init__(self, parameters):
        self.parameters = parameters
    
    def run_simulation(self):
        pass 

class PeleSimulation(SimulationRunner):
    
    def createSymbolicLinks(self):
        if not os.path.islink("Data"):
            os.system("ln -s " + self.parameters.DATA_FOLDER + " Data")
        if not os.path.islink("Documents"):
            os.system("ln -s " + self.parameters.DOCUMENTS_FOLDER + " Documents")

    def run_simulation(self):
        createSymbolicLinks()

        toRun = ["mpirun -np " + str(self.parameters.processors), self.parameters.executable, self.parameters.runningControlfilename]
        toRun = " ".join(toRun)
        print toRun
        startTime = time.time() 
        proc = subprocess.Popen(toRun, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        print out
        if err: print err

        endTime = time.time() 
        print "PELE took %.2f sec" % (endTime - startTime)

class TestSimulation(SimulationRunner):
    
    def run_simulation(self):
        pass

class RunnerBuilder:

    def build(self, simulationRunnerBlock):
        simulationType = simulationRunnerBlock[blockNames.SIMULATION_TYPE.type]
        paramsBlock = simulationRunnerBlock[blockNames.SIMULATION_PARAMS.params]
        params = SimulationParameters()
        if simulationType == blockNames.SIMULATION_TYPE.PELE:
            params.processors = paramsBlock[blockNames.SIMULATION_PARAMS.processors] 
            params.Datafolder = paramsBlock.get(blockNames.SIMULATION_PARAMS.Datafolder, constants.DATA_FOLDER) 
            params.Documentsfolder = paramsBlock.get(blockNames.SIMULATION_PARAMS.Documentsfolder, constants.DOCUMENTS_FOLDER) 
            params.executable = paramsBlock.get(blockNames.SIMULATION_PARAMS.executable, constants.PELE_EXECUTABLE) 
            params.runningControlfilename = paramsBlock[blockNames.SIMULATION_PARAMS.runningControlfilename]
            SimulationRunner = PeleSimulation(params)
        elif simulationType == blockNames.SIMULATION_TYPE.MD:
            pass
        elif simulationType == blockNames.SIMULATION_TYPE.TEST:
            SimulationRunner = TestSimulation(params)
        else:
            sys.exit("Unknown simulation type! Choices are: " + str(blockNames.SIMULATION_TYPE_TO_STRING_DICTIONARY.values()))
        return SimulationRunner

