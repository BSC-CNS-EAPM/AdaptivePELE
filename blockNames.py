class STRING_SPAWNING_TYPES:
    sameWeight = "sameWeight"
    inverselyProportional = "inverselyProportional"
    epsilon = "epsilon"
    FAST = "FAST"
    
class SPAWNING_PARAMS:
    EPSILON = "epsilon"
    TEMPERATURE = "T"
    THRESHOLD = "threshold"
    REPORT_FILENAME = "reportFilename"
    REPORT_COL = "metricColumnInReport"

class SIMULATION_BLOCK:
    blockname = "simulation"

class SIMULATION_TYPE:
    type = "type"
    PELE = "pele"
    MD = "md"
    TEST  = "testing"

class SIMULATION_PARAMS:
    params = "params"
    processors = "processors"
    executable = "executable"
    runningControlfilename = "controlfile"
    Datafolder = "data"
    Documentsfolder = "documents"

""" MOVE TO CONSTANTS """
import startingConformationsCalculator
import SimulationRunner

SIMULATION_TYPE_TO_STRING_DICTIONARY = {
    SimulationRunner.SIMULATION_TYPE.PELE:SIMULATION_TYPE.PELE
    SimulationRunner.SIMULATION_TYPE.MD:SIMULATION_TYPE.MD
    SimulationRunner.SIMULATION_TYPE.TEST:SIMULATION_TYPE.TEST
    }

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    startingConformationsCalculator.SPAWNING_TYPES.sameWeight:STRING_SPAWNING_TYPES.sameWeight, 
    startingConformationsCalculator.SPAWNING_TYPES.inverselyProportional:STRING_SPAWNING_TYPES.inverselyProportional, 
    startingConformationsCalculator.SPAWNING_TYPES.epsilon:STRING_SPAWNING_TYPES.epsilon,
    startingConformationsCalculator.SPAWNING_TYPES.FAST:STRING_SPAWNING_TYPES.FAST
}

""" END MOVE TO CONSTANTS """
