class STRING_SPAWNING_TYPES:
    type = 'type'
    sameWeight = "sameWeight"
    inverselyProportional = "inverselyProportional"
    epsilon = "epsilon"
    FAST = "FAST"
    
class SPAWNING_PARAMS:
    params = "params"
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
    TEST  = "test"

class SIMULATION_PARAMS:
    params = "params"
    processors = "processors"
    executable = "executable"
    templetizedControlFile = "controlFile"
    dataFolder = "data"
    documentsFolder = "documents"
    destination = "destination"
    origin = "origin"
    seed = "seed"
    peleSteps = "peleSteps"
    iterations = "iterations"
