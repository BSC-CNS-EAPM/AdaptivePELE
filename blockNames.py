class StringSpawningTypes:
    type = 'type'
    sameWeight = "sameWeight"
    inverselyProportional = "inverselyProportional"
    epsilon = "epsilon"
    fast = "FAST"
    
class SpawningParams:
    params = "params"
    epsilon = "epsilon"
    temperature = "T"
    threshold = "threshold"
    report_filename = "reportFilename"
    report_col = "metricColumnInReport"

class SimulationType:
    type = "type"
    pele = "pele"
    md = "md"
    test  = "test"

class SimulationParams:
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

class ControlFileParams:
    restart = 'restart'
    spawningBlockname = 'spawning'
    outputPath = 'outputPath'
    initialStructures = 'initialStructures'
    ligandResname = 'ligandResname'
    debug = 'debug'
    simulationBlockname = "simulation"
