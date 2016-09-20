class ClusteringTypes:
    type = "type"
    contacts = "contacts"
    contactMap = "contactMap"
    agglomerative = "agglomerative"
    thresholdCalculator = "thresholdCalculator"
    ligandResname = "ligandResname"
    nclusters = "nclusters"

class ThresholdCalculator:
    type = "type"
    heaviside = "heaviside"
    constant = "constant"
    conditions = "conditions"
    values = "values"
    value = "value"

class DensityCalculator:
    type = "type"
    heaviside = "heaviside"
    null = "null"
    constant = "constant"
    conditions = "conditions"
    values = "values"

class StringSpawningTypes:
    type = "type"
    sameWeight = "sameWeight"
    inverselyProportional = "inverselyProportional"
    epsilon = "epsilon"
    fast = "FAST"
    ### New parameters for variable epsilon(experimental)
    variableEpsilon = "variableEpsilon"
    
class SpawningParams:
    params = "params"
    epsilon = "epsilon"
    temperature = "T"
    threshold = "threshold"
    report_filename = "reportFilename"
    report_col = "metricColumnInReport"
    ### New parameters for variable epsilon(experimental)
    varEpsilonType = "varEpsilonType"
    maxEpsilon = "maxEpsilon"
    minEpsilon = "minEpsilon"
    variationWindow = "variationWindow" # Last epoch of variable epsilon,if
    # current epoch > than variation Window, set epsilon to minEpsilon
    maxEpsilonWindow = "maxEpsilonWindow"
    period = "period" # Only useful for periodic epsilon modes
    density = "density"

class SpawningDensity:
    values = "values"
    conditions = "conditions"

class VariableEpsilonTypes:
    linearVariation = "linearVariation"

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
    generalParams = "generalParams"
    spawningBlockname = "spawning"
    simulationBlockname = "simulation"
    clusteringBlockname = "clustering"

class GeneralParams:
    restart = "restart"
    outputPath = "outputPath"
    initialStructures = "initialStructures"
    debug = "debug"
