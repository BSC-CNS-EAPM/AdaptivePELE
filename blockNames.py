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

""" MOVE TO CONSTANTS """
import startingConformationsCalculator

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    startingConformationsCalculator.SPAWNING_TYPES.sameWeight:STRING_SPAWNING_TYPES.sameWeight, 
    startingConformationsCalculator.SPAWNING_TYPES.inverselyProportional:STRING_SPAWNING_TYPES.inverselyProportional, 
    startingConformationsCalculator.SPAWNING_TYPES.epsilon:STRING_SPAWNING_TYPES.epsilon,
    startingConformationsCalculator.SPAWNING_TYPES.FAST:STRING_SPAWNING_TYPES.FAST
}

""" END MOVE TO CONSTANTS """
