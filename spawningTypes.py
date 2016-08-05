import blockNames

class SPAWNING_TYPES:
    sameWeight, inverselyProportional, epsilon, simulatedAnnealing, FAST = range(5)

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    SPAWNING_TYPES.sameWeight:blockNames.StringSpawningTypes.sameWeight, 
    SPAWNING_TYPES.inverselyProportional:blockNames.StringSpawningTypes.inverselyProportional, 
    SPAWNING_TYPES.epsilon:blockNames.StringSpawningTypes.epsilon,
    SPAWNING_TYPES.FAST:blockNames.StringSpawningTypes.fast
}
