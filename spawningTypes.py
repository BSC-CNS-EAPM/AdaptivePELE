import blockNames

class SPAWNING_TYPES:
    sameWeight, inverselyProportional, epsilon, simulatedAnnealing, FAST = range(5)

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    SPAWNING_TYPES.sameWeight:blockNames.STRING_SPAWNING_TYPES.sameWeight, 
    SPAWNING_TYPES.inverselyProportional:blockNames.STRING_SPAWNING_TYPES.inverselyProportional, 
    SPAWNING_TYPES.epsilon:blockNames.STRING_SPAWNING_TYPES.epsilon,
    SPAWNING_TYPES.FAST:blockNames.STRING_SPAWNING_TYPES.FAST
}
