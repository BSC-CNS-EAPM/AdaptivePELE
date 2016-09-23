import blockNames

class SPAWNING_TYPES:
    sameWeight, inverselyProportional, epsilon, simulatedAnnealing, FAST, variableEpsilon = range(6)

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    SPAWNING_TYPES.sameWeight:blockNames.StringSpawningTypes.sameWeight,
    SPAWNING_TYPES.inverselyProportional:blockNames.StringSpawningTypes.inverselyProportional,
    SPAWNING_TYPES.epsilon:blockNames.StringSpawningTypes.epsilon,
    SPAWNING_TYPES.FAST:blockNames.StringSpawningTypes.fast,
    SPAWNING_TYPES.variableEpsilon:blockNames.StringSpawningTypes.variableEpsilon
    SPAWNING_TYPES.simulatedAnnealing:blockNames.StringSpawningTypes.simulatedAnnealing,
}
class EPSILON_VARIATION_TYPES:
    linearVariation, = range(1)

EPSILON_VARIATION_TYPE_TO_STRING_DICTIONARY = {
    EPSILON_VARIATION_TYPES.linearVariation:blockNames.VariableEpsilonTypes.linearVariation,
}
