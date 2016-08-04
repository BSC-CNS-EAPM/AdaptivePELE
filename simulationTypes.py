import blockNames

class SIMULATION_TYPE:
    PELE, MD, TEST = range(3)

SIMULATION_TYPE_TO_STRING_DICTIONARY = {
    SIMULATION_TYPE.PELE:blockNames.SimulationType.pele,
    SIMULATION_TYPE.MD:blockNames.SimulationType.md,
    SIMULATION_TYPE.TEST:blockNames.SimulationType.test
}
