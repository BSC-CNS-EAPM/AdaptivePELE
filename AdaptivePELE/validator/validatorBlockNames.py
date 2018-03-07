class ControlFileParams:
    generalParams = "GeneralParams"
    spawningBlockname = "SpawningParams"
    simulationBlockname = "SimulationParams"
    clusteringBlockname = "clusteringTypes"


class GeneralParams:
    mandatory = {
        "restart": "bool",
        "outputPath": "str",
        "initialStructures": "list"
    }
    params = {
        "restart": "bool",
        "outputPath": "str",
        "initialStructures": "list",
        "debug": "bool",
        "writeAllClusteringStructures": "bool",
        "nativeStructure": "str"
    }


class SpawningParams:
    params = {
        "epsilon": "numbers.Real",
        "T": "numbers.Real",
        "reportFilename": "str",
        "metricColumnInReport": "numbers.Real",
        "varEpsilonType": "str",
        "maxEpsilon": "numbers.Real",
        "minEpsilon": "numbers.Real",
        "variationWindow": "numbers.Real",
        "maxEpsilonWindow": "numbers.Real",
        "period": "numbers.Real",
        "alpha": "numbers.Real",
        "metricWeights": "str",
        "metricsInd": "list",
        "condition": "str",
        "n": "numbers.Real"
    }
    types = {
        "sameWeight": {
            "reportFilename": "str"
        },
        "independent": {
            "reportFilename": "str"
        },
        "inverselyProportional": {
            "reportFilename": "str"
        },
        "null": {
            "reportFilename": "str"
        },
        "epsilon": {
            "epsilon": "numbers.Real",
            "reportFilename": "str",
            "metricColumnInReport": "numbers.Real",
        },
        "FAST": {
            "epsilon": "numbers.Real",
            "reportFilename": "str",
            "metricColumnInReport": "numbers.Real",
        },
        "variableEpsilon": {
            "epsilon": "numbers.Real",
            "reportFilename": "str",
            "metricColumnInReport": "numbers.Real",
            "varEpsilonType": "str",
            "maxEpsilon": "numbers.Real"
        },
        "UCB": {
            "reportFilename": "str",
            "metricColumnInReport": "numbers.Real"
        },
        "REAP": {
            "reportFilename": "str",
            "metricColumnInReport": "numbers.Real"
        }
    }
    density = {
        "types": {
            "heaviside": "str",
            "null": "str",
            "constant": "str",
            "exitContinuous": "str",
            "continuous": "str"
        },
        "params": {
            "heaviside": "str",
            "null": "str",
            "constant": "str",
            "values": "list",
            "conditions": "list",
            "exitContinuous": "str",
            "continuous": "str"
        }
    }


class SimulationParams:
    types = {
        "pele": {
            "processors": "numbers.Real",
            "controlFile": "str",
            "seed": "numbers.Real",
            "peleSteps": "numbers.Real",
            "iterations": "numbers.Real"
            },
        "test": {
            "destination": "str",
            "origin": "str",
            "processors": "numbers.Real",
            "seed": "numbers.Real",
            "peleSteps": "numbers.Real",
            "iterations": "numbers.Real"
            },
        "md": {}}
    params = {
        "executable": "str",
        "data": "str",
        "documents": "str",
        "destination": "str",
        "origin": "str",
        "processors": "numbers.Real",
        "controlFile": "str",
        "seed": "numbers.Real",
        "peleSteps": "numbers.Real",
        "iterations": "numbers.Real",
        "modeMovingBox": "str",
        "boxCenter": "list",
        "boxRadius": "numbers.Real",
        "runEquilibration": "bool",
        "equilibrationMode": "str",
        "equilibrationLength": "numbers.Real",
        "numberEquilibrationStructures": "numbers.Real",
        "useSrun": "bool",
        "exitCondition": "dict"
    }
    exitCondition = {
        "types": {
            "metric": "str",
            "clustering": "str",
            "metricMultipleTrajectories": "str"
        },
        "params": {
            "metricCol": "int",
            "exitValue": "numbers.Real",
            "condition": "str",
            "numTrajs": "numbers.Real"
        }
    }


class clusteringTypes:
    types = {
        "rmsd": {
        },
        "contactMap": {
            "similarityEvaluator": "str"
        },
        "lastSnapshot": {
        }
    }
    params = {
        "rmsd": "str",
        "contactMap": "str",
        "contactThresholdDistance": "numbers.Real",
        "lastSnapshot": "str",
        "ligandResname": "str",
        "ligandResnum": "numbers.Real",
        "ligandChain": "str",
        "similarityEvaluator": "str",
        "symmetries": "list",
        "alternativeStructure": "bool",
        "nclusters": "numbers.Real"
    }
    thresholdCalculator = {
        "types": {
            "heaviside": "str",
            "constant": "str"
        },
        "params": {
            "conditions": "list",
            "values": "list",
            "value": "numbers.Real",
            "heaviside": "str",
            "constant": "str"
        }
    }
