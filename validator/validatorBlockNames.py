class ControlFileParams:
    generalParams = "GeneralParams"
    spawningBlockname = "SpawningParams"
    simulationBlockname = "SimulationParams"
    clusteringBlockname = "clusteringTypes"


class GeneralParams:
    mandatory = {
        "restart": "bool",
        "outputPath": "unicode",
        "initialStructures": "list",
        "debug": "bool",
        "writeAllClusteringStructures": "bool"
    }
    params = {
        "restart": "bool",
        "outputPath": "unicode",
        "initialStructures": "list",
        "debug": "bool",
        "writeAllClusteringStructures": "bool",
        "nativeStructure": "unicode"
    }


class SpawningParams:
    params = {
        "epsilon": "numbers.Real",
        "T": "numbers.Real",
        "reportFilename": "unicode",
        "metricColumnInReport": "numbers.Real",
        "varEpsilonType": "unicode",
        "maxEpsilon": "numbers.Real",
        "minEpsilon": "numbers.Real",
        "variationWindow": "numbers.Real",
        "maxEpsilonWindow": "numbers.Real",
        "period": "numbers.Real",
        "metricWeights": "unicode"
    }
    types = {
        "sameWeight": {},
        "inverselyProportional": {},
        "epsilon": {
            "epsilon": "numbers.Real",
            "T": "numbers.Real",
            "reportFilename": "unicode",
            "metricColumnInReport": "numbers.Real",
        },
        "FAST": {
            "epsilon": "numbers.Real",
            "T": "numbers.Real",
            "reportFilename": "unicode",
            "metricColumnInReport": "numbers.Real",
        },
        "variableEpsilon": {
            "epsilon": "numbers.Real",
            "T": "numbers.Real",
            "reportFilename": "unicode",
            "metricColumnInReport": "numbers.Real",
            "varEpsilonType": "unicode",
            "maxEpsilon": "numbers.Real",
            "minEpsilon": "numbers.Real",
            "variationWindow": "numbers.Real",
            "maxEpsilonWindow": "numbers.Real",
            "period": "numbers.Real",
        }
    }
    density = {
        "types": {
            "heaviside": "unicode",
            "null": "unicode",
            "constant": "unicode",
            "continuous": "unicode"
        },
        "params": {
            "heaviside": "unicode",
            "null": "unicode",
            "constant": "unicode",
            "values": "list",
            "conditions": "list",
            "continuous": "unicode"
        }
    }


class SimulationParams:
    types = {"pele": {
                "processors": "numbers.Real",
                "controlFile": "unicode",
                "seed": "numbers.Real",
                "peleSteps": "numbers.Real",
                "iterations": "numbers.Real"
    },
             "test": {
                "destination": "unicode",
                "origin": "unicode",
                "processors": "numbers.Real",
                "seed": "numbers.Real",
                "peleSteps": "numbers.Real",
                "iterations": "numbers.Real"
             },
             "md": {}}
    params = {
        "executable": "unicode",
        "data": "unicode",
        "documents": "unicode",
        "destination": "unicode",
        "origin": "unicode",
        "processors": "numbers.Real",
        "controlFile": "unicode",
        "seed": "numbers.Real",
        "peleSteps": "numbers.Real",
        "iterations": "numbers.Real",
        "exitCondition": "dict"
    }
    exitCondition = {
        "types": {
            "metric": "unicode"
        },
        "params": {
            "metricCol": "int",
            "exitValue": "numbers.Real"
        }
    }


class clusteringTypes:
    types = {
        "contacts": {
            "ligandResname": "unicode",
            "contactThresholdDistance": "numbers.Real",
        },
        "contactMapAffinity": {
            "ligandResname": "unicode",
            "contactThresholdDistance": "numbers.Real",
        },
        "contactMapAgglomerative": {
            "ligandResname": "unicode",
            "contactThresholdDistance": "numbers.Real",
            "nclusters": "numbers.Real"
        },
        "contactMapAccumulative": {
            "ligandResname": "unicode",
            "contactThresholdDistance": "numbers.Real",
            "similarityEvaluator": "unicode"
        },
        "lastSnapshot": {
            "ligandResname": "unicode",
            "contactThresholdDistance": "numbers.Real",
        }
    }
    params = {
        "contacts": "unicode",
        "contactMapAffinity": "unicode",
        "contactMapAgglomerative": "unicode",
        "contactMapAccumaltive": "unicode",
        "lastSnapshot": "unicode",
        "ligandResname": "unicode",
        "contactThresholdDistance": "numbers.Real",
        "similarityEvaluator": "unicode",
        "symmetries": "list",
        "nclusters": "numbers.Real"
    }
    thresholdCalculator = {
        "types": {
            "heaviside": "unicode",
            "constant": "unicode"
        },
        "params": {
            "conditions": "list",
            "values": "list",
            "value": "numbers.Real",
            "heaviside": "unicode",
            "constant": "unicode"
        }
    }
