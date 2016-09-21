class ControlFileParams:
    generalParams = "GeneralParams"
    spawningBlockname = "SpawningParams"
    simulationBlockname = "SimulationParams"
    clusteringBlockname = "clusteringTypes"

class GeneralParams:
    params = {
        "restart" : "bool",
        "outputPath" : "unicode",
        "initialStructures" : "list",
        "debug" : "bool"
    }

class SpawningParams:
    params = {
        "epsilon" : "float",
        "T" : "int",
        "reportFilename" : "unicode",
        "metricColumnInReport" : "int",
        "varEpsilonType" : "unicode",
        "maxEpsilon" : "float",
        "minEpsilon": "float",
        "variationWindow" : "int",
        "maxEpsilonWindow" : "int",
        "period" : "int",
    }
    types = {
        "sameWeight" : {},
        "inverselyProportional" : {},
        "epsilon" : {
            "epsilon" : "float",
            "T" : "int",
            "reportFilename" : "unicode",
            "metricColumnInReport" : "int",
        },
        "FAST" : {
            "epsilon" : "float",
            "T" : "int",
            "reportFilename" : "unicode",
            "metricColumnInReport" : "int",
        },
        "variableEpsilon" : {
            "epsilon" : "float",
            "T" : "int",
            "reportFilename" : "unicode",
            "metricColumnInReport" : "int",
            "varEpsilonType" : "unicode",
            "maxEpsilon" : "float",
            "minEpsilon": "float",
            "variationWindow" : "int",
            "maxEpsilonWindow" : "int",
            "period" : "int",
        }
    }
    density = {
        "types" : {
            "heaviside" : "unicode",
            "null": "unicode",
            "constant" : "unicode"
        },
        "params" : {
            "heaviside" : "unicode",
            "null": "unicode",
            "constant" : "unicode",
            "values": "list",
            "conditions" : "list"
        }
    }

class SimulationParams:
    types = {"pele":{
                "processors" : "int",
                "controlFile" : "unicode",
                "seed" : "int",
                "peleSteps" : "int",
                "iterations" : "int"
    },
             "test":{
                "destination" : "unicode",
                "origin" :"unicode",
                "processors" : "int",
                "seed" : "int",
                "peleSteps" : "int",
                "iterations" : "int"
             },
             "md":{}}
    params = {
        "executable": "unicode",
        "data":"unicode",
        "documents": "unicode",
        "destination" : "unicode",
        "origin" :"unicode",
        "processors" : "int",
        "controlFile" : "unicode",
        "seed" : "int",
        "peleSteps" : "int",
        "iterations" : "int"
    }


class clusteringTypes:
    types = {
        "contacts" : {
            "ligandResname" : "unicode",
            "contactThresholdDistance" : "int",
        },
        "contactMapAffinity" : {
            "ligandResname" : "unicode",
            "contactThresholdDistance" : "int",
        },
        "contactMapAgglomerative" : {
            "ligandResname" : "unicode",
            "contactThresholdDistance" : "int",
            "nclusters" : "int"
        }
    }
    params = {
        "contacts" : "unicode",
        "contactMapAffinity" : "unicode",
        "contactMapAgglomerative" : "unicode",
        "ligandResname" : "unicode",
        "contactThresholdDistance" : "int",
        "nclusters" : "int"
    }
    thresholdCalculator = {
        "types" : {
            "heaviside": "unicode",
            "constant" : "unicode"
        },
        "params" : {
            "conditions" : "list",
            "values" : "list",
            "heaviside": "unicode",
            "constant" : "unicode"
        }
    }
