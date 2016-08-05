import math
import numpy as np
import clustering
import random
import blockNames
import spawningTypes

class StartingConformationBuilder:
    def buildSpawningCalculator(self, spawningBlock):
        spawningTypeString = spawningBlock[blockNames.StringSpawningTypes.type]
        if spawningTypeString == blockNames.StringSpawningTypes.sameWeight:
            spawningCalculator = SameWeightDegeneracyCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.inverselyProportional:
            spawningCalculator = InverselyProportionalToPopulationCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.epsilon:
            spawningCalculator = EpsilonDegeneracyCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.fast:
            spawningCalculator = FASTDegeneracyCalculator()
        else:
            sys.exit("Unknown spawning type! Choices are: " + str(blockNames.SPAWNING_TYPE_TO_STRING_DICTIONARY.values()))
        return spawningCalculator


class SpawningParams:
    def __init__(self):
        self.epsilon = None
        self.temperature = None
        self.threshold = None
        self.reportFilename = None
        self.reportCol = None
        self.decrement = None

    def buildSpawningParameters(self,spawningBlock):
        spawningParamsBlock = spawningBlock[blockNames.SpawningParams.params]
        spawningType = spawningBlock[blockNames.StringSpawningTypes.type]
        if spawningType == blockNames.StringSpawningTypes.epsilon:
            self.epsilon = spawningParamsBlock[blockNames.SpawningParams.epsilon]
            self.reportFilename = spawningParamsBlock[blockNames.SpawningParams.report_filename]
            self.reportCol = spawningParamsBlock[blockNames.SpawningParams.report_col]
            self.temperature = spawningParamsBlock[blockNames.SpawningParams.temperature]

class SpawningAlgorithmBuilder:
    def build(self, spawningBlock):
        spawningCalculatorBuilder = StartingConformationBuilder()
        spawningCalculator = spawningCalculatorBuilder.buildSpawningCalculator(spawningBlock)

        spawningParams = SpawningParams()
        spawningParams.buildSpawningParameters(spawningBlock)

        return spawningCalculator, spawningParams

from abc import ABCMeta, abstractmethod
class StartingConformationsCalculator:
    def __init__(self):
        self.type = "BaseClass" #change for abstract attribute
        pass

    @abstractmethod
    def calculate(self, clusters, trajToDivide, clusteringParams, currentEpoch=None):
        pass

    @abstractmethod
    def log(self):
        pass

    def divideTrajAccordingToWeights(self, weights, trajToDistribute):
        """
            weights must be normalized (i.e. sum(weights) = 1)
        """
        degeneracy = []
        for i, weight in enumerate(weights):
            degeneracy.append(int(weight*trajToDistribute)) 

        #divide remaining traj to distribute according to decimal part
        decimalPart = []
        decimalPart = [math.modf(weight*trajToDistribute)[0] for weight in weights]
        sortedDecimals = np.argsort(decimalPart)
        sortedDecimals = sortedDecimals[::-1] #flip list

        leftProcessors = trajToDistribute-sum(degeneracy)
        for i in range(leftProcessors):
            degeneracy[sortedDecimals[i]] += 1

        return degeneracy

    def divideProportionalToArray(self, array, trajToDistribute):
        """
            Divides "trajToDistribute" proportionally to the values in "array"
        """
        if isinstance(array, list): array = np.array(array)
        weights = array/sum(array)
        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

    def divideInverselyProportionalToArray(self, array, trajToDistribute):
        """
            Divides "trajToDistribute" inversely proportional to the values in "array"
        """
        if isinstance(array, list): array = np.array(array) #it should always be true
        weights = 1./array
        weights /= sum(weights)
        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

    def getMetrics(self, clusters):
        metrics = []
        for cluster in clusters:
            metrics.append(cluster.getMetric())
        return metrics

    def getSizes(self, clusters):
        sizes = []
        for cluster in clusters:
            sizes.append(cluster.elements)
        return sizes


class SameWeightDegeneracyCalculator(StartingConformationsCalculator):
    def __init__(self):
        self.type = spawningTypes.SPAWNING_TYPES.sameWeight
        self.startingConformationsCalculator = StartingConformationsCalculator()

    def calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch=None):
        """
            We back on PELE to split equal trajectories
        """
        numClusters = len(clusters)
        trajToDistribute = min(trajToDistribute, numClusters)
        samples = random.sample(range(numClusters), trajToDistribute)
        degeneracy = [0] * len(clusters)
        for sample in samples:
            degeneracy[sample] = 1 
        return degeneracy

    def log(self):
        pass

class InverselyProportionalToPopulationCalculator(StartingConformationsCalculator):
    def __init__(self):
        self.type = spawningTypes.SPAWNING_TYPES.inverselyProportional
        self.startingConformationsCalculator = StartingConformationsCalculator()

    def log(self):
        pass

    def calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch=None):
        sizes = self.getSizes(clusters)
        return self.divideInverselyProportionalToArray(sizes, trajToDistribute)


class EpsilonDegeneracyCalculator(InverselyProportionalToPopulationCalculator):
    """
        It uses epsilon * numTraj trajectories proportional to their energy and the rest inversely proportional to each cluster's population
    """
    def __init__(self):
        self.type = spawningTypes.SPAWNING_TYPES.epsilon
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator()
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None

    #TODO add possibility for different pipes
    def log(self):
        if self.degeneracyTotal != None:
            print "[SpawningLog] Total: %s" % str(self.degeneracyTotal)
        if self.degeneracyInverselyProportional != None:
            print "[SpawningLog] Inversely prop: %s" % str(self.degeneracyInverselyProportional)
        if self.degeneracyMetricProportional != None:
            print "[SpawningLog] Metric prop:    %s" % str(self.degeneracyMetricProportional)

    def calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch=None):
        trajToEnergyProportional = int(clusteringParams.epsilon * trajToDistribute)
        trajToInverselyProportional = trajToDistribute - trajToEnergyProportional

        self.degeneracyInverselyProportional = self.inverselyProportionalCalculator.calculate(clusters, trajToInverselyProportional, clusteringParams)
        self.degeneracyMetricProportional = self.divideProcessorsMetricProportional(clusters, trajToEnergyProportional, clusteringParams.temperature)

        self.degeneracyTotal = np.array(self.degeneracyInverselyProportional) + np.array(self.degeneracyMetricProportional)
        return self.degeneracyTotal

    def divideProcessorsMetricProportional(self, clusters, trajToDistribute, T):
        metrics = self.getMetrics(clusters)
        if isinstance(metrics, list): metrics = np.array(metrics)

        """
            Shift so that differences become larger.
            Also, we can now merge positive & negative values
            Alternatives: Boltzmann weights
        """
        maximumValue = np.max(metrics)
        shiftedMetrics = np.subtract(metrics, maximumValue)

        #all shiftedMetrics <= 0, sum(shiftedMetrics) < 0 => weights >= 0
        weights = (1.*shiftedMetrics)/sum(shiftedMetrics)
        """


        minimumValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, minimumValue)

        kbT = 0.001987*T
        weights = np.exp(-shiftedMetrics/kbT)
        weights /= sum(weights)
        """



        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

class SimulatedAnnealingCalculator(StartingConformationsCalculator):
    def __init__(self):
        self.type = spawningTypes.SPAWNING_TYPES.simulatedAnnealing
        self.startingConformationsCalculator = StartingConformationsCalculator()

    def log(self):
        pass

    def computeTemperature(self, params, epoch):
        T = params.temperature - params.decrement*epoch
        if T < 300:
            return 300
        else:
            return T

    def calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch):
        metrics = getMetric()

        minimumValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, minimumValue)

        T = computeTemperature(clusteringParams, currentEpoch)
        kbT = 0.001987*T
        weights = np.exp(-shiftedMetrics/kbT)
        weights /= sum(weights)

        return self.divideProportionalToArray(weights, trajToDistribute)

class FASTDegeneracyCalculator(StartingConformationsCalculator):
    def __init__(self):
        self.type = spawningTypes.SPAWNING_TYPES.FAST
        self.startingConformationsCalculator = StartingConformationsCalculator()

    def calculate(self, clusters, trajToDivide, clusteringParams, currentEpoch=None):
        sizes = self.getSizes(clusters)
        maximumValue = float(np.max(sizes))
        minimumValue = float(np.min(sizes))
        normalisedSizes = np.subtract(maximumValue,sizes)
        if maximumValue-minimumValue > 1e-5:
            normalisedSizes /= (maximumValue - minimumValue)

        metrics = self.getMetrics(clusters)
        if isinstance(metrics, list): metrics = np.array(metrics)
        maximumValue = float(np.max(metrics))
        minimumValue = float(np.min(metrics))
        normalisedMetrics = np.subtract(maximumValue,metrics)
        if maximumValue-minimumValue > 1e-5:
            normalisedMetrics /= (maximumValue - minimumValue)

        weight = normalisedSizes + 4.*normalisedMetrics

        return self.divideProportionalToArray(weight, trajToDivide)

    def log(self):
        pass

