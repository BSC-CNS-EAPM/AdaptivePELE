import math
import numpy as np
import clustering
import random
import blockNames
import spawningTypes


def return_sign(i, m, n):
    """ Helper function, creates a three-piece step function"""
    if i < n/2-m/2:
        return 1
    elif i <= n/2+m/2:
        return 0
    else:
        return -1

class SpawningAlgorithmBuilder:
    def build(self, spawningBlock):
        spawningCalculatorBuilder = SpawningBuilder()
        spawningCalculator = spawningCalculatorBuilder.buildSpawningCalculator(spawningBlock)

        spawningParams = SpawningParams()
        spawningParams.buildSpawningParameters(spawningBlock)

        return spawningCalculator, spawningParams

class SpawningBuilder:
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
        elif spawningTypeString == blockNames.StringSpawningTypes.variableEpsilon:
            spawningCalculator = VariableEpsilonDegeneracyCalculator()
        else:
            sys.exit("Unknown spawning type! Choices are: " + str(spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY.values()))
        return spawningCalculator


class SpawningParams:
    def __init__(self):
        self.epsilon = None
        self.temperature = None
        self.threshold = None
        self.reportFilename = None
        self.reportCol = None
        self.decrement = None
        self.varEpsilonType = None
        self.maxEpsilon = None
        self.minEpsilon = None
        self.variationWindow = None
        self.maxEpsilonWindow = None

    def buildSpawningParameters(self,spawningBlock):
        spawningParamsBlock = spawningBlock[blockNames.SpawningParams.params]
        spawningType = spawningBlock[blockNames.StringSpawningTypes.type]
        if spawningType == blockNames.StringSpawningTypes.epsilon or \
                spawningType == blockNames.StringSpawningTypes.variableEpsilon:
            self.epsilon = spawningParamsBlock[blockNames.SpawningParams.epsilon]
            self.reportFilename = spawningParamsBlock[blockNames.SpawningParams.report_filename]
            self.reportCol = spawningParamsBlock[blockNames.SpawningParams.report_col]
            self.temperature = spawningParamsBlock[blockNames.SpawningParams.temperature]
        if spawningType == blockNames.StringSpawningTypes.variableEpsilon:
            self.varEpsilonType = spawningParamsBlock[blockNames.SpawningParams.varEpsilonType]
            self.maxEpsilon = spawningParamsBlock[blockNames.SpawningParams.maxEpsilon]
            self.minEpsilon = spawningParamsBlock[blockNames.SpawningParams.minEpsilon]
            self.variationWindow = spawningParamsBlock[blockNames.SpawningParams.variationWindow]
            self.maxEpsilonWindow = spawningParamsBlock[blockNames.SpawningParams.maxEpsilonWindow]
            self.period = spawningParamsBlock.get(blockNames.SpawningParams.period, self.variationWindow) 
            self.period += np.sign(np.abs(self.variationWindow-self.period))
            # Add one epoch to the total lenght of the variation in the case of periodic variation to leave a step between variation periods


from abc import ABCMeta, abstractmethod
class SpawningCalculator:
    """
        The purpose of this abstract class is to contain the behaviour of the different strategies for the spawning.
        Spawning is the way in which we split the different explorers at the begining of each epoch.
    """
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


class SameWeightDegeneracyCalculator(SpawningCalculator):
    def __init__(self):
        self.spawningCalculator = SpawningCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.sameWeight

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

class InverselyProportionalToPopulationCalculator(SpawningCalculator):
    def __init__(self):
        self.spawningCalculator = SpawningCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.inverselyProportional

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
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.epsilon
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

class VariableEpsilonDegeneracyCalculator(EpsilonDegeneracyCalculator):

    def __init__(self):
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.variableEpsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None

    def linearVariation(self, clusteringParams, currentEpoch):
        if currentEpoch == 0:
            clusteringParams.epsilon = clusteringParams.minEpsilon
            return

        middleWindow = int(clusteringParams.period/2)
        leftWindow = int(clusteringParams.maxEpsilonWindow/2)
        rightWindow = leftWindow+middleWindow

        rateEpsilonVariation = [(clusteringParams.maxEpsilon-clusteringParams.minEpsilon)/(middleWindow-leftWindow-1), (clusteringParams.maxEpsilon-clusteringParams.minEpsilon)/(clusteringParams.period-rightWindow-1)]

        clusteringParams.epsilon += return_sign(currentEpoch,clusteringParams.maxEpsilonWindow,clusteringParams.period)*rateEpsilonVariation[currentEpoch>middleWindow]

    def calculateEpsilonValue(self, clusteringParams, currentEpoch):
        if currentEpoch is None or clusteringParams.variationWindow < currentEpoch:
            clusteringParams.epsilon = clusteringParams.minEpsilon
            return
        if clusteringParams.varEpsilonType == blockNames.VariableEpsilonTypes.linearVariation:
            self.linearVariation(clusteringParams,
                                 (currentEpoch%clusteringParams.period))
        else:
            sys.exit("Unknown epsilon variation type! Choices are: " +
                     str(spawningTypes.EPSILON_VARIATION_TYPE_TO_STRING_DICTIONARY.values()))

    def calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch=None):
        self.calculateEpsilonValue(clusteringParams, currentEpoch)
        return EpsilonDegeneracyCalculator.calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch)

class SimulatedAnnealingCalculator(SpawningCalculator):
    def __init__(self):
        self.spawningCalculator = SpawningCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.simulatedAnnealing

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

class FASTDegeneracyCalculator(SpawningCalculator):
    def __init__(self):
        self.type = spawningTypes.SPAWNING_TYPES.FAST
        self.spawningCalculator = SpawningCalculator()

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

