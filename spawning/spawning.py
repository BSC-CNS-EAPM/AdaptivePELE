import math
import numpy as np
from clustering import clustering
import random
import blockNames
import spawningTypes
import densitycalculator


def return_sign(i, m, n, r):
    """ Helper function, creates a three-piece step function"""
    if i <= n-m:
        return 1
    elif i <= r:
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

        densityBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityBuilder.build(spawningBlock)

        spawningTypeString = spawningBlock[blockNames.StringSpawningTypes.type]
        if spawningTypeString == blockNames.StringSpawningTypes.sameWeight:
            spawningCalculator = SameWeightDegeneracyCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.inverselyProportional:
            spawningCalculator = InverselyProportionalToPopulationCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.epsilon:
            spawningCalculator = EpsilonDegeneracyCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.fast:
            spawningCalculator = FASTDegeneracyCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.variableEpsilon:
            spawningCalculator = VariableEpsilonDegeneracyCalculator(densityCalculator)
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
        #Params general to all
        #Params specific to epsilon related spawning
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
        metrics = np.zeros(len(clusters))
        for i,cluster in enumerate(clusters):
            metrics[i] = cluster.getMetric()
        return metrics

    def getSizes(self, clusters):
        sizes = np.zeros(len(clusters))
        for i,cluster in enumerate(clusters):
            sizes[i] = cluster.elements
        return sizes

from abc import ABCMeta, abstractmethod
class DensitySpawningCalculator(SpawningCalculator):
    """
        Subclass of Spawning calculator, that ensures the definition of a density calculator.
    """
    def __init__(self, densityCalculator = densitycalculator.NullDensityCalculator()):
        SpawningCalculator.__init__(self)
        self.type = "BaseDensityClass" #change for abstract attribute
        self.densityCalculator = densityCalculator

    def calculateDensities(self, clusters):
        densities = np.zeros(len(clusters))
        for i,cluster in enumerate(clusters):
            contacts = cluster.getContacts()
            densities[i] = self.densityCalculator.calculate(contacts)
        return densities

class SameWeightDegeneracyCalculator(SpawningCalculator):
    def __init__(self):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.sameWeight

    #TODO: Don't back on pele, so that density calculator can be used
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

class InverselyProportionalToPopulationCalculator(DensitySpawningCalculator):
    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.inverselyProportional

    def log(self):
        pass

    def calculate(self, clusters, trajToDistribute, clusteringParams, currentEpoch=None):
        sizes = self.getSizes(clusters)
        densities = self.calculateDensities(clusters)

        weights = sizes/densities
        return self.divideInverselyProportionalToArray(sizes, trajToDistribute)

class EpsilonDegeneracyCalculator(DensitySpawningCalculator):
    """
        It uses epsilon * numTraj trajectories proportional to their energy and the rest inversely proportional to each cluster's population
    """
    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator(densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.epsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None

    #TODO add possibility for different pipes
    def log(self):
        if not self.degeneracyTotal is None:
            print "[SpawningLog] Total: %s" % str(self.degeneracyTotal)
        if not self.degeneracyInverselyProportional is None:
            print "[SpawningLog] Inversely prop: %s" % str(self.degeneracyInverselyProportional)
        if not self.degeneracyMetricProportional is None:
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
        if abs(shiftedMetrics.sum()) < 1e-8:
            weights = np.ones(len(metrics))/len(metrics)
        else:
            weights = (1.*shiftedMetrics)/sum(shiftedMetrics)

        """
        minimumValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, minimumValue)

        kbT = 0.001987*T
        weights = np.exp(-shiftedMetrics/kbT)
        weights /= sum(weights)
        """

        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

#TODO: change for composition rather than for inheritance
class VariableEpsilonDegeneracyCalculator(DensitySpawningCalculator):

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.epsilonDegeneracyCalculator = EpsilonDegeneracyCalculator(densityCalculator)
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

        rateEpsilonVariation = [(clusteringParams.maxEpsilon-clusteringParams.minEpsilon)/(middleWindow-leftWindow), (clusteringParams.maxEpsilon-clusteringParams.minEpsilon)/(clusteringParams.period-rightWindow-1)]

        clusteringParams.epsilon += return_sign(currentEpoch,leftWindow,middleWindow, rightWindow)*rateEpsilonVariation[currentEpoch>middleWindow]

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
        return self.epsilonDegeneracyCalculator.calculate(clusters, trajToDistribute, clusteringParams, currentEpoch)

class SimulatedAnnealingCalculator(SpawningCalculator):
    def __init__(self):
        SpawningCalculator.__init__(self)
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

class FASTDegeneracyCalculator(DensitySpawningCalculator):
    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.FAST
        self.densityCalculator = densityCalculator

    def normaliseArray(self, array):
        maxValue = float(np.max(array))
        minValue = float(np.min(array))

        if maxValue-minValue > 1e-5:
            normalisedArray = np.subtract(maxValue, array)/(maxValue - minValue)
        else:
            #all elements are equal
            n = len(array)
            normalisedSizes = np.ones(n)/float(n)
        return normalisedArray

    def calculateNormalisedSizes(self, clusters):
        sizes = self.getSizes(clusters)

        densities = self.calculateDensities(clusters)
        weightedSizes = sizes/densities

        return self.normaliseArray(weightedSizes)

    def calculateNormalisedMetrics(self, clusters):
        metrics = self.getMetrics(clusters)
        return self.normaliseArray(metrics)

    def calculate(self, clusters, trajToDivide, clusteringParams, currentEpoch=None):
        normalisedSizes = self.calculateNormalisedSizes(clusters)
        normalisedMetrics = self.calculateNormalisedMetrics(clusters)

        weight = normalisedSizes + 1.*normalisedMetrics

        return self.divideProportionalToArray(weight, trajToDivide)

    def log(self):
        pass

