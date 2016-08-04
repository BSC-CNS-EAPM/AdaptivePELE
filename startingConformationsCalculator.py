import math
import numpy as np
import clustering
import random
import blockNames
import spawningTypes

class StartingConformationBuilder:
    def buildSpawningCalculator(self, spawningBlock):
        spawningTypeString = spawningBlock[blockNames.STRING_SPAWNING_TYPES.type]
        if spawningTypeString == blockNames.STRING_SPAWNING_TYPES.sameWeight:
            spawningCalculator = SameWeightDegeneracyCalculator()
        elif spawningTypeString == blockNames.STRING_SPAWNING_TYPES.inverselyProportional:
            spawningCalculator = InverselyProportionalToPopulationCalculator()
        elif spawningTypeString == blockNames.STRING_SPAWNING_TYPES.epsilon:
            spawningCalculator = EpsilonDegeneracyCalculator()
        elif spawningTypeString == blockNames.STRING_SPAWNING_TYPES.FAST:
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
        spawningParamsBlock = spawningBlock[blockNames.SPAWNING_PARAMS.params]
        spawningType = spawningBlock[blockNames.STRING_SPAWNING_TYPES.type]
        if spawningType == spawningTypes.SPAWNING_TYPES.epsilon:
            spawningParams.epsilon = spawningParamsBlock[blockNames.SPAWNING_PARAMS.EPSILON]
            spawningParams.reportFilename = spawningParamsBlock[blockNames.SPAWNING_PARAMS.REPORT_FILENAME]
            spawningParams.reportCol = spawningParamsBlock[blockNames.SPAWNING_PARAMS.REPORT_COL]
            spawningParams.temperature = spawningParamsBlock[blockNames.SPAWNING_PARAMS.TEMPERATURE]

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

if __name__ == "__main__":
    """ TESTS, all in one place, to be moved """
    startingConfCalculator = StartingConformationsCalculator() 
    weights1 = [0.5, 0.2, 0.2, 0.1]
    trajToDistribute1 = 12
    degeneracy1 = startingConfCalculator.divideTrajAccordingToWeights(weights1, trajToDistribute1)

    golden1 = [6, 2, 3, 1]

    test1Passed = True
    if degeneracy1 != golden1: 
        print "ERROR IN TEST 1"
        test1Passed = False



    weights2 = np.array([0.5, 0.2, 0.2, 0.1])
    trajToDistribute2 = 12
    degeneracy2 = startingConfCalculator.divideInverselyProportionalToArray(weights2, trajToDistribute2)

    golden2 = [1, 3, 3, 5]

    test2Passed = True
    if degeneracy2 != golden2: 
        print "ERROR IN TEST 2"
        test2Passed = False

    #test 3
    clusters = clustering.Clusters()
    sizes = [6,2,3,1]
    for size in sizes:
        cluster = clustering.Cluster(None, None, None, None)
        cluster.elements = size
        clusters.addCluster(cluster)

    inverselyProp = InverselyProportionalToPopulationCalculator()
    clusteringParams = None
    trajs = 10
    degeneracy3 = inverselyProp.calculate(clusters.clusters, trajs, clusteringParams)
    golden3 = [1,2,2,5]

    test3Passed = True
    if degeneracy3 != golden3: 
        print "ERROR IN TEST 3"
        test3Passed = False


    #test 4
    epsilon = EpsilonDegeneracyCalculator()
    params = SpawningParams()
    params.epsilon = 0.5

    clusters = clustering.Clusters()
    sizes = [6,2,3,1]
    energies = [-4,-2,-2,-1]
    for size, energy in zip(sizes, energies):
        cluster = clustering.Cluster(None, None, None, None)
        cluster.elements = size
        cluster.metric = energy
        clusters.addCluster(cluster)

    trajs = 20
    degeneracy4 = epsilon.calculate(clusters.clusters, trajs, params)
    golden4 = np.array([7,4,4,5])

    test4Passed = True 
    if not np.array_equal(degeneracy4,golden4): 
        print "ERROR IN TEST 4"
        test4Passed = False


    #test 5
    sameWeightDegCalculator = SameWeightDegeneracyCalculator()
    params = None

    clusters = clustering.Clusters()
    sizes = [6,2,3,1]
    for size in sizes:
        cluster = clustering.Cluster(None, None, None, None)
        cluster.elements = size
        clusters.addCluster(cluster)

    trajs = 10
    degeneracy5 = sameWeightDegCalculator.calculate(clusters.clusters, trajs, params)
    golden5 = [1, 1, 1, 1]

    test5Passed = True
    if degeneracy5 != golden5: 
        print "ERROR IN TEST 5"
        test5Passed = False



    if test1Passed and test2Passed and test3Passed and test4Passed:
        print "All tests passed!"
    else:
        print "Not all tests passed!"
