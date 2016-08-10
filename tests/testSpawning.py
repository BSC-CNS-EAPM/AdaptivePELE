import unittest
import numpy as np
import spawning.spawning as spawning
import clustering

class TestSpawningCalculator(unittest.TestCase):
    def testDivideTrajAccordingToWeights(self):
        spawningCalculator = spawning.SpawningCalculator() 
        weights1 = [0.5, 0.2, 0.2, 0.1]
        trajToDistribute1 = 12
        degeneracy1 = spawningCalculator.divideTrajAccordingToWeights(weights1, trajToDistribute1)

        golden1 = [6, 2, 3, 1]
        self.assertEqual(degeneracy1, golden1)

        falsegolden1 = [6, 2, 3, 2]
        self.assertNotEqual(degeneracy1, falsegolden1)


    def testDivideInverselyProportionalToArray(self):
        spawningCalculator = spawning.SpawningCalculator() 
        weights2 = np.array([0.5, 0.2, 0.2, 0.1])
        trajToDistribute2 = 12
        degeneracy2 = spawningCalculator.divideInverselyProportionalToArray(weights2, trajToDistribute2)

        golden2 = [1, 3, 3, 5]

        test2Passed = True
        self.assertEqual(degeneracy2, golden2)

        falsegolden2 = [1, 3, 3, 6]
        self.assertNotEqual(degeneracy2, falsegolden2)

    def testInverselyProportionalToPopulationCalculator(self):
        #test 3
        clusters = clustering.Clusters()
        sizes = [6,2,3,1]
        for size in sizes:
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            clusters.addCluster(cluster)

        inverselyProp = spawning.InverselyProportionalToPopulationCalculator()
        clusteringParams = None
        trajs = 10
        degeneracy3 = inverselyProp.calculate(clusters.clusters, trajs, clusteringParams)
        golden3 = [1,2,2,5]

        self.assertEqual(degeneracy3, golden3)

        falsegolden3 = [1,2,2,6]
        self.assertNotEqual(degeneracy3, falsegolden3)


    def testEpsilonCalculator(self):
        epsilon = spawning.EpsilonDegeneracyCalculator()
        params = spawning.SpawningParams()
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
        np.testing.assert_array_equal(degeneracy4, golden4)

        #test 5
    def testSameWeightDegeneracyCalculator(self):
        sameWeightDegCalculator = spawning.SameWeightDegeneracyCalculator()
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

        self.assertEqual(degeneracy5, golden5)

        falseGolden5 = [1, 1, 1, 2]
        self.assertNotEqual(degeneracy5, falseGolden5)

    def testVariableEpsilonCalculator(self):
        variable_epsilon = spawning.VariableEpsilonDegeneracyCalculator()
        params = spawning.SpawningParams()
        params.epsilon = 0.5
        params.varEpsilonType = "linearVariation"
        params.maxEpsilon = 0.75
        params.minEpsilon = 0.5
        params.variationWindow = 8
        params.maxEpsilonWindow = 2
        rateVariation= 0.25/2
        clusters = clustering.Clusters()
        sizes = [6,2,3,1]
        energies = [-4,-2,-2,-1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            cluster.metric = energy
            clusters.addCluster(cluster)

        trajs = 20
        
        degeneracy6 = variable_epsilon.calculate(clusters.clusters, trajs, params, 0)
        golden6 = np.array([7,4,4,5])
        np.testing.assert_array_equal(degeneracy6, golden6)
        self.assertAlmostEqual(params.epsilon,params.minEpsilon)
        degeneracy7 = variable_epsilon.calculate(clusters.clusters, trajs, params, 1)
        #TODO: check degeneracy after next steps
        self.assertAlmostEqual(params.epsilon, params.minEpsilon+rateVariation)
        degeneracy8 = variable_epsilon.calculate(clusters.clusters, trajs, params, 2)
        self.assertAlmostEqual(params.epsilon, params.maxEpsilon)
        degeneracy9 = variable_epsilon.calculate(clusters.clusters, trajs, params, 9)
        self.assertAlmostEqual(params.epsilon, params.minEpsilon)



def main():
    return unittest.main(exit=False)

if __name__ == '__main__':
    main()
