import clustering
import unittest

class clusteringTest(unittest.TestCase):
    def test_cluster(self):
        #preparation
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringInstance = clusteringBuilder.buildClustering("contacts","AIN", "ain_report", 3)

        trajNames = ["tests/data/aspirin_data/traj*"]

        #function to test
        clusteringInstance.cluster(trajNames)

        #assertion
        allClusters = clusteringInstance.clusters.clusters
        goldenNumberOfClusters = 2
        goldenEnergyCluster1 = -8424.8
        goldenEnergyCluster2 = -8453.29
        goldenElementsCluster1 = 2
        goldenElementsCluster2 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].metric, goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].metric, goldenEnergyCluster2, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)

    def test_cluster_sklearn(self):
         #preparation
         clusteringBuilder = clustering.ClusteringBuilder()
         clusteringInstance = clusteringBuilder.buildClustering("contactmap","AIN", "ain_report", 3)

         trajNames = ["tests/data/aspirin_data/traj*"]

         #function to test
         clusteringInstance.cluster(trajNames)

         #assertion
         # Testing is diffrent now, previous assertions won't work
         # TODO: create new assertion for the new algorithm
         allClusters = clusteringInstance.clusters.clusters
         clusteringInstance.clusters.printClusters()
         goldenNumberOfClusters = 2
         goldenEnergyCluster1 = -8424.8
         goldenEnergyCluster2 = -8453.29
         goldenElementsCluster1 = 2
         goldenElementsCluster2 = 1

         self.assertEqual(len(allClusters), goldenNumberOfClusters)
         self.assertAlmostEqual(allClusters[0].metric, goldenEnergyCluster1, 2)
         self.assertAlmostEqual(allClusters[1].metric, goldenEnergyCluster2, 2)
         self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
         self.assertEqual(allClusters[1].elements, goldenElementsCluster2)
