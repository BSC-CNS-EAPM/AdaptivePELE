from clustering import clustering
import unittest


class clusteringTest(unittest.TestCase):
    def testCluster(self):
        # preparation
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringParams = {"type": "contacts",
                            "params": {"ligandResname": "AIN",
                                       "contactThresholdDistance": 8}}
        clusteringInstance = clusteringBuilder.buildClustering(clusteringParams,
                                                               "ain_report", 3)

        trajNames = ["tests/data/aspirin_data/traj*"]

        # function to test
        clusteringInstance.cluster(trajNames)

        # assertion
        allClusters = clusteringInstance.clusters.clusters
        goldenNumberOfClusters = 2
        goldenEnergyCluster1 = -8424.8
        goldenEnergyCluster2 = -8453.29
        goldenElementsCluster1 = 2
        goldenElementsCluster2 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].getMetric(), goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].getMetric(), goldenEnergyCluster2, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)

    def test_cluster_sklearn_affinity(self):
        # preparation
        clusteringParams = {"type": "contactMapAffinity",
                            "params": {"ligandResname": "AIN",
                                       "contactThresholdDistance": 8}}
        clusteringBuilder = clustering.ClusteringBuilder()
        reportCol = 3
        clusteringInstance = clusteringBuilder.buildClustering(clusteringParams,
                                                               "ain_report", reportCol)

        trajNames = ["tests/data/aspirin_data/traj*"]

        # function to test
        clusteringInstance.cluster(trajNames)

        # assertion
        # Testing is diffrent now, previous assertions won't work
        # TODO: create new assertion for the new algorithm
        allClusters = clusteringInstance.clusters.clusters
        goldenNumberOfClusters = 2
        goldenEnergyCluster1 = -8424.8
        goldenEnergyCluster2 = -8453.29
        goldenElementsCluster1 = 2
        goldenElementsCluster2 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].metrics[reportCol], goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].metrics[reportCol], goldenEnergyCluster2, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)

    def test_cluster_sklearn_agglomerative(self):
        # preparation
        nclusters = 2
        clusteringParams = {"type": "contactMapAgglomerative",
                            "params": {"ligandResname": "AIN",
                                       "contactThresholdDistance": 8,
                                       "nclusters": nclusters}}
        clusteringBuilder = clustering.ClusteringBuilder()
        reportCol = 3
        clusteringInstance = clusteringBuilder.buildClustering(clusteringParams,
                                                               "ain_report", reportCol)

        trajNames = ["tests/data/aspirin_data/traj*"]

        # function to test
        clusteringInstance.cluster(trajNames)

        # assertion
        # Testing is diffrent now, previous assertions won't work
        # TODO: create new assertion for the new algorithm
        allClusters = clusteringInstance.clusters.clusters
        goldenNumberOfClusters = nclusters
        goldenEnergyCluster1 = -8424.8
        goldenEnergyCluster2 = -8453.29
        goldenElementsCluster1 = 2
        goldenElementsCluster2 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].metrics[reportCol], goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].metrics[reportCol], goldenEnergyCluster2, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)

    def test_cluster_accumulative(self):
        # preparation
        clusteringParams = {"type": "contactMapAccumulative",
                            "params": {"ligandResname": "AIN",
                                       "contactThresholdDistance": 8,
                                       "similarityEvaluator": "correlation"},
                            "thresholdCalculator": {
                                "type": "constant",
                                "params": {
                                    "value": 0.15
                                }
                            }}
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringInstance = clusteringBuilder.buildClustering(clusteringParams,
                                                               "ain_report", 3)

        trajNames = ["tests/data/aspirin_data/traj*"]

        # function to test
        clusteringInstance.cluster(trajNames)

        # assertion
        # Testing is diffrent now, previous assertions won't work
        # TODO: create new assertion for the new algorithm
        allClusters = clusteringInstance.clusters.clusters
        goldenNumberOfClusters = 3
        goldenEnergyCluster1 = -8421.5
        goldenEnergyCluster2 = -8424.8
        goldenEnergyCluster3 = -8453.29
        goldenElementsCluster1 = 1
        goldenElementsCluster2 = 1
        goldenElementsCluster3 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].getMetric(), goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].getMetric(), goldenEnergyCluster2, 2)
        self.assertAlmostEqual(allClusters[2].getMetric(), goldenEnergyCluster3, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)
        self.assertEqual(allClusters[2].elements, goldenElementsCluster3)
