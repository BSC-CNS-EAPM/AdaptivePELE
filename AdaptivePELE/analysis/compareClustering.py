from AdaptivePELE.clustering import clustering, thresholdcalculator
import time
import sys
# import pickle
# import pdb as debug
import os
from AdaptivePELE.spawning import spawning, densitycalculator

nclusters = 80
thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
thresholdCalculator = thresholdCalculatorBuilder.build({
        "thresholdCalculator": {
            "type": "heaviside",
            "params": {
                "values": [2, 3, 4, 5],
                "conditions": [1, 0.75, 0.5]
            }
        }
})
# thresholdCalculator = thresholdCalculatorBuilder.build({})
# Distance index
# similarityEvaluator = clustering.differenceDistanceEvaluator()
# thresholdCalculatorAcc = thresholdCalculatorBuilder.build({
#     "thresholdCalculator": {
#         "type": "heaviside",
#         "params": {
#             "conditions": [1.2, 1.0, 0.5, 0.0],
#             "values": [0.2, 0.4, 0.5, 0.8]
#         }
#     }
# })
# Jaccard index
similarityEvaluator = clustering.JaccardEvaluator()
# thresholdCalculatorAcc = thresholdCalculatorBuilder.build({
#     "thresholdCalculator": {
#         "type" : "heaviside",
#         "params" : {
#             "conditions" : [1, 0.75, 0.5],
#             "values" : [0.025, 0.03, 0.04, 0.05]
#         }
#     }
# })
thresholdCalculatorAcc = thresholdCalculatorBuilder.build({
    "thresholdCalculator": {
        "type" : "heaviside",
        "params" : {
            "values" : [0.2, 0.3, 0.5, 0.8],
            "conditions" : [1.0, 0.75, 0.5]
        }
    }
})
# correlation
# similarityEvaluator = clustering.correlationEvaluator()
# thresholdCalculatorAcc = thresholdCalculatorBuilder.build({
#     "thresholdCalculator": {
#         "type": "constant",
#         "params": {
#             "value": 0.15
#         }
#     }
# })
densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
# densityCalculator = densityCalculatorBuilder.build({
#        "density": {
#             "type": "heaviside",
#             "params": {
#                  "conditions": [1.5, 1.0, 0.0],
#                  "values": [8.0, 2.37, 1.0, 0.5]
#              }
#         }
#    }
# )
# densityCalculator = densityCalculatorBuilder.build({})
densityCalculator = densityCalculatorBuilder.build({
    "density": {
        "type": "continuous"
    }
})
spawnParams = spawning.SpawningParams()
spawnParams.buildSpawningParameters({
    "type": "epsilon",
    "params": {
        "epsilon": 0.0,
        "T": 1000,
        "reportFilename": "report",
        "metricColumnInReport": 5
        }
    })
contactThresholdDistance = 8
resname = "STR"
nEpochs = 62
altSel = False
ntrajs = 64
ClDouble = clustering.ContactMapClustering(resname=resname,
                                           reportBaseFilename="report",
                                           columnOfReportFile=5,
                                           contactThresholdDistance=contactThresholdDistance)
ClAgg = clustering.ContactMapAgglomerativeClustering(nclusters, resname=resname,
                                                     reportBaseFilename="report",
                                                     columnOfReportFile=5,
                                                     contactThresholdDistance=contactThresholdDistance)
ClCont = clustering.ContactsClustering(thresholdCalculator, resname=resname,
                                       reportBaseFilename="report",
                                       columnOfReportFile=6,
                                       contactThresholdDistance=contactThresholdDistance,
                                       symmetries=[], altSelection=altSel)
ClAcc = clustering.ContactMapAccumulativeClustering(thresholdCalculatorAcc,
                                                    similarityEvaluator,
                                                    resname=resname,
                                                    reportBaseFilename="report",
                                                    columnOfReportFile=4,
                                                    contactThresholdDistance=contactThresholdDistance,
                                                    altSelection=altSel)
# spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
# spawningObject = spawning.UCBCalculator(densityCalculator)
spawningObject = spawning.EpsilonDegeneracyCalculator(densityCalculator)
processorMapping = [0 for i in xrange(ntrajs-1)]
for i in range(nEpochs):
    # path =["trajs/%d/run_traj*"%i]
    # paths_report = ["trajs/%d/run_report*"%i]
    # path = ["/home/bsc72/bsc72021/simulations/PR/testCM_4_32/simulation/PRprog_CM_variabExtra_UCB_5/%d/traj*" % i]
    # paths_report = ["/home/bsc72/bsc72021/simulations/PR/testCM_4_32/simulation/PRprog_CM_variabExtra_UCB_5//%d/report*" % i]
    # path = ["/gpfs/scratch/bsc72/bsc72021/AdaptiveCM/simulation/PRprog_4_64CMExtraSubset_prova_SASA3/%d/traj*"%i]
    # paths_report = ["/gpfs/scratch/bsc72/bsc72021/AdaptiveCM/simulation/PRprog_4_64CMExtraSubset_prova_SASA3/%d/report*"%i]
    path = ["/home/jgilaber/PR/PR_simulation_network/%d/traj*"%i]
    paths_report = ["/home/jgilaber/PR/PR_simulation_network/%d/report*"%i]
    # path = ["/home/jgilaber/4DAJ/4DAJ_4_64CMExtraSubsetSingle/%d/traj*"%i]
    # paths_report = ["/home/jgilaber/4DAJ/4DAJ_4_64CMExtraSubsetSingle/%d/report*"%i]
    trajs = clustering.getAllTrajectories(paths_report)
    total_snapshots = 0
    for traj in trajs:
        for line in open(traj, "r"):
            total_snapshots += 1
        total_snapshots -= 1
    sys.stderr.write("Total snapsthots for epoch %d: %d\n" % (i, total_snapshots))
    # startTimeCont = time.time()
    # ClCont.cluster(path, processorMapping)
    # endTimeCont = time.time()
    # sys.stderr.write("Total time of clustering contacts, epoch %d: %.6f\n"%(i,endTimeCont-startTimeCont))
    # sys.stderr.write("Number of clusters contacts epoch %d: %d\n"%(i,len(ClCont.clusters.clusters)))
    # degeneraciesCont = spawningObject.calculate(ClCont.clusters.clusters, ntrajs-1, spawnParams)
    # nProc = 0
    # clusterList = processorMapping[:]
    # for icl in xrange(len(ClCont.clusters.clusters)):
    #     for j in range(int(degeneraciesCont[icl])):
    #         clusterList[nProc] = icl
    #         nProc += 1
    # assert nProc == ntrajs-1
    # processorMapping = clusterList[1:]+[clusterList[0]]
    # with open("mappings/mapping%d.txt"%i, "w") as f:
    #     f.write(','.join(map(str, processorMapping)))
    # ClCont.writeOutput("clsummary",degeneraciesCont,"ClCont.pkl", False)
    # os.rename("clsummary/summary.txt", "results/summary_ClCont.txt")

    startTimeAcc = time.time()
    ClAcc.cluster(path, processorMapping)
    endTimeAcc = time.time()
    sys.stderr.write("Total time of clustering accumulative, epoch %d: %.6f\n" % (i,endTimeAcc-startTimeAcc))
    sys.stderr.write("Number of clusters accumulative epoch %d: %d\n" % (i,len(ClAcc.clusters.clusters)))
    degeneraciesAcc = spawningObject.calculate(ClAcc.clusters.clusters, ntrajs-1, spawnParams)
    # summaryFolder = "%d/clustering" % i
    # if not os.path.exists(summaryFolder):
    #     os.makedirs(summaryFolder)
    nProc = 0
    clusterList = processorMapping[:]
    for icl in xrange(len(ClAcc.clusters.clusters)):
        for j in range(int(degeneraciesAcc[icl])):
            clusterList[nProc] = icl
            nProc += 1
    assert nProc == ntrajs-1
    ClAcc.writeOutput("clsummary",degeneraciesAcc,"ClAcc.pkl", False)
    os.rename("clsummary/summary.txt", "results/summary_ClAcc.txt")
    processorMapping = clusterList[1:]+[clusterList[0]]
    with open("mappings/mapping%d.txt"%i, "w") as f:
        f.write(','.join(map(str, processorMapping)))
    # for i, element in enumerate(degeneraciesCont):
    #     if element > 2 and ClCont.clusters.clusters[i].elements > 1:
    #         for j in xrange(element):
    #             ClCont.clusters.clusters[i].writeSpawningStructure("/dev/null")
    # for i, cluster in enumerate(ClAcc.clusters.clusters):
    #     print "Cluster", i
    #     print cluster.elements
    #     print cluster.altStructure.altStructPQ

# ClCont.writeConformationNetwork("conformationNetwork.edgelist")
# ClCont.writeFDT("FDT.edgelist")
# ClCont.writeConformationNodeMetric("nodesMetrics.txt", 4)
# ClCont.writeConformationNodePopulation("nodesPopulation.txt")
# ClCont.writePathwayOptimalCluster("pathwayFDT.pdb")

# ClAcc.writeConformationNetwork("conformationNetwork.edgelist")
# ClAcc.writeFDT("FDT.edgelist")
# ClAcc.writeConformationNodeMetric("nodesMetrics.txt", 4)
# ClAcc.writeConformationNodePopulation("nodesPopulation.txt")
# ClAcc.writePathwayOptimalCluster("pathwayFDT.pdb")
