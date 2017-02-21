from ApativePELE.clustering import clustering, thresholdcalculator
import time
# import pickle
# import pdb as debug
import os
from AdaptivePELE.spawning import spawning, densitycalculator

ntrajs = 31
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
thresholdCalculatorAcc = thresholdCalculatorBuilder.build({
    "thresholdCalculator": {
        "type" : "constant",
        "params" : {
            "value" : 0.03
        }
    }
})
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
            "conditions" : [1.0, 0.75, 0.5],
            "values" : [0.4, 0.45, 0.5, 0.8]
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
densityCalculator = densityCalculatorBuilder.build({
       "density": {
            "type": "heaviside",
            "params": {
                 "conditions": [1.5, 1.0, 0.0],
                 "values": [8.0, 2.37, 1.0, 0.5]
             }
        }
   }
)
# densityCalculator = densityCalculatorBuilder.build({})
contactThresholdDistance = 8
resname = "STR"
nEpochs = 5
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
                                       columnOfReportFile=4,
                                       contactThresholdDistance=contactThresholdDistance,
                                       symmetries=[])
ClAcc = clustering.ContactMapAccumulativeClustering(thresholdCalculatorAcc,
                                                    similarityEvaluator,
                                                    resname=resname,
                                                    reportBaseFilename="report",
                                                    columnOfReportFile=5,
                                                    contactThresholdDistance=contactThresholdDistance)
# spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
spawningObject = spawning.UCBCalculator(densityCalculator)
for i in range(nEpochs):
    # path =["trajs/%d/run_traj*"%i]
    # paths_report = ["trajs/%d/run_report*"%i]
    path = ["/home/bsc72/bsc72021/simulations/PR/testCM_4_32/simulation/PRprog_CM_variabExtra_UCB_5/%d/traj*" % i]
    paths_report = ["/home/bsc72/bsc72021/simulations/PR/testCM_4_32/simulation/PRprog_CM_variabExtra_UCB_5//%d/report*" % i]
    # path = ["/gpfs/scratch/bsc72/bsc72021/AdaptiveCM/simulation/PRprog_4_64CMExtra/automate_CM_7/%d/traj*"%i]
    # paths_report = ["/gpfs/scratch/bsc72/bsc72021/AdaptiveCM/simulation/PRprog_4_64CMExtra/automate_CM_7/%d/report*"%i]
    # path = ["/home/jgilaber/PR/PRprog_4_64CMExtraSubset_thres/automate_CM_3/%d/traj*"%i]
    # paths_report = ["/home/jgilaber/PR/PRprog_4_64CMExtraSubset_thres/automate_CM_3/%d/report*"%i]
    # path = ["/home/jgilaber/4DAJ/4DAJ_4_64CMExtraSubsetSingle/%d/traj*"%i]
    # paths_report = ["/home/jgilaber/4DAJ/4DAJ_4_64CMExtraSubsetSingle/%d/report*"%i]
    trajs = clustering.getAllTrajectories(paths_report)
    total_snapshots = 0
    for traj in trajs:
        for line in open(traj, "r"):
            total_snapshots += 1
        total_snapshots -= 1
    print "Total snapsthots for epoch %d: %d" % (i, total_snapshots)
    # startTimeCont = time.time()
    # ClCont.cluster(path)
    # endTimeCont = time.time()
    # print "Total time of clustering contacts, epoch %d: %.6f"%(i,endTimeCont-startTimeCont)
    # print "Number of clusters contacts epoch %d: %d"%(i,len(ClCont.clusters.clusters))
    # degeneraciesCont = spawningObject.calculate(ClCont.clusters.clusters, ntrajs, {})
    # ClCont.writeOutput("clsummary",degeneraciesCont,"object_Cont.pkl", False)
    # os.rename("clsummary/summary.txt", "results/summary_ClCont.txt")
    # startTimeDouble = time.time()
    # ClDouble.cluster(path)
    # endTimeDouble = time.time()
    # print "Total time of clustering double, epoch %d: %.6f" % (i, endTimeDouble-startTimeDouble)
    # print "Number of clusters double epoch %d: %d" % (i, len(ClDouble.clusters.clusters))
    # degeneraciesDouble = spawningObject.calculate(ClDouble.clusters.clusters, ntrajs, {})
    # ClDouble.writeOutput("results",degeneraciesDouble,"ClDouble.pkl", False)
    # os.rename("results/summary.txt", "results/summmary_ClDouble.txt")
    # startTimeAgg = time.time()
    # ClAgg.cluster(path)
    # endTimeAgg = time.time()
    # print "Total time of clustering agglomerative, epoch %d: %.6f"%(i,endTimeAgg-startTimeAgg)
    # print "Number of clusters agglomerative epoch %d: %d"%(i,len(ClAgg.clusters.clusters))
    # print ""
    # degeneraciesAgg = spawningObject.calculate(ClAgg.clusters.clusters, ntrajs, {})
    # ClAgg.writeOutput("clsummary",degeneraciesAgg,"ClAgg.pkl", False)
    # os.rename("clsummary/summary.txt", "results/summary_ClAgg.txt")
    startTimeAcc = time.time()
    ClAcc.cluster(path)
    endTimeAcc = time.time()
    print "Total time of clustering accumulative, epoch %d: %.6f" % (i,endTimeAcc-startTimeAcc)
    print "Number of clusters accumulative epoch %d: %d" % (i,len(ClAcc.clusters.clusters))
    print ""
    degeneraciesAcc = spawningObject.calculate(ClAcc.clusters.clusters, ntrajs, {})
    summaryFolder = "%d/clustering" % i
    if not os.path.exists(summaryFolder):
        os.makedirs(summaryFolder)
    ClAcc.writeOutput(summaryFolder, degeneraciesAcc, "ClAcc.pkl", False)
    # os.remove(os.path.join(summaryFolder, "ClAcc.pkl"))
    # os.rename("clsummary/summary.txt", "%d/clustering/summary_ClAcc_set1.txt" % i)
