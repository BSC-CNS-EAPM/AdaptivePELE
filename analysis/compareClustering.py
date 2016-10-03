from clustering import clustering, thresholdcalculator
import time
import pickle
import pdb as debug
from spawning import spawning, densitycalculator

ntrajs = 63
nclusters = 80
thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
#thresholdCalculator = thresholdCalculatorBuilder.build({
#        "thresholdCalculator" : {
#            "type" : "constant",
#            "params" : {
#                "value" : 2
#            }
#        }
#})
thresholdCalculator = thresholdCalculatorBuilder.build({})
densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
densityCalculator = densityCalculatorBuilder.build({
       "density" : {
            "type" : "heaviside",
            "params" : {
                 "conditions" : [1.5, 1.0, 0.0],
                 "values" : [8.0, 2.37, 1.0, 0.5]
             }
        }
   }
)
# densityCalculator = densityCalculatorBuilder.build({})
contactThresholdDistance = 8
ClDouble = clustering.ContactMapClustering(resname="ALJ",
                                           reportBaseFilename="report",
                                           columnOfReportFile=5,
                                           contactThresholdDistance=contactThresholdDistance)
ClAgg = clustering.ContactMapAgglomerativeClustering(nclusters, resname="ALJ",
                                                     reportBaseFilename="report",
                                                     columnOfReportFile=5,
                                                     contactThresholdDistance=contactThresholdDistance)
ClCont = clustering.ContactsClustering(thresholdCalculator, resname="ALJ",
                                       reportBaseFilename="report",
                                       columnOfReportFile=5,
                                       contactThresholdDistance=contactThresholdDistance)
spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)

for i in range(1):
    # path=["trajs/%d/run_traj*"%i]
    # paths_report=["trajs/%d/run_report*"%i]
    # path=["../3ptb_4_64_inversely_1/%d/traj*" % i]
    # paths_report=["../3ptb_4_64_inversely_1/%d/report*" % i]
    path=["/home/bsc72/bsc72021/simulations/5ALJ/5ALJ_evp_agg_linear_80clusters/simulations/5ALJ_evp_agg/%d/traj*"%i]
    paths_report=["/home/bsc72/bsc72021/simulations/5ALJ/5ALJ_evp_agg_linear_80clusters/simulations/5ALJ_evp_agg/%d/report*"%i]
    trajs=clustering.getAllTrajectories(paths_report)
    total_snapshots = 0
    for traj in trajs:
        for line in open(traj, "r"):
            total_snapshots += 1
        total_snapshots -= 1
    print "Total snapsthots for epoch %d: %d" % (i, total_snapshots)
    startTimeCont = time.time()
    ClCont.cluster(path)
    endTimeCont = time.time()
    print "Total time of clustering contacts, epoch %d: %.6f"%(i,endTimeCont-startTimeCont)
    print "Number of clusters contacts epoch %d: %d"%(i,len(ClCont.clusters.clusters))
    degeneraciesCont = spawningObject.calculate(ClCont.clusters.clusters, ntrajs, {})
    ClCont.writeOutput("results",degeneraciesCont,"bla.pkl", False)
    # for index in range(len(degeneracies)):
    #     #print index, ClCont.clusters.clusters[index].elements, degeneracies[index]
    #     ClCont.clusters.clusters[index].metric = degeneracies[index]
    #startTimeDouble = time.time()
    #ClDouble.cluster(path)
    #endTimeDouble = time.time()
    #print "Total time of clustering double, epoch %d: %.6f" % (i, endTimeDouble-startTimeDouble)
    #print "Number of clusters double epoch %d: %d" % (i, len(ClDouble.clusters.clusters))
    startTimeAgg = time.time()
    ClAgg.cluster(path)
    endTimeAgg = time.time()
    print "Total time of clustering agglomerative, epoch %d: %.6f"%(i,endTimeAgg-startTimeAgg)
    print "Number of clusters agglomerative epoch %d: %d"%(i,len(ClAgg.clusters.clusters))
    print ""
    degeneraciesAgg = spawningObject.calculate(ClAgg.clusters.clusters, ntrajs, {})
    ClAgg.writeOutput("results",degeneraciesAgg,"bla.pkl", False)

with open("ClCont.pkl", "w") as f:
    pickle.dump(ClCont.clusters, f, -1)

#with open("ClDouble.pkl", "w") as f:
#    pickle.dump(ClDouble.clusters, f, -1)

with open("ClAgg.pkl", "w") as f:
    pickle.dump(ClAgg.clusters, f, -1)
