from __future__ import absolute_import, division, print_function, unicode_literals
import os
import sys
import time
import argparse
from builtins import range
from AdaptivePELE.utilities import utilities
from AdaptivePELE.clustering import clustering, thresholdcalculator
from AdaptivePELE.spawning import spawning, densitycalculator


def parseArguments():
    desc = "Create a new clustering with different parameters"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("path", type=str, help="Path to the simulation")
    parser.add_argument("resname", type=str, help="Resname in the pdb")
    parser.add_argument("nEpochs", type=int, help="Number of epochs to cluster")
    parser.add_argument("ntrajs", type=int, help="Number of trajectories per epoch")
    parser.add_argument("--altSel", action="store_true", help="Whether to use alternative selection")
    parser.add_argument("--writeClusters", action="store_true", help="Whether to write pdb structures for the clusters")
    args = parser.parse_args()
    return args.path, args.resname, args.nEpochs, args.ntrajs, args.altSel, args.writeClusters


if __name__ == "__main__":
    simulation_path, resname, nEpochs, ntrajs, altSel, writeClusters = parseArguments()
    contactThresholdDistance = 8
    altSel = False
    top_object = utilities.readClusteringObject("%s/topologies/topologies.pkl" % simulation_path)

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
    thresholdCalculator = thresholdCalculatorBuilder.build({
        "thresholdCalculator": {
            "type": "constant",
            "params": {
                "value": 1.0
            }
        }
    })
    # thresholdCalculator = thresholdCalculatorBuilder.build({})
    # Distance index
    # similarityEvaluator = clustering.CMSimilarityEvaluator("differenceDistance")
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
    similarityEvaluator = clustering.CMSimilarityEvaluator("Jaccard")
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
            "type": "heaviside",
            "params": {
                "values": [0.2, 0.3, 0.5, 0.8],
                "conditions": [1.0, 0.75, 0.5]
            }
        }
    })
    # correlation
    # similarityEvaluator = clustering.CMSimilarityEvaluator("correlation")
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
    densityCalculator = densityCalculatorBuilder.build({})
    densityCalculator = densityCalculatorBuilder.build({
        "density": {
            "type": "exitContinuous"
        }
    })
    spawnParams = spawning.SpawningParams()
    spawnParams.buildSpawningParameters({
        "type": "inverselyProportional",
        "params": {
            "epsilon": 0.0,
            "T": 1000,
            "reportFilename": "report",
            "metricColumnInReport": 3
            }
        })
    # spawnParams.buildSpawningParameters({
    #     "type": "epsilon",
    #     "params": {
    #         "epsilon": 0.0,
    #         "T": 1000,
    #         "reportFilename": "report",
    #         "metricColumnInReport": 4
    #         }
    #     })
    # spawnParams.buildSpawningParameters({
    #     "type": "REAP",
    #     "params": {
    #         "epsilon": 0.0,
    #         "T": 1000,
    #         "reportFilename": "report",
    #         "metricsInd": -1,
    #         "metricColumnInReport": 4
    #         }
    #     })

    ClCont = clustering.ContactsClustering(thresholdCalculator, resname=resname,
                                           reportBaseFilename="report",
                                           columnOfReportFile=4,
                                           contactThresholdDistance=contactThresholdDistance,
                                           symmetries=[], altSelection=altSel)

    ClAcc = clustering.ContactMapAccumulativeClustering(thresholdCalculatorAcc,
                                                        similarityEvaluator,
                                                        resname=resname,
                                                        reportBaseFilename="report",
                                                        columnOfReportFile=4,
                                                        contactThresholdDistance=contactThresholdDistance,
                                                        altSelection=altSel)
    spawningObject = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
    # spawningObject = spawning.UCBCalculator(densityCalculator)
    # spawningObject = spawning.EpsilonDegeneracyCalculator(densityCalculator)
    # spawningObject = spawning.REAPCalculator()
    # ClAcc.clusterInitialStructures(["/home/jgilaber/PR/PR_prog_initial_adaptive.pdb"])
    # ClCont.clusterInitialStructures(["/home/jgilaber/4DAJ/4DAJ_initial_adaptive.pdb"])
    # processorMapping = [0 for i in range(ntrajs-1)]
    if not os.path.exists("results"):
        os.makedirs("results")

    for i in range(nEpochs):
        path = ["%s/%d/traj*" % (simulation_path, i)]
        paths_report = ["%s/%d/report*" % (simulation_path, i)]
        trajs = clustering.getAllTrajectories(paths_report)
        total_snapshots = 0
        for traj in trajs:
            for line in open(traj, "r"):
                total_snapshots += 1
            total_snapshots -= 1
        sys.stderr.write("Total snapsthots for epoch %d: %d\n" % (i, total_snapshots))
        startTimeCont = time.time()
        # ClCont.cluster(path, processorMapping)
        ClCont.cluster(path, topology=top_object)
        endTimeCont = time.time()
        sys.stderr.write("Total time of clustering contacts, epoch %d: %.6f\n" % (i, endTimeCont-startTimeCont))
        sys.stderr.write("Number of clusters contacts epoch %d: %d\n" % (i, len(ClCont.clusters.clusters)))
        # invTrajs = (ntrajs - 1)/2
        # centTrajs = ntrajs - 1- invTrajs
        invTrajs = ntrajs-1
        degeneraciesCont = spawningObject.calculate(ClCont.clusters.clusters, invTrajs, spawnParams)
        spawningObject.log()
        nProc = 0
        clusterList = []
        for icl in range(len(ClCont.clusters.clusters)):
            for j in range(int(degeneraciesCont[icl])):
                clusterList.append(ClCont.clusters.clusters[icl].trajPosition)
                nProc += 1
        assert nProc == ntrajs-1
        processorMapping = clusterList[1:]+[clusterList[0]]
        if i < nEpochs - 1:
            with open("%d/processorMapping.txt" % (i+1), "w") as f:
                f.write(':'.join(map(str, processorMapping)))
        ClCont.writeOutput("results", degeneraciesCont, "results/ClCont.pkl", writeClusters)
        # os.rename("clsummary/summary.txt", "results/summary_ClCont.txt")
        # os.rmdir("clsummary")
