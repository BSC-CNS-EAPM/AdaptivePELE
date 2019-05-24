from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import shutil
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pyemma.msm as msm
import pyemma.coordinates as coor
from AdaptivePELE.freeEnergies import cluster
plt.style.use("ggplot")


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Plot information related to an MSM"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--lagtimes", type=int, nargs=2, required=True, help="Lagtimes to analyse, specifying lower and upper bound")
    parser.add_argument("-c", "--clusters", type=int, nargs=2, required=True, help="Number of clusters to analyse, specifying lower and upper bound")
    parser.add_argument("--lagtime_step", type=int, default=20, help="Resolution in the lagtimes to explore, default is 20")
    parser.add_argument("--cluster_step", type=int, default=20, help="Resolution in the clusters to explore, default is 20")
    parser.add_argument("-m", type=int, default=6, help="Number of eigenvalues to sum in the GMRQ")
    parser.add_argument("-n", type=int, default=5, help="Number of iterations for the cross validation")
    parser.add_argument("--skip_steps", type=int, default=0, help="Number of initial steps to skip")
    parser.add_argument("--tica", action="store_true", help="Whether to use TICA before clustering")
    parser.add_argument("--tica_lag", type=int, default=30, help="Lagtime for the TICA estimation")
    parser.add_argument("--out_path", type=str, default="", help="Path to store the output")
    args = parser.parse_args()
    return args.lagtimes, args.m, args.tica, args.tica_lag, args.out_path, args.clusters, args.lagtime_step, args.cluster_step, args.n, args.skip_steps


def lengthVsNtrajs(data, nruns, lagtime, clusters, outputFilename, cache, m, stride):
    nClusters = len(clusters)
    nLags = len(lagtime)
    results = np.zeros((nClusters, nLags))
    results_cv = np.zeros((nClusters, nLags))
    for i, cl in enumerate(clusters):
        clustering = coor.cluster_kmeans(data=data, k=cl, max_iter=500, stride=stride)
        for j, lag in enumerate(lagtime):
            if (cl, lag) in cache:
                print("Loading cached computation for %d clusters and %d lagtime" % (cl, lag))
                results[i][j], results_cv[i][j] = cache[(cl, lag)]
                with open(outputFilename, 'a') as f:
                    f.write("%d %d %f %f\n" % (cl, lag, results[i][j], results_cv[i][j]))
                continue
            print("Computing for %d clusters and %d lagtime" % (cl, lag))
            try:
                MSM = msm.estimate_markov_model(clustering.dtrajs, lag)
                print("MSM estimated on %d states" % MSM.nstates)
            except Exception:
                print("Estimation error in %d clusters, %d lagtime" % (cl, lag))
                results[i][j] = 0.0
                results_cv[i][j] = 0.0
                continue
            try:
                results[i][j] = np.mean(MSM.score(MSM.dtrajs_full, score_k=m))
            except Exception:
                print("Estimation error in %d clusters, %d lagtime" % (cl, lag))
                results[i][j] = 0.0
                results_cv[i][j] = 0.0
                continue
            try:
                results_cv[i][j] = np.mean(MSM.score_cv(MSM.dtrajs_full, score_k=m, n=nruns))
            except Exception:
                print("Estimation error in %d clusters, %d lagtime" % (cl, lag))
                results_cv[i][j] = 0.0

            with open(outputFilename, 'a') as f:
                f.write("%d %d %f %f\n" % (cl, lag, results[i][j], results_cv[i][j]))
    return results, results_cv


def saveResultsFileBckp(outputFilename):
    i = 1
    bckpFilename = outputFilename+".%d.bckp"
    while os.path.isfile(bckpFilename % i):
        i += 1
    try:
        shutil.copy(outputFilename, bckpFilename % i)
    except IOError:
        pass


def plotIsocostLines(extent, minCost, maxCost, steps=10):
    # d = (maxCost - minCost) / steps
    # costs = np.arange(minCost, maxCost, d)
    minCost = np.log10(minCost)
    maxCost = np.log10(maxCost)
    costs = np.logspace(minCost, maxCost, num=steps)
    for cost in costs:
        x = np.arange(extent[0], extent[1], 1)
        y = cost / x
        plt.plot(x, y, color="black")


def main(lagtimes, clusters, m, tica_lag, tica, output_path, lag_res, cl_res, nruns, skipFirstSnaphots):
    trajectoryFolder = "allTrajs"
    trajectoryBasename = "traj*"
    stride = 1
    if output_path and not os.path.exists(output_path):
        os.makedirs(output_path)
    scores_path = os.path.join(output_path, "scores")
    if not os.path.exists(scores_path):
        os.makedirs(scores_path)
    data, _ = cluster.loadTrajFiles(trajectoryFolder, trajectoryBasename)
    data = [da[skipFirstSnaphots:] for da in data]
    if tica:
        tica_obj = coor.tica(data, lag=tica_lag, var_cutoff=0.9, kinetic_map=True)
        print('TICA dimension ', tica_obj.dimension())
        data = tica_obj.get_output()

    cl_values = list(range(clusters[0], clusters[1], cl_res))
    if not cl_values:
        raise ValueError("Cluster values wrongly specified, ensure that the clusters are provided in the format lower bound - upper bound")
    lag_values = list(range(lagtimes[0], lagtimes[1], lag_res))
    if not lag_values:
        raise ValueError("Lagtime values wrongly specified, ensure that the clusters are provided in the format lower bound - upper bound")
    outputFilename = os.path.join(output_path, "results.txt")
    cache = {}
    if os.path.exists(outputFilename):
        with open(outputFilename) as fr:
            for line in fr:
                contents = line.rstrip().split()
                if not contents[0].isdigit():
                    continue
                cl, lag, score, score_cv = map(float, contents)
                cache[(int(cl), int(lag))] = (score, score_cv)

    saveResultsFileBckp(outputFilename)
    results, results_cv = lengthVsNtrajs(data, nruns, lag_values, cl_values, outputFilename, cache, m, stride)

    extent = [lag_values[0]-lag_res/2, lag_values[-1]-lag_res/2, cl_values[0]-cl_res/2, cl_values[-1]+cl_res/2]  # +1 for aesthetical purposes
    plt.figure(1)
    # plotIsocostLines(extent, (ilengths+dlengths)*(itrajs+dtrajs), (flengths-dlengths)*(ftrajs-dtrajs), 6)
    plt.imshow(results, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
    plt.title("Variational score using all data")
    plt.ylabel("Number of clusters")
    plt.xlabel("Lagtime")
    plt.colorbar()
    plt.savefig(os.path.join(output_path, "scores_training.png"))
    plt.figure(2)
    plt.imshow(results_cv, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.title("Variational score using cross validation data")
    plt.ylabel("Number of clusters")
    plt.xlabel("Lagtime")
    plt.savefig(os.path.join(output_path, "scores_crossvalidation.png"))
    plt.figure(3)
    plt.imshow(results, interpolation="bilinear", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.title("Variational score using all data")
    plt.ylabel("Number of clusters")
    plt.xlabel("Lagtime")
    plt.savefig(os.path.join(output_path, "scores_training_finer.png"))
    plt.figure(4)
    # plotIsocostLines(extent, (ilengths+dlengths)*(itrajs+dtrajs), (flengths-dlengths)*(ftrajs-dtrajs), 6)
    plt.imshow(results_cv, interpolation="bilinear", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.title("Variational score using cross validation data")
    plt.ylabel("Number of clusters")
    plt.xlabel("Lagtime")
    plt.savefig(os.path.join(output_path, "scores_crossvalidation_finer.png"))
    plt.show()


if __name__ == "__main__":
    lags, GMRQ, use_tica, lag_tica, out_path, cluster_list, lagtime_step, cluster_step, n, skipSteps = parse_arguments()
    main(lags, cluster_list, GMRQ, lag_tica, use_tica, out_path, lagtime_step, cluster_step, n, skipSteps)
