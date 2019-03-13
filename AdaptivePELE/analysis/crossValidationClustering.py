from __future__ import absolute_import, division, print_function, unicode_literals
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pyemma.coordinates as coor
import pyemma.msm as msm
from AdaptivePELE.freeEnergies import cluster
plt.switch_backend("pdf")
plt.style.use("ggplot")


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Plot information related to an MSM"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--lagtimes", type=int, nargs="*", help="Lagtimes to analyse")
    parser.add_argument("-c", "--clusters", type=int, nargs="*", help="Number of clusters to analyse")
    parser.add_argument("-m", type=int, default=6, help="Number of eigenvalues to sum in the GMRQ")
    parser.add_argument("--tica", action="store_true", help="Whether to use TICA before clustering")
    parser.add_argument("--tica_lag", type=int, default=30, help="Lagtime for the TICA estimation")
    parser.add_argument("--out_path", type=str, default="", help="Path to store the output")
    args = parser.parse_args()
    return args.lagtimes, args.m, args.tica, args.tica_lag, args.out_path, args.clusters


def main(lagtimes, clusters, m, tica_lag, tica, output_path):
    trajectoryFolder = "allTrajs"
    trajectoryBasename = "traj*"
    stride = 1
    if output_path and not os.path.exists(output_path):
        os.makedirs(output_path)
    scores_path = os.path.join(output_path, "scores")
    if not os.path.exists(scores_path):
        os.makedirs(scores_path)
    data, _ = cluster.loadTrajFiles(trajectoryFolder, trajectoryBasename)
    if tica:
        tica_obj = coor.tica(data, lag=tica_lag, var_cutoff=0.9, kinetic_map=True)
        print('TICA dimension ', tica_obj.dimension())
        data = tica_obj.get_output()
    for tau in lagtimes:
        scores = []
        scores_cv = []
        print("Estimating MSM with %d lagtime" % tau)
        for k in clusters:
            print("Calculating scores with %d clusters" % k)
            # cluster data
            cl = coor.cluster_kmeans(data=data, k=k, max_iter=500, stride=stride)
            try:
                MSM = msm.estimate_markov_model(cl.dtrajs, tau)
                print("MSM estimated on %d states" % MSM.nstates)
            except Exception:
                print("Estimation error in %d clusters, %d lagtime" % (k, tau))
                scores.append(0)
                scores_cv.append(np.array([0, 0, 0, 0, 0]))
                continue
            try:
                scores.append(MSM.score(MSM.dtrajs_full, score_k=m))
            except Exception:
                print("Estimation error in %d clusters, %d lagtime" % (k, tau))
                scores.append(0)
                scores_cv.append(np.array([0, 0, 0, 0, 0]))
                continue
            try:
                scores_cv.append(MSM.score_cv(MSM.dtrajs_full, score_k=m, n=5))
            except Exception:
                print("Estimation error in %d clusters, %d lagtime" % (k, tau))
                scores_cv.append(np.array([0, 0, 0, 0, 0]))
        np.save(os.path.join(scores_path, "scores_lag_%d.npy" % tau), scores)
        np.save(os.path.join(scores_path, "scores_cv_lag_%d.npy" % tau), scores_cv)
        mean_scores = [sc.mean() for sc in scores_cv]
        std_scores = [sc.std() for sc in scores_cv]
        plt.figure()
        plt.plot(clusters, scores, label="Training")
        plt.errorbar(clusters, mean_scores, yerr=std_scores, fmt='k', label="Testing")
        plt.xlabel("Number of states")
        plt.ylabel("Score")
        plt.legend()
        plt.savefig(os.path.join(output_path, "scores_cv_lag_%d.png" % tau))

if __name__ == "__main__":
    lags, GMRQ, use_tica, lag_tica, out_path, cluster_list = parse_arguments()
    main(lags, cluster_list, GMRQ, lag_tica, use_tica, out_path)
