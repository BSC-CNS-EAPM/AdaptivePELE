from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import shutil
import argparse
import numpy as np
import pyemma.msm as msm
import pyemma.coordinates as coor
import matplotlib.pyplot as plt
from AdaptivePELE.freeEnergies import cluster, computeDeltaG as compute
plt.style.use('ggplot')


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Compute dG, stdDG for different lengths and number of trajectories"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--lagtime", type=int, default=100, help="Lagtimes to use, default 100")
    parser.add_argument("-c", "--cluster", type=int, default=100, help="Number of clusters to use, default 100")
    parser.add_argument("--lengths", type=int, nargs=2, required=True, help="Lengths to analyse, specifying lower and upper bound")
    parser.add_argument("--trajs", type=int, nargs=2, required=True, help="Number of trajectories to analyse, specifying lower and upper bound")
    parser.add_argument("--length_step", type=int, default=20, help="Resolution in the lengths to explore, default is 20")
    parser.add_argument("--trajs_step", type=int, default=20, help="Resolution in the trajs to explore, default is 20")
    parser.add_argument("-n", type=int, default=10, help="Number of iterations for the cross validation")
    parser.add_argument("--skip_steps", type=int, default=0, help="Number of initial steps to skip")
    parser.add_argument("--out_path", type=str, default="", help="Path to store the output")
    parser.add_argument("--cluster_each_iteration", action="store_true", help="Whether to cluster at each iteration, slower but more accurate results")
    args = parser.parse_args()
    return args.lagtime, args.out_path, args.cluster, args.length_step, args.trajs_step, args.n, args.skip_steps, args.lengths, args.trajs, args.cluster_each_iteration


def select_iteration_data(data, ntraj):
    return np.random.choice(range(len(data)), ntraj)


def estimateDG(data, nruns, cl, lag, ntraj, len_traj, skipFirstSnaphots, cluster_each_iteration):
    deltaG = []
    if not cluster_each_iteration:
        clustering = coor.cluster_kmeans(data=data, k=cl, max_iter=500, stride=1)
    for _ in range(nruns):
        data_it = select_iteration_data(data, ntraj)
        data_it = [data[j][skipFirstSnaphots:len_traj] for j in data_it]
        if cluster_each_iteration:
            clustering = coor.cluster_kmeans(data=data_it, k=cl, max_iter=500, stride=1)
            dtrajs = clustering.dtrajs
        else:
            dtrajs = clustering.assign(data_it)
        try:
            MSM = msm.estimate_markov_model(dtrajs, lag)
            print("MSM estimated on %d states" % MSM.nstates)
        except Exception:
            print("Estimation error in %d clusters, %d lagtime, %d trajectories of %d steps" % (cl, lag, ntraj, len_traj))
            continue
        pi, cl_centers = compute.ensure_connectivity(MSM, clustering.clustercenters)
        d = 0.75
        bins = compute.create_box(cl_centers, data_it, d)
        microstateVolume = compute.calculate_microstate_volumes_new(cl_centers, data_it, bins, d)
        _, string = compute.calculate_pmf(microstateVolume, pi)
        value = float(string.split()[1])
        deltaG.append(value)
    return np.mean(deltaG), np.std(deltaG)


def lengthVsNtrajs(data, nruns, lagtime, cl_num, lengths, ntrajs, outputFilename, cache, skipFirstSnaphots, cluster_each_iteration):
    nLengths = len(lengths)
    nNtrajs = len(ntrajs)
    results = np.zeros((nLengths, nNtrajs))
    stdDev = np.zeros((nLengths, nNtrajs))
    for i, length in enumerate(lengths):
        for j, ntraj in enumerate(ntrajs):
            if (length, ntraj) in cache:
                print("Loading cached computation for length:%d and ntrajs:%d" % (length, ntraj))
                results[i][j], stdDev[i][j] = cache[(length, ntraj)]
                with open(outputFilename, 'a') as f:
                    f.write("%d %d %f %f\n" % (length, ntraj, results[i][j], stdDev[i][j]))
                continue
            print("Computing for length:%d and ntrajs:%d" % (length, ntraj))
            results[i][j], stdDev[i][j] = estimateDG(data, nruns, cl_num, lagtime, ntraj, length, skipFirstSnaphots, cluster_each_iteration)
            with open(outputFilename, 'a') as f:
                f.write("%d %d %f %f\n" % (length, ntraj, results[i][j], stdDev[i][j]))
    return results, stdDev


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


def main(lagtime, clusters, output_path, dlengths, dtrajs, nruns, skipFirstSnaphots, lengths, trajs, cluster_each_iteration):
    ilengths, flengths = lengths
    lengths = list(range(ilengths, flengths+1, dlengths))
    itrajs, ftrajs = trajs
    ntrajs = list(range(itrajs, ftrajs+1, dtrajs))
    if not lengths:
        raise ValueError("Length values wrongly specified, ensure that the lengths are provided in the format lower bound - upper bound")
    if not ntrajs:
        raise ValueError("Trajectory values wrongly specified, ensure that the trajectory numbers are provided in the format lower bound - upper bound")
    trajectoryFolder = "allTrajs"
    trajectoryBasename = "traj*"
    data, _ = cluster.loadTrajFiles(trajectoryFolder, trajectoryBasename)
    data = [da[skipFirstSnaphots:] for da in data]

    cache = {}
    if output_path and not os.path.exists(output_path):
        os.makedirs(output_path)
    outputFilename = os.path.join(output_path, "results.txt")
    if os.path.exists(outputFilename):
        with open(outputFilename) as fr:
            for line in fr:
                contents = line.rstrip().split()
                if not contents[0].isdigit():
                    continue
                l, traj, dg, stdDG = map(float, contents)
                cache[(int(l), int(traj))] = (dg, stdDG)

    saveResultsFileBckp(outputFilename)

    with open(outputFilename, 'w') as f:
        f.write("Computing DG, stdDG for different lengths and number of trajectories\n")
        f.write("Lengths: %s\n" % lengths)
        f.write("Ntrajs: %s\n" % ntrajs)
        f.write("Skipping first: %d snapshots of each trajectory\n" % skipFirstSnaphots)
    results, stdDev = lengthVsNtrajs(data, nruns, lagtime, clusters, lengths, ntrajs, outputFilename, cache, skipFirstSnaphots, cluster_each_iteration)
    np.save("results.npy", results)
    np.save("stdDev.npy", stdDev)

    extent = [itrajs+1-dtrajs/2, ftrajs+dtrajs/2, ilengths-dlengths/2, flengths-dlengths/2]  # +1 for aesthetical purposes
    plt.figure(1)
    plotIsocostLines(extent, (ilengths+dlengths)*(itrajs+dtrajs), (flengths-dlengths)*(ftrajs-dtrajs), 6)
    plt.imshow(results, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
    # plt.imshow(results, interpolation="nearest", origin="lower", aspect="auto", extent=extent, vmin=-7, vmax=-5)
    plt.colorbar()
    plt.xlabel("Number of trajectories")
    plt.ylabel("Length of trajectories")
    plt.savefig(os.path.join(output_path, "dgGrid_%d_%d_rough.png" % (ftrajs, flengths)))
    plt.figure(2)
    plotIsocostLines(extent, (ilengths+dlengths)*(itrajs+dtrajs), (flengths-dlengths)*(ftrajs-dtrajs), 6)
    plt.imshow(results, interpolation="bilinear", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.xlabel("Number of trajectories")
    plt.ylabel("Length of trajectories")
    plt.savefig(os.path.join(output_path, "dgGrid_%d_%d_finer.png" % (ftrajs, flengths)))
    plt.show()


if __name__ == "__main__":
    lags, out_path, cl_val, lengths_step, traj_step, n, skipSteps, lengths_vals, trajs_vals, cl_iteration = parse_arguments()
    main(lags, cl_val, out_path, lengths_step, traj_step, n, skipSteps, lengths_vals, trajs_vals, cl_iteration)
