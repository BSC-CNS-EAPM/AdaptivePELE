import matplotlib
matplotlib.use('Agg')
import sys
import os
import glob
from pylab import rcParams
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from AdaptivePELE.utilities import utilities
from AdaptivePELE.freeEnergies import cluster, extractCoords
from AdaptivePELE.analysis import splitTrajectory

def parseArgs():
    parser = argparse.ArgumentParser(description="Script that reclusters the Adaptive clusters")
    parser.add_argument('nClusters', type=int, help="Number of clusters")
    parser.add_argument('crit1', type=int, help="Metric1 to calculate from clusters")
    parser.add_argument('crit2', type=int, help="Metric2 to calculate from clusters")
    parser.add_argument("ligand_resname", type=str, help="Name of the ligand in the PDB")
    parser.add_argument("-atomId", nargs="*", default="", help="Atoms to use for the coordinates of the conformation, if not specified use the center of mass")
    parser.add_argument('--o', type=str, help="Output folder", default="Cluster_analisis")
    parser.add_argument('--top', type=str, help="Topology file", default=None)
    parser.add_argument('--cpus', type=int, help="Cpus to use", default=1)
    parser.add_argument('--report', type=str, help="Report filenames i.e. run_report_", default="report_")
    parser.add_argument('--traj', type=str, help="Trajectory filenames i.e. run_trajectory_", default="trajectory_")
    parser.add_argument('--use_pdb', action="store_true", help="To use when having pdb files with .xtc extension")
    args = parser.parse_args()
    return args.nClusters, args.crit1, args.crit2, args.ligand_resname, args.atomId, args.o, args.top, args.cpus, args.report, args.traj, args.use_pdb


def plotClusters(fields1, fields2, crit1, crit2, output):
    labels = ["cluster_{}".format(i) for i in np.arange(len(fields1))]
    fig, ax = plt.subplots()
    ax.scatter(fields1, fields2, label=labels)
    ax.set_title('RMSD Cluster {} vs {}'.format(crit1, crit2))
    ax.set_xlabel(crit1)
    ax.set_ylabel(crit2)
    print("Plotting")
    fig.savefig(os.path.join(output, "ClusterMap.pdf"), format='pdf', dpi=1200)

def writePDB(pmf_xyzg, title="clusters.pdb"):
    templateLine = "HETATM%s  H%sCLT L 502    %s%s%s  0.75%s           H\n"

    content = ""
    for j, line in enumerate(pmf_xyzg):
        number = str(j).rjust(5)
        number3 = str(j).ljust(3)
        x = ("%.3f" % line[0]).rjust(8)
        y = ("%.3f" % line[1]).rjust(8)
        z = ("%.3f" % line[2]).rjust(8)
        g = 0
        content += templateLine % (number, number3, x, y, z, g)

    with open(title, 'w') as f:
        f.write(content)


def writeInitialStructures(field1, field2, crit1, crit2, centers_info, filename_template, traj, topology=None, use_pdb=False):
    for cluster_num, field1, field2 in zip(centers_info, field1, field2):
        epoch_num, traj_num, snap_num = map(int, centers_info[cluster_num]['structure'])
        trajectory = "{}/{}{}.xtc".format(epoch_num, traj, traj_num) if topology else "{}/{}{}.pdb".format(epoch_num, traj, traj_num)
        snapshots = utilities.getSnapshots(trajectory, topology=topology, use_pdb=use_pdb)
        filename = filename_template.format(cluster_num, crit1, field1, crit2, field2)
        if not topology:
            with open(filename, "w") as fw:
                fw.write(snapshots[snap_num])
        else:
            splitTrajectory.main("", [trajectory, ], topology, [snap_num+1,],template=filename, use_pdb=use_pdb)


def get_centers_info(trajectoryFolder, trajectoryBasename, num_clusters, clusterCenters):
    centersInfo = {x: {"structure": None, "minDist": 1e6, "center": None} for x in xrange(num_clusters)}

    trajFiles = glob.glob(os.path.join(trajectoryFolder, trajectoryBasename))
    for traj in trajFiles:
        _, epoch, iTraj = os.path.splitext(traj)[0].split("_", 3)
        trajCoords = np.loadtxt(traj)
        if len(trajCoords.shape) < 2:
            trajCoords = [trajCoords]
        for snapshot in trajCoords:
            nSnap = snapshot[0]
            snapshotCoords = snapshot[1:]
            dist = np.sqrt(np.sum((clusterCenters-snapshotCoords)**2, axis=1))
            for clusterInd in xrange(num_clusters):
                if dist[clusterInd] < centersInfo[clusterInd]['minDist']:
                    centersInfo[clusterInd]['minDist'] = dist[clusterInd]
                    centersInfo[clusterInd]['structure'] = (epoch, int(iTraj), nSnap)
                    centersInfo[clusterInd]['center'] = snapshotCoords
    return centersInfo


def get_metric(criteria, epoch_num, traj_num, snap_num, report):
    report = os.path.join(str(epoch_num), "{}{}".format(report, traj_num))
    report_data = pd.read_csv(report, sep='    ', engine='python')
    value = report_data.iloc[snap_num, criteria-1]
    header = list(report_data)[criteria-1]
    return value, header


def main(num_clusters, criteria1, criteria2, output_folder, ligand_resname, atom_ids, cpus=2, topology=None, report="report_", traj="trajectory_", use_pdb=False):
    if not glob.glob("*/extractedCoordinates/coord_*"):
    	extractCoords.main(lig_resname=ligand_resname, non_Repeat=True, atom_Ids=atom_ids, nProcessors=cpus, parallelize=False, topology=topology, use_pdb=use_pdb)
    trajectoryFolder = "allTrajs"
    trajectoryBasename = "*traj*"
    stride = 1
    clusterCountsThreshold = 0
    folders = utilities.get_epoch_folders(".")
    folders.sort(key=int)

    clusteringObject = cluster.Cluster(num_clusters, trajectoryFolder,
                                       trajectoryBasename, alwaysCluster=True,
                                       stride=stride)
    clusteringObject.clusterTrajectories()
    clusteringObject.eliminateLowPopulatedClusters(clusterCountsThreshold)
    clusterCenters = clusteringObject.clusterCenters
    centersInfo = get_centers_info(trajectoryFolder, trajectoryBasename, num_clusters, clusterCenters)
    COMArray = [centersInfo[i]['center'] for i in xrange(num_clusters)]

    fields1 = []
    fields2 = []
    for cluster_num in centersInfo:
        epoch_num, traj_num, snap_num = map(int, centersInfo[cluster_num]['structure'])
        field1, crit1_name = get_metric(criteria1, epoch_num, traj_num, snap_num, report)
        field2, crit2_name = get_metric(criteria2, epoch_num, traj_num, snap_num, report)
	fields1.append(field1)	
	fields2.append(field2)	

    if output_folder is not None:
        outputFolder = os.path.join(output_folder, "")
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)
    else:
        outputFolder = ""
    writePDB(COMArray, outputFolder+"clusters_%d_KMeans_allSnapshots.pdb" % num_clusters)
    writeInitialStructures(fields1, fields2, crit1_name, crit2_name, centersInfo, outputFolder+"cluster_{}_{}_{}_{}_{}.pdb", traj, topology=topology, use_pdb=use_pdb) 
    plotClusters(fields1, fields2, crit1_name, crit2_name, outputFolder)

if __name__ == "__main__":
    n_clusters, criteria1, criteria2, lig_name, atom_id, output, top, cpus, report, traj, use_pdb = parseArgs()
    main(n_clusters, criteria1, criteria2, output, lig_name, atom_id, cpus, top, report, traj, use_pdb)
