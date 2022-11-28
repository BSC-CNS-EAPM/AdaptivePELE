import os
import argparse
from AdaptivePELE.utilities import utilities
from AdaptivePELE.clustering import clustering
from AdaptivePELE.clustering import thresholdcalculator


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Cluster a simulation run"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("path", type=str, help="Path to the structures")
    parser.add_argument("resname", type=str, help="resname of the ligand")
    parser.add_argument("out_path", type=str, help="Path to the save the cluster structures")
    parser.add_argument("--radius", type=float, default=2, help="Size of the clusters")
    parser.add_argument("--top", type=str, default=None, help="Path to the topology")
    args = parser.parse_args()
    return args.path, args.resname, args.out_path, args.radius, args.top


def main(path_to_structures, resname, output_path, topology, cluster_size):
    thresholdCalculator = thresholdcalculator.ThresholdCalculatorConstant(cluster_size)
    topology_obj = utilities.getTopologyObject(topology)
    clusteringObj = clustering.ContactsClustering(thresholdCalculator, resname=resname)
    paths = [path_to_structures]
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    clusteringObj.cluster(paths, topology=topology_obj)
    for i, cl in enumerate(clusteringObj):
        cl.writePDB(os.path.join(output_path, "cluster_%d.pdb" % i))

if __name__ == "__main__":
    path, lig_resname, out, radius, top = parse_arguments()
    main(path, lig_resname, out, top, radius)
