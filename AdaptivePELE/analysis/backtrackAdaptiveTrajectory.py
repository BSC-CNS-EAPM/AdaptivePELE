"""
    Recreate the trajectory fragments to the led to the discovery of a snapshot,
    specified by the tuple (epoch, trajectory, snapshot) and write as a pdb file
"""
from __future__ import print_function
import os
import sys
import argparse
import glob
import itertools
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import atomset
try:
    basestring
except NameError:
    basestring = str


def parseArguments():
    """
        Parse the command-line options

        :returns: :py:class:`.Clustering`, int, int, int, str -- Clustering
            object, number of trajectory, number of snapshot, number of epoch,
            output path where to write the files
    """
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("epoch", type=str, help="Path to the epoch to search the snapshot")
    parser.add_argument("trajectory", type=int, help="Trajectory number")
    parser.add_argument("snapshot", type=int, help="Snapshot to select (in accepted steps)")
    parser.add_argument("-o", type=str, default=None, help="Output path where to write the files")
    parser.add_argument("--name", type=str, default="pathway.pdb", help="Name of the pdb to write the files")
    parser.add_argument("--top", type=str, default=None, help="Name of the pdb topology for loading non-pdb trajectories")
    args = parser.parse_args()
    return args.trajectory, args.snapshot, args.epoch, args.o, args.name, args.top


def main(trajectory, snapshot, epoch, outputPath, out_filename, topology):
    if outputPath is not None:
        outputPath = os.path.join(outputPath, "")
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    else:
        outputPath = ""
    if topology is not None:
        topology_contents = utilities.getTopologyFile(topology)
    else:
        topology_contents = None
    if os.path.exists(outputPath+out_filename):
        # If the specified name exists, append a number to distinguish the files
        name, ext = os.path.splitext(out_filename)
        out_filename = "".join([name, "_%d", ext])
        i = 1
        while os.path.exists(outputPath+out_filename % i):
            i += 1
        out_filename %= i
    pathway = []
    # Strip out trailing backslash if present
    pathPrefix, epoch = os.path.split(epoch.rstrip("/"))
    sys.stderr.write("Creating pathway...\n")
    while True:
        filename = glob.glob(os.path.join(pathPrefix, epoch, "*traj*_%d.*" % trajectory))
        snapshots = utilities.getSnapshots(filename[0], topology=topology)
        if not isinstance(snapshots[0], basestring):
            new_snapshots = []
            for i in range(snapshot+1):
                snapshot = snapshots.slice(i, copy=False)
                PDB = atomset.PDB()
                PDB.initialise(snapshot, topology=topology_contents)
                new_snapshots.append(PDB.get_pdb_string())
            snapshots = new_snapshots
        else:
            snapshots = snapshots[:snapshot+1]
        pathway.insert(0, snapshots)
        if epoch == '0':
            # Once we get to epoch 0, we just need to append the trajectory
            # where the cluster was found and we can break out of the loop
            break
        procMapping = open(os.path.join(pathPrefix, epoch, "processorMapping.txt")).read().rstrip().split(':')
        epoch, trajectory, snapshot = map(int, procMapping[trajectory-1][1:-1].split(','))
        epoch = str(epoch)
    sys.stderr.write("Writing pathway...\n")
    with open(outputPath+out_filename, "a") as f:
        f.write("ENDMDL\n".join(itertools.chain.from_iterable(pathway)))


if __name__ == "__main__":
    traj, num_snapshot, num_epoch, output_path, output_filename, top = parseArguments()
    main(traj, num_snapshot, num_epoch, output_path, output_filename, top)
