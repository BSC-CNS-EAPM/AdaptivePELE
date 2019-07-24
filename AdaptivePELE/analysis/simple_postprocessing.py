from __future__ import print_function
import os
import glob
import argparse
import mdtraj as md
import multiprocessing as mp
from AdaptivePELE.utilities import utilities


def parse_arguments():
    desc = "Program that performs simple postprocessing of MD simulations."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--dont_image", action="store_false", help="Flag to set whether trajectories should be imaged before the alignment (if not specfied performs the imaging)")
    parser.add_argument("--processors", type=int, default=4, help="Number of cpus to use")
    parser.add_argument("--trajectory_name", default="trajectory", help="Name of the trajectory files")
    parser.add_argument("--path", default=".", help="Path to the simulation files")
    args = parser.parse_args()
    return args.trajectory_name, args.path, args.processors, args.dont_image


def process_traj(inputs):
    top_ind, traj_name, epoch, traj_num, imaging = inputs
    ext = os.path.splitext(traj_name)[1]
    utilities.print_unbuffered("Processing trajectory", traj_name)
    top = md.load("topologies/topology_%s.pdb" % top_ind)
    atoms = top.top.select("backbone")
    t = md.load(traj_name, top="topologies/system_%s.prmtop" % top_ind)
    if imaging:
        t.image_molecules(inplace=True)
    t.superpose(top, atom_indices=atoms)
    t.save(os.path.join(epoch, "trajectory_postprocessed_%d%s" % (traj_num, ext)))


def main(trajectory_name, path, n_processors, imaging):
    epochs = utilities.get_epoch_folders(path)
    to_process = []
    pool = mp.Pool(n_processors)
    trajectory_glob = trajectory_name + "_*"
    for epoch in epochs:
        with open(os.path.join(epoch, "topologyMapping.txt")) as f:
            top_map = f.read().rstrip().split(":")
        for traj in glob.glob(os.path.join(path, epoch, trajectory_glob)):
            traj_num = utilities.getTrajNum(traj)
            to_process.append((top_map[traj_num-1], traj, epoch, traj_num, imaging))

    pool.map(process_traj, to_process)
    pool.close()
    pool.terminate()

if __name__ == "__main__":
    traj_names, path_files, proc_num, image = parse_arguments()
    main(traj_names, path_files, proc_num, image)
