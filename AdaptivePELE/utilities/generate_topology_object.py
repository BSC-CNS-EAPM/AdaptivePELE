import os
import glob
import argparse
from AdaptivePELE.utilities import utilities


def parseArguments():
    """
        Parse the command-line options

        :returns: object -- Object containing the options passed
    """
    desc = "Program that writes a pkl for the topology object."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("path", type=str, help="Path where the topology files are stored")
    args = parser.parse_args()

    return args.path


def main(top_path):
    sim_folder = os.path.abspath(os.path.join(top_path, os.path.pardir))
    epochs = utilities.get_epoch_folders(sim_folder)
    top = utilities.Topology(top_path)
    topology_files = glob.glob(os.path.join(top_path, "topology*.pdb"))
    topology_files.sort(key=utilities.getTrajNum)
    top.setTopologies(topology_files)
    for epoch in epochs:
        top.readMappingFromDisk(os.path.join(sim_folder, epoch), int(epoch))
    top.writeTopologyObject()


if __name__ == "__main__":
    path = parseArguments()
    main(path)
