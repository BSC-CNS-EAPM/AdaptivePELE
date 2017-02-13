from utilities import clusteringUtilities
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description="Write the requested cluster "
                                     "structures from a clustering object")
    parser.add_argument('clObject', type=str)
    parser.add_argument('outputPath', type=str,
                        help="Path where to write the structures, including "
                        "name of the files, i.e output/path/cluster.pdb")
    parser.add_argument("structures", nargs='*', type=list, default=None,
                        help="Structures to write")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parseArgs()
    clusteringUtilities.writeStructures(args.clObject, args.structures,
                                        args.outputPath)
