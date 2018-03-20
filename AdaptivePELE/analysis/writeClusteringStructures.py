"""
    Write specified cluster representative structures to pdb
"""
from AdaptivePELE.utilities import clusteringUtilities
import argparse


def parseArgs():
    """
        Parse command line arguments

        :returns: object -- Object containing command line options
    """
    parser = argparse.ArgumentParser(description="Write the requested cluster "
                                     "structures from a clustering object")
    parser.add_argument('clObject', type=str, help="Path to the clustering object")
    parser.add_argument('outputPath', type=str,
                        help="Path where to write the structures, including "
                        "name of the files, i.e output/path/cluster.pdb")
    parser.add_argument("structures", nargs='*', type=list, default=None,
                        help="Structures to write")
    parser.add_argument("--threshold", type=float, default=None,
                        help="Only print those structures with mathcing threshold")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    arguments = parseArgs()
    if arguments.threshold is not None:
        condition = lambda x: abs(x.threshold-arguments.threshold) < 0.01
    else:
        condition = None
    clusteringUtilities.writeStructures(arguments.clObject, arguments.structures, condition,
                                        arguments.outputPath)
