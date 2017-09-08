import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from AdaptivePELE.utilities import utilities



def parseArgs():
    desc = "Write the conformation network to a PDB file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("clusteringObject", type=str, help="Path to the clustering object")
    parser.add_argument("metricCol", type=int, help="Column of the metric of interest")
    parser.add_argument("filename", type=str, help="Name of the pdb file")
    args = parser.parse_args()
    return args.clusteringObject, args.metricCol, args.filename


def writePDB(X, Y, Z, metric, f):
    templateLine = "HETATM%s  H%sCLT L 502    %s%s%s  0.75%s           H\n"

    content = ""
    for i, xn in enumerate(X):
        number = str(i).rjust(5)
        number3 = str(i).ljust(3)
        yn = Y[i]
        zn = Z[i]
        x = ("%.3f"%xn).rjust(8)
        y = ("%.3f"%yn).rjust(8)
        z = ("%.3f"%zn).rjust(8)
        g = ("%.3f"%metric[i]).rjust(8)

        content += templateLine%(number, number3, x, y, z, g)
    f.write(content)
    return f


def writeConnectionsPDB(network, f):
    template = "CONECT%s%s%s%s%s%s%s%s%s%s              \n"
    for node in network.nodes_iter():
        contents = [node]
        for source, target in network.out_edges_iter(node):
            contents.append(target)
            if len(contents) == 10:
                f.write(template % tuple(map(lambda x: str(x).rjust(5), contents)))
                contents = [node]
        if contents:
            contents += [""]*(10-len(contents))
            f.write(template % tuple(map(lambda x: str(x).rjust(5), contents)))
    return f


def readClustering(clusteringPath, metricCol):
    print "Reading clustering object..."
    clustering = utilities.readClusteringObject(clusteringPath)
    network = clustering.conformationNetwork.network
    metrics = [cl.metrics[metricCol] for cl in clustering.clusterIterator()]
    return clustering, network, metrics


def getCoords(clustering):
    print "Extracting nodes coordinates..."
    coords = [cl.pdb.getCOM() for cl in clustering.clusterIterator()]
    # Attempt to delete clustering object from memory, since it is not
    # needed and is very heavy
    del clustering
    coords = np.array(coords)
    Xn = coords[:, 0]
    Yn = coords[:, 1]
    Zn = coords[:, 2]
    return Xn, Yn, Zn

if __name__ == "__main__":
    clusteringPath, metricCol, filename = parseArgs()

    clustering, network, metrics = readClustering(clusteringPath, metricCol)
    Xn, Yn, Zn = getCoords(clustering)

    # Write info in pdb file
    f = open(filename, "w")
    f = writePDB(Xn, Yn, Zn, metrics, f)
    f = writeConnectionsPDB(network, f)
    f.close()

# Plot using matplotlib

# Xe = []
# Ye = []
# Ze = []
# for e in network.edges_iter():
#     Xe += [[coords[e[0], 0], coords[e[1], 0]]]
#     Ye += [[coords[e[0], 1], coords[e[1], 1]]]
#     Ze += [[coords[e[0], 2], coords[e[1], 2]]]
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # ax.plot(Xe, Ye, Xe)
# for i, xl in enumerate(Xe):
#     ax.plot(xl, Ye[i], Ze[i], 'k', alpha=0.1, linewidth=0.5)
# ax.scatter(Xn, Yn, Zn, c=metrics)
# # plt.colorbar()
# plt.show()
