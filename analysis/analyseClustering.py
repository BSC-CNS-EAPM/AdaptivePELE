import pickle
import numpy as np
import matplotlib.pyplot as plt
import atomset
from mpl_toolkits.mplot3d import Axes3D
import pdb as debug


def extractCOMMatrix(clusters, resname):
    n = len(clusters)
    cluster_matrix = np.zeros((n, 3))
    metrics = np.zeros(n)
    population = np.zeros(n)
    total_elements = 0
    contacts = np.zeros(n)
    for index, cluster in enumerate(clusters):
        metrics[index] = cluster.metric
        contacts[index] = cluster.contacts
        ligandPDB = atomset.PDB()
        ligandPDB.initialise(cluster.pdb.pdb, resname=resname)
        cluster_matrix[index, :] = ligandPDB.extractCOM()
        population[index] = cluster.elements
        total_elements += cluster.elements
    return cluster_matrix, metrics, total_elements, population, contacts


def plotClusters2D(cluster_matrix, metrics, title):
    ccx = cluster_matrix[:, 0]
    ccy = cluster_matrix[:, 1]
    ccz = cluster_matrix[:, 2]
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex='col',
                             sharey='row')
    fig.suptitle(title)
    scatter1 = axes[0][0].scatter(ccx, ccy, c=metrics)  # ,label="Set %d"%index)
    axes[0][1].scatter(ccz, ccy, c=metrics)  # ,label="Set %d"%index)
    axes[1][0].scatter(ccx, ccz, c=metrics)  # ,label="Set %d"%index)
    fig.colorbar(scatter1)

#    axes[1][0].legend(loc='center right', bbox_to_anchor=[1.8,0.5])
    axes[0][0].set_ylabel('y')
    axes[1][0].set_ylabel('z')
    axes[1][0].set_xlabel('x')
    axes[1][1].axis('off')
    axes[0][1].set_xticks(axes[1][1].get_xticks())
    axes[0][1].set_xticklabels(axes[1][1].get_xticklabels())
    axes[0][1].set_xlabel('z')
    return fig


def plotClusters(cluster_matrix, metrics, title):
    fig = plt.figure()
    ax = Axes3D(fig)
    ccx = cluster_matrix[:, 0]
    ccy = cluster_matrix[:, 1]
    ccz = cluster_matrix[:, 2]
    print title
    print ccx.size
    fig.suptitle(title)
    scatter1 = ax.scatter(ccx, ccy, zs=ccz, c=metrics)
    fig.colorbar(scatter1)
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    return fig

def plotClusteringData(pklObjectFilename, resname, metricPlotFilename="", populationPlotFilename="", contactsPlotFilename=""):
    with open(pklObjectFilename, "r") as f:
        clObject = pickle.load(f)

    comCoord, metrics, totalElements, population, contacts = extractCOMMatrix(clObject.clusters.clusters, resname)

    plot = plotClusters(comCoord, metrics, 'Clusters Contacts')
    if metricPlotFilename: plot.savefig(metricPlotFilename)

    plotContpop = plotClusters(comCoord, population, 'Clusters Contacts')
    if populationPlotFilename: plotContpop.savefig(populationPlotFilename)

    plotContpop = plotClusters(comCoord, contacts, 'Clusters Contacts')
    if contactsPlotFilename: plotContpop.savefig(contactsPlotFilename)

    print "Number of elements", totalElements

if __name__ == "__main__":
    resname = "ALJ"

    # with open("ClDouble.pkl","r") as f:
    #    ClDouble = pickle.load(f)
    # matrixDouble, metricsDouble,totalElementsDouble, popDouble = extractCOMMatrix(ClDouble.clusters.clusters,resname)
    # plotDouble = plotClusters(matrixDouble, metricsDouble, 'Clusters ContactMap Double')
    # plotDouble.savefig('results/contactmapDouble.png')
    # plotDoublepop = plotClusters(matrixDouble, popDouble, 'Clusters ContactMap Double')
    # plotDouble.savefig('results/contactmapDoublepop.png')
    # print "Number of elements", totalElementsDouble

    # with open("ClAgg.pkl", "r") as f:
    #     ClAgg = pickle.load(f)
    # contactsAgg = []
    # for cluster in ClAgg.clusters.clusters:
    #         contactsAgg.append(cluster.contacts)
    # fig = plt.figure()
    # fig.suptitle("Contacts of Agglomerative clustering")
    # plt.hist(np.array(contactsAgg))
    # fig.savefig('results/aggContacts.png')

    # matrixAgg, metricsAgg,totalElementsAgg, popAgg = extractCOMMatrix(ClAgg.clusters.clusters, resname)
    # plotAgg = plotClusters(matrixAgg, metricsAgg, 'Clusters ContactMap Agg')
    # plotAgg.savefig('results/contactmapAgg.png')
    # plotAggpop = plotClusters(matrixAgg, popAgg, 'Clusters ContactMap Agg')
    # plotAgg.savefig('results/contactmapAggpop.png')
    # print "Number of elements", totalElementsAgg

    pklObjectFilename = "17/clustering/object.pkl"
    resname = "K5Y"
    metricPlotFilename = ""#"results/contactClusters.png"
    populationPlotFilename = ""#"results/contactClusterspop.png"
    contactsPlotFilename = ""#"results/contactClustersContacts.png"

    plotClusteringData(pklObjectFilename, resname, metricPlotFilename, populationPlotFilename, contactsPlotFilename)
    plt.show()
    sys.exit()


    with open("ClCont.pkl", "r") as f:
        ClCont = pickle.load(f)
    contactsCont = []
    for cluster in ClCont.clusters.clusters:
            contactsCont.append(cluster.contacts)
    fig = plt.figure()
    fig.suptitle("Contacts of Contacts clustering")
    plt.hist(np.array(contactsCont))
    fig.savefig('results/contContacts.png')
    matrixCont, metricsCont, totalElementsCont, popCont = extractCOMMatrix(ClCont.clusters.clusters,resname)
    plotCont = plotClusters(matrixCont, metricsCont, 'Clusters Contacts')
    plotCont.savefig('results/contactClusters.png')
    plotContpop = plotClusters(matrixCont, popCont, 'Clusters Contacts')
    plotContpop.savefig('results/contactClusterspop.png')
    print "Number of elements", totalElementsCont

    with open("ClAcc.pkl", "r") as f:
        ClAcc = pickle.load(f)
    contactsAcc = []
    for cluster in ClAcc.clusters.clusters:
            contactsAcc.append(cluster.contacts)
    fig = plt.figure()
    fig.suptitle("Contacts of Accumulative clustering")
    plt.hist(np.array(contactsAcc))
    fig.savefig('results/accContacts.png')
    matrixAcc, metricsAcc, totalElementsAcc, popAcc = extractCOMMatrix(ClAcc.clusters.clusters, resname)
    plotAcc = plotClusters(matrixAcc, metricsAcc, 'Clusters Accumulative')
    plotAcc.savefig('results/accClusters.png')
    plotAccpop = plotClusters(matrixAcc, popAcc, 'Clusters Accumulative')
    plotAccpop.savefig('results/accClusterspop.png')
    print "Number of elements", totalElementsAcc

    # debug.set_trace()
    plt.show()
