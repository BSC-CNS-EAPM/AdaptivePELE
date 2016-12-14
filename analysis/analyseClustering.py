import pickle
import numpy as np
import matplotlib.pyplot as plt
from atomset import atomset
from mpl_toolkits.mplot3d import Axes3D
import pdb as debug


def extractCOMMatrix(clusters, resname):
    """ Extract a matrix contaning the coordinates of the center of mass of
    the ligand for each cluster structure

        clusters [In] List of clusters
        resname [In] Residue name of the ligand in the pdb
    """
    n = len(clusters)
    cluster_matrix = np.zeros((n, 3))
    metrics = np.zeros(n)
    population = np.zeros(n)
    total_elements = 0
    contacts = np.zeros(n)
    for index, cluster in enumerate(clusters):
        metrics[index] = cluster.metrics[cluster.metricCol]
        contacts[index] = cluster.contacts
        ligandPDB = atomset.PDB()
        ligandPDB.initialise(cluster.pdb.pdb, resname=resname)
        cluster_matrix[index, :] = ligandPDB.extractCOM()
        population[index] = cluster.elements
        total_elements += cluster.elements
    return cluster_matrix, metrics, total_elements, population, contacts


def plotClusters2D(cluster_matrix, metrics, title):
    """ Create all combination of xyz projections in 2D of the scatter plot
    of the center of mass of the ligand with a colormap given by a certain
    quantity (usually a metric or the clusters population)

        cluster_matrix [In] matrix contaning the coordinates of the center of
        mass of the ligand for each cluster structure
        metrics [In] Array with the quantity that will be used to create the
        colormap
        title [In] Title for the plot figure
    """
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
    """ Create a 3D scatter plot of the center of mass of the ligand with a
    colormap given by a certain quantity
    (usually a metric or the clusters population)

        cluster_matrix [In] matrix contaning the coordinates of the center of
        mass of the ligand for each cluster structure
        metrics [In] Array with the quantity that will be used to create the
        colormap
        title [In] Title for the plot figure
    """
    fig = plt.figure()
    ax = Axes3D(fig)
    ccx = cluster_matrix[:, 0]
    ccy = cluster_matrix[:, 1]
    ccz = cluster_matrix[:, 2]
    fig.suptitle(title)
    scatter1 = ax.scatter(ccx, ccy, zs=ccz, c=metrics)
    fig.colorbar(scatter1)
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    return fig


def plotClusteringData(pklObjectFilename, resname, titlemetric, titlepopulation,
                       titlecontacts, metricPlotFilename="",
                       populationPlotFilename="", contactsPlotFilename=""):

    with open(pklObjectFilename, "r") as f:
        clObject = pickle.load(f)

    comCoord, metrics, totalElements, population, contacts = extractCOMMatrix(clObject.clusters.clusters, resname)

    plot = plotClusters(comCoord, metrics, titlemetric)
    if metricPlotFilename:
        plot.savefig(metricPlotFilename)

    plotContpop = plotClusters(comCoord, population, titlepopulation)
    if populationPlotFilename:
        plotContpop.savefig(populationPlotFilename)

    plotContcont = plotClusters(comCoord, contacts, titlecontacts)
    if contactsPlotFilename:
        plotContcont.savefig(contactsPlotFilename)

    print "Number of elements", totalElements

if __name__ == "__main__":
    # resname = "ALJ"
    resname = "STR"

    # # Cont
    # pklObjectFilename = "ClCont.pkl"
    # metricPlotFilename = ""  # "results/contactClusters.png"
    # populationPlotFilename = ""  # "results/contactClusterspop.png"
    # contactsPlotFilename = ""  # "results/contactClustersContacts.png"
    # titlemetric = "Metrics Contacts"
    # titlepopulation = "Population Contacts"
    # titlecontacts = "Number of contacts Contacts"

    # plotClusteringData(pklObjectFilename, resname, titlemetric, titlepopulation,
    #                    titlecontacts, metricPlotFilename,
    #                    populationPlotFilename, contactsPlotFilename)

    # Acc
    pklObjectFilename = "ClAcc_PR_heav.pkl"
    metricPlotFilename = "results/metricplotAcc_acc_PR_heav.png"
    populationPlotFilename = "results/populationAcc_acc_PR_heav.png"
    contactsPlotFilename = "results/contactsplotAcc_acc_PR_heav.png"
    titlemetric = "Metrics Accumulative"
    titlepopulation = "Population Accumulative"
    titlecontacts = "Number of contacts Accumulative"

    plotClusteringData(pklObjectFilename, resname, titlemetric, titlepopulation,
                       titlecontacts, metricPlotFilename,
                       populationPlotFilename, contactsPlotFilename)

    plt.show()
