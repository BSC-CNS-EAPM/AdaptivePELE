import pickle
import numpy as np
import matplotlib.pyplot as plt
import atomset
from mpl_toolkits.mplot3d import Axes3D
import pdb as debug

def extractCOMMatrix(clusters):
    n = len(clusters)
    cluster_matrix = np.zeros((n,3))
    metrics = np.zeros(n)
    population = np.zeros(n)
    total_elements = 0
    for index,cluster in enumerate(clusters):
        metrics[index] = cluster.metric
        ligandPDB = atomset.PDB()
        ligandPDB.initialise(cluster.pdb.pdb, resname="AEN")
        cluster_matrix[index,:] = ligandPDB.extractCOM()
        population[index] = cluster.elements
        total_elements += cluster.elements
    return cluster_matrix, metrics, total_elements, population

def plotClusters2D(cluster_matrix, metrics,title):
    ccx = cluster_matrix[:,0]
    ccy = cluster_matrix[:,1]
    ccz = cluster_matrix[:,2]
    fig,axes = plt.subplots(nrows=2, ncols=2, sharex='col',
                            sharey='row')
    fig.suptitle(title)
    scatter1 = axes[0][0].scatter(ccx,ccy,c=metrics) #,label="Set %d"%index)
    scatter2 = axes[0][1].scatter(ccz,ccy,c=metrics) #,label="Set %d"%index)
    scatter3 = axes[1][0].scatter(ccx,ccz,c=metrics) #,label="Set %d"%index)
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

def plotClusters(cluster_matrix, metrics,title):
    fig = plt.figure()
    ax = Axes3D(fig)
    ccx = cluster_matrix[:,0]
    ccy = cluster_matrix[:,1]
    ccz = cluster_matrix[:,2]
    print title
    print ccx.size
    fig.suptitle(title)
    scatter1 = ax.scatter(ccx,ccy,zs = ccz,c=metrics) #,label="Set %d"%index)
    fig.colorbar(scatter1)
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    return fig

#with open("ClDouble.pkl","r") as f:
#    ClDouble = pickle.load(f)
#matrixDouble, metricsDouble,totalElementsDouble, popDouble = extractCOMMatrix(ClDouble.clusters)
#plotDouble = plotClusters(matrixDouble, metricsDouble, 'Clusters ContactMap Double')
#plotDouble.savefig('results/contactmapDouble.png')
#plotDoublepop = plotClusters(matrixDouble, popDouble, 'Clusters ContactMap Double')
#plotDouble.savefig('results/contactmapDoublepop.png')
#print "Number of elements", totalElementsDouble

with open("ClAgg.pkl","r") as f:
    ClAgg = pickle.load(f)
contactsAgg = []
for cluster in ClAgg.clusters:
        contactsAgg.append(cluster.contacts)
fig = plt.figure()
fig.suptitle("Contacts of Agglomerative clustering")
plt.hist(np.array(contactsAgg))

# matrixAgg, metricsAgg,totalElementsAgg, popAgg = extractCOMMatrix(ClAgg.clusters)
# plotAgg = plotClusters(matrixAgg, metricsAgg, 'Clusters ContactMap Agg')
# plotAgg.savefig('results/contactmapAgg.png')
# plotAggpop = plotClusters(matrixAgg, popAgg, 'Clusters ContactMap Agg')
# plotAgg.savefig('results/contactmapAggpop.png')
# print "Number of elements", totalElementsAgg


with open("ClCont.pkl","r") as f:
    ClCont = pickle.load(f)
contactsCont = []
for cluster in ClCont.clusters:
        contactsCont.append(cluster.contacts)
fig = plt.figure()
fig.suptitle("Contacts of Contacts clustering")
plt.hist(np.array(contactsCont))
# matrixCont, metricsCont, totalElementsCont, popCont = extractCOMMatrix(ClCont.clusters)
# plotCont = plotClusters(matrixCont, metricsCont, 'Clusters Contacts')
# plotCont.savefig('results/contactClusters.png')
# plotContpop = plotClusters(matrixCont, popCont, 'Clusters Contacts')
# plotContpop.savefig('results/contactClusterspop.png')
# print "Number of elements", totalElementsCont
# debug.set_trace()
plt.show()
