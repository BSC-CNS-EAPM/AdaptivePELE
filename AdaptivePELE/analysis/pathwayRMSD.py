import networkx as nx
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import RMSDCalculator

def weight(path, conf):
    w = 0
    for i in xrange(1, len(path)):
        w += conf[path[i-1]][path[i]]['metric']
    return w

metricCol = 4
RMSDCalc = RMSDCalculator.RMSDCalculator()
cl = utilities.readClusteringObject("ClCont.pkl")
nodeFin = cl.getOptimalMetric(column=metricCol)
conf = nx.DiGraph()
nx.read_edgelist("conformationNetwork_4DAJ.edgelist", create_using=conf,
                 data=True, nodetype=int)
for source, target in conf.edges_iter():
    clusterSource = cl.getCluster(source)
    clusterTarget = cl.getCluster(target)
    conf[source][target]['metric'] = RMSDCalc.computeRMSD(clusterSource.pdb, clusterTarget.pdb)

# paths = nx.all_simple_paths(conf, 0, nodeFin)
paths = nx.shortest_simple_paths(conf, 0, nodeFin, weight='metric')
for i, path in enumerate(paths, start=1):
    if i > 50:
        break
    elif i < 3:
        cl.writePathwayTrajectory(path, "pathway_%d.pdb" % (i-1))
    print "Path ", i, "length ", weight(path, conf)
