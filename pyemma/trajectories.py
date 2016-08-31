import pyemma.coordinates as coor
import glob

def loadCoordinates(path, trajectories_basename='*traj*.pdb'):
    """ Load the coordinates from the simulation into an object that can be
    used with the rest of PyEMMA tools.
    Returns an array with the trajectories"""

    path_to_file = os.path.join(path,trajectories_basename)
    files = glob.glob(path_to_file)
    return coor.load(files)

#More preprocessing staff?

def clusterTrajectories(trajectories, numClusters):
    """ Cluster the trajectories into numClusters clusters using kmeans
    algorithm.
    Returns a KmeansClusteringObject
    """
    return coor.cluster_kmeans(data=trajectories, k=numClusters)
