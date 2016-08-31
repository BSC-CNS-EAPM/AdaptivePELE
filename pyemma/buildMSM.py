import trajectories
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyemma.plots as plots
import msm
import os
import numpy as np


def writeClusterCenters(cl, outputFilename):
    np.savetxt(outputFilename, cl.clustercenters)

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

def main():
    ### parameters
    trajectoryFolder = "test/MSM2" 
    trajectoryBasename = "*traj_*" 

    numClusters = 100

    lagtimes = [1, 2, 5, 10, 20, 50, 100, 200, 500]
    itsOutput = "its.png"
    numberOfITS = -1
    itsErrors=None #'bayes'


    ### constants
    discretizedFolder = "discretized"
    clusterCentersFile = os.path.join(discretizedFolder, "clusterCenters.dat")
    discTraj = os.path.join(discretizedFolder, "%s.disctraj")



    #program

    x = trajectories.loadCOMFiles(trajectoryFolder, trajectoryBasename)

    #cluster & assign
    cl = trajectories.clusterTrajectories(x, numClusters)

    #write output
    makeFolder(discretizedFolder)
    writeClusterCenters(cl, clusterCentersFile)

    """
    fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    ax.scatter(clX, clY, clZ)
    plt.show()
    """

    #connectivity

    #its
    its_object = msm.calculateITS(cl.dtrajs, lagtimes, itsErrors)
    msm.plotITS(its_object, itsOutput, numberOfITS)

    #estimation
    MSM_object = msm.estimateMSM(cl.dtrajs, 100)

    kbt = 0.0019872041*300
    pi = MSM_object.stationary_distribution


    """
    for each cluster
        for each disctraj
            these_frames = np.argwhere(disc_traj==i)
            theta_i = cont_traj[these_frames]
        rho_i, _ = np.histogram(theta_i, range=myrange, normed=True)

    pi_total = rho_i * pi_i, donde pi son los counts relativos de cada centro.

    Si el mismo "myrange" para para las distintas bases (regspace y kmeans),
    la Delta_G deberia ser igual!

    _total = \Sigma_i rho_i * pi_i,  donde pi son los counts relativos de cada
    Ojo, que para que esto funcione, cada rho_i tiene que ser una densidad
    verdadera (rho_i.sum() * dx = 1) y pi tambien \Sigma_i pi_i = 0).
    """

    Gpmf = -kbt*np.log(pi/np.max(pi))

    print Gpmf.shape
    Gpmf = np.expand_dims(Gpmf, axis=1)
    print Gpmf.shape
    output = np.hstack((cl.clustercenters, Gpmf))
    np.savetxt("output.txt", output)

    

if __name__ == "__main__":
    main()
