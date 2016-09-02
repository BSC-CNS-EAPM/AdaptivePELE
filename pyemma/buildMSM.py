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
    #TODO: Define blocks with tasks to make the program more modular
    #TODO: Add lines saying what the program is doing

    ### parameters
    trajectoryFolder = "test/MSM2" 
    trajectoryBasename = "*traj_*" 

    numClusters = 100

    lagtimes = [1, 2, 5, 10, 20, 50, 100, 200, 500]
    itsOutput = "its.png"
    numberOfITS = -1
    itsErrors=None #'bayes'
    nsetsCK = 2
    error_estimationCK=False
    membershipsCK=None
    numPCCA = 7

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

    plot_its = None #just a placeholder so it doesn't go out of scope
    is_converged = False
    #its
    while not is_converged:
        
        its_object = msm.calculateITS(cl.dtrajs, lagtimes, itsErrors)
        plot_its = msm.plotITS(its_object, itsOutput, numberOfITS)
        plt.show()
        while True:
            convergence_answer = raw_input("Has the ITS plot converged?[y/n] ")
            convergence_answer.rstrip()
            convergence_answer = convergence_answer or "y" #Making yes the default
            #answer
            if convergence_answer.lower() == "y" or convergence_answer.lower() == "yes":
                is_converged = True
                lagtime_str = raw_input("Please input the lagtime to construct the MSM: ")
                lagtime = int(lagtime_str.rstrip())
                break
            elif convergence_answer.lower() == "n" or convergence_answer.lower() == "no":
                break
            else:
                print "Answer not valid. Please answer yes or no"
        if not is_converged:
            new_lagtimes = raw_input("Do you want to define new lagtimes or add to the previous?[add(a)/new(n)] ")
            new_lagtimes.rstrip()
            if new_lagtimes.lower() == "add" or new_lagtimes.lower() == "a":
                lag_list = raw_input("Please input the lagtimes you want to add separated by a space: ")
                lag_list.rstrip()
                lagtimes.extend(map(int,lag_list.split(" ")))
            elif new_lagtimes.lower() == "new" or new_lagtimes.lower() == "n":
                lag_list = raw_input("Please input the new lagtimes separated by a space: ")
                lag_list.rstrip()
                lagtimes = map(int,lag_list.split(" "))
            lagtimes.sort()
    #estimation
    MSM_object = msm.estimateMSM(cl.dtrajs, lagtime)

    #connectivity
    if msm.is_connected(MSM_object):
        print "The MSM estimated is fully connected"
    else:
        print "The MSM estimated is not fully connected"
        unconnected_sets = get_connected_sets(MSM_object)
        print "The MSM estimated has %d connected sets with sizes:" % len(unconnected_sets)
        for index, uncon_set in enumerate(unconnected_sets):
            print "Set %d has %d elements" % (index, uncon_set.size)

    #PCCA
    MSM_object = msm.calculatePCCA(MSM_object, numPCCA)

    #Chapman-Kolgomorov validation
    #TODO: add decision-making loop for validity of CK test
    nsetsCK = len(MSM_object.metastable_sets)
    membershipsCK = MSM_object.metastable_memberships
    CKObject = msm.ChapmanKolmogorovTest(MSM_object,
                                         nsetsCK,memberships=membershipsCK,
                                         error_estimation=error_estimationCK)
    msm.plotChapmanKolmogorovTest(CKObject)
    plt.show()

    #TODO: TPT block
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

    Gpmf = np.expand_dims(Gpmf, axis=1)
    output = np.hstack((cl.clustercenters, Gpmf))
    np.savetxt("output.txt", output)


if __name__ == "__main__":
    main()
