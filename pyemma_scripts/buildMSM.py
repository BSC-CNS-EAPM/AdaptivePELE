import MSMblocks
import numpy as np
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description="Build a MSM using PyEMMA "
                                     "package")
    parser.add_argument('controlFile', type=str)
    args = parser.parse_args()
    return args


def readParams(control_file):
    params = MSMblocks.readParams(control_file)
    trajectoryFolder = params["trajectoryFolder"]
    trajectoryBasename = params["trajectoryBasename"]
    numClusters = params["numClusters"]
    lagtimes = params["lagtimes"]
    numPCCA = params["numPCCA"]
    itsOutput = params["itsOutput"]
    numberOfITS = params["numberOfITS"]
    itsErrors = params["itsErrors"]
    error_estimationCK = params["error_estimationCK"]
    state_labels = params["state_labels"]
    if state_labels is None:
        state_labels = 'auto'  # json returns string as
    # unicode, and this breaks some code in pyemma
    outfile_fluxTPT = params["outfile_fluxTPT"]
    return trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, itsErrors, error_estimationCK, state_labels, outfile_fluxTPT


def main(control_file):

    # parameters
    trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, itsErrors, error_estimationCK, state_labels, outfile_fluxTPT = readParams(control_file)

    # program
    prepareMSM = MSMblocks.PrepareMSM(numClusters, trajectoryFolder, trajectoryBasename)
    cl = prepareMSM.getClusteringObject()
    calculateMSM = MSMblocks.MSM(cl, lagtimes, numPCCA, itsOutput, numberOfITS,
                                 itsErrors, error_estimationCK)
    calculateMSM.estimate()
    MSM_object = calculateMSM.getMSM_object()
    TPTinstance = MSMblocks.TPT(MSM_object, cl, outfile_fluxTPT, state_labels)
    TPT_Object = TPTinstance.getTPTObject()
    coarseTPT_Object = TPTinstance.getCoarseTPTObject()

    # Free energy estimation
    print "Calculating free energies..."
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
    args = parseArgs()
    main(args.controlFile)
