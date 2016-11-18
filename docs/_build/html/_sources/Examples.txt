Adaptive Control File
=====================

The control file needed to run AdaptiveSampling is a json document that
contains 4 sections: general parameters, simulation parameters, clustering
parameters and spawning parameters.

general Params
--------------

The general parameters block has five mandatory fields:

* restart: boolean (true or false) This parameter specifies wether you want to
  continue a previous simulation or not

* debug: boolean Run adaptive in debug mode, without calling
  any propagation algorithm

* outputPath: string The path where the results of the simulation will be
  written

* initialStructures: list The path(s) to the intial structure(s)  

* writeAllClusteringStructures: boolean Wether to write all the structures of
  the clusters 

Additionaly, it can also have a nativeStructure parameter, a string containing
the path to the native structure. This structure will only be used to correct
the RMSD in case of symmetries. The symmetries will also have to be specified
(see Clustering section)


