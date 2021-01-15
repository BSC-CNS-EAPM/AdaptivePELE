environments = {
        "mn4": ["module load intel mkl impi ANACONDA/5.0.1 openmp 2> /dev/null", "module load intel mkl impi python/2.7.13 boost/1.64.0_py2 gcc openmp 2> /dev/null"],
        "nord3": ["module load gcc/6.1.0 MKL GTK+3/3.2.4 intel/2017.0.098 python/2.7.12 2> /dev/null",
                  'module load gcc/8.4.0 impi/2017.4 MKL/2017.4 OPENSSL/1.1.1h python/3.6.9-pip boost/1.74.0 rdkit/2020.09.1 2> /dev/null; python() { python3 "$@"; }; export -f python'],
        "nvidia": ["module load intel/16.0.2 amber/16 python/2.7.2 2> /dev/null"],
        "power": ["module load gcc/6.4.0 python/3.6.5 openmpi/3.0.0 cuda/9.1 ambertools/18 2> /dev/null"]
        }
