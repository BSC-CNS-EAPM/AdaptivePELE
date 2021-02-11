#!/bin/bash
#SBATCH --job-name="AdaptiveTest"
#SBATCH -D .
#SBATCH --output=AdaptiveTest.out
#SBATCH --error=AdaptiveTest.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=01:45:00
#SBATCH --qos=debug


module purge #2> /dev/null
module load intel mkl impi python/2.7.13 boost/1.64.0_py2 gcc openmp 2> /dev/null
export PYTHONPATH="/home/bsc72/bsc72021/adaptiveSampling/AdaptiveSampling"
export PYTHONPATH="/gpfs/projects/bsc72/MSM_PELE/bin:/home/bsc72/bsc72021/python-tools:$PYTHONPATH"
export PYTHONPATH="$PYTHONPATH:/gpfs/projects/bsc72/lib/site-packages"
python runAllTests.py --exclude MD_CUDA MD R

module purge 2> /dev/null
module load intel mkl impi  2> /dev/null #boost/1.64.0_py2 gcc 2> /dev/null
module load ANACONDA/5.0.1 2> /dev/null
module load ambertools/17 2> /dev/null
module load openmp 
export PYTHONPATH="/home/bsc72/bsc72021/adaptiveSampling/AdaptiveSampling"
export PYTHONPATH="/gpfs/projects/bsc72/MSM_PELE/bin:/home/bsc72/bsc72021/python-tools:$PYTHONPATH"
export PYTHONPATH="$PYTHONPATH:/gpfs/projects/bsc72/lib/python3"

python runAllTests.py --run MD R
