#!/bin/bash
#SBATCH --job-name="AdaptiveTest_CUDA"
#SBATCH -D .
#SBATCH --output=AdaptiveTest_CUDA.out
#SBATCH --error=AdaptiveTest_CUDA.err
#SBATCH --ntasks=2
##SBATCH --cpus-per-task=5
#SBATCH --cpus-per-task=80
#SBATCH --time=01:20:00
#SBATCH --qos=debug
#SBATCH --gres gpu:2


module purge 2> /dev/null
module load gcc/6.4.0 2> /dev/null   
module load openmpi/3.0.0 2> /dev/null
module load cuda/9.1 2> /dev/null
module load python/3.6.5 2> /dev/null
module load ambertools/18 2> /dev/null
export PYTHONPATH="/home/bsc72/bsc72021/AdaptiveCTEP/adaptivePELE/:/gpfs/projects/bsc72/lib/site-packages-cte"
python runAllTests.py --exclude MD_CUDA Ad
srun python runAllTests.py --run MD_CUDA
