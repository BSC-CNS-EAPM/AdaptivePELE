#!/bin/bash
#SBATCH --job-name="AdaptiveTest"
#SBATCH -D .
#SBATCH --output=AdaptiveTest.out
#SBATCH --error=AdaptiveTest.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=01:25:00
#SBATCH --qos=debug


python runAllTests.py --exclude MD_CUDA
