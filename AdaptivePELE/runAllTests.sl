#!/bin/bash
#SBATCH --job-name="AdaptiveTest"
#SBATCH -D .
#SBATCH --output=AdaptiveTest.out
#SBATCH --error=AdaptiveTest.err
#SBATCH --ntasks=4
#SBATCH --time=00:40:00
#SBATCH --qos=debug


python runAllTests.py
