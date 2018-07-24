#!/bin/bash
#BSUB -J "AdaptiveTest"
#BSUB -o AdaptiveTest.out
#BSUB -e AdaptiveTest.err
#BSUB -n 5
#BSUB -W 02:00
#BSUB -q bsc_debug


python runAllTests.py
