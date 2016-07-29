import numpy as np
import math

spawning = "INVERSELY_PROPORTIONAL"
numberOfProcessors = 64

sizes = np.array([64,59,22,19,17,17,11,11,9,5,5,5,4,2,2,2,1,1,1,1,1,1,1,1,1,1,1])

weights = 1./sizes
weights /= sum(weights)

degeneracy = []
if spawning == "INVERSELY_PROPORTIONAL":
    trajToDistribute = numberOfProcessors-1
    for i, weight in enumerate(weights):
        degeneracy.append(int(weight*trajToDistribute)) 

elif spawning == "AT_LEAST_ONE":
    numberOfInitialStructures = len(weights)
    trajToDistribute = numberOfProcessors-1-numberOfInitialStructures
    for i, weight in enumerate(weights):
        degeneracy.append(1 + int(weight*trajToDistribute)) #at least one per cluster


decimalPart = []
for i, weight in enumerate(weights):
    decimalPart.append(math.modf(weight*trajToDistribute)[0])
sortedDecimals = np.argsort(decimalPart)
sortedDecimals = sortedDecimals[::-1]

leftProcessors = numberOfProcessors-1-sum(degeneracy)
for i in range(leftProcessors):
    degeneracy[sortedDecimals[i]] += 1

print len(sizes)
print sizes
print degeneracy
