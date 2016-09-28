import blockNames
import densitycalculatortypes
import sys

class DensityCalculatorBuilder():
    def build(self, spawningBlock):
        try:
            densityBlock = spawningBlock[blockNames.SpawningParams.density]
        except KeyError:
            print "Using null density calculator (no preference for any cluster)"
            return NullDensityCalculator()

        try:
            type = densityBlock[blockNames.DensityCalculator.type]
        except KeyError:
            sys.exit("Density calculator must have a type")

        if type == blockNames.DensityCalculator.null or type == blockNames.DensityCalculator.constant:
            print "Using constant density"
            return NullDensityCalculator()
        elif type == blockNames.DensityCalculator.heaviside:
            try:
                paramsBlock = densityBlock[blockNames.DensityCalculator.params]
                values = paramsBlock[blockNames.DensityCalculatorParams.values]
                conditions = paramsBlock[blockNames.DensityCalculatorParams.conditions]
                return DensityCalculatorHeaviside(conditions, values)
            except KeyError:
                print "Using default parameters for Heaviside density calculator"
                return DensityCalculatorHeaviside()
        else:
            sys.exit("Unknown density calculator type! Choices are: " + str(densitycalculatortypes.DENSITY_CALCULATOR_TYPE_TO_STRING_DICTIONARY.values()))

from abc import ABCMeta, abstractmethod
class DensityCalculator():
    def __init__(self):
        self.type = "BaseClass"

    @abstractmethod
    def calculate(self, contacts):
        pass

class DensityCalculatorHeaviside(DensityCalculator):
    #Mostly duplicated code with threshold calculator
    def __init__(self, conditions=[], values=[1.]):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside

        if len(values) != len(conditions) and len(values) != len(conditions) + 1:
            raise ValueError('The number of values must be equal or one more, than the number of conditions')

        self.conditions = conditions
        self.values = values

    def calculate(self, contacts):
        for i in range(len(self.conditions)):
            #change, so that whole condition is in array
            if contacts > self.conditions[i]:
                return self.values[i]
        #the way it's built, it makes more sense to return this value, but, should check that len(value) = len(conditions) + 1 in order to return the "else" value
        return self.values[-1]

class NullDensityCalculator(DensityCalculator):
    def __init__(self):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null

    def calculate(self, contacts):
        return 1.
