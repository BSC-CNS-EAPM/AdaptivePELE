import blockNames
import thresholdcalculatortypes

#make test
class ThresholdCalculatorBuilder():
    def build(self, clusteringBlock):
        try:
            thresholdCalculatorBlock = clusteringBlock[blockNames.ClusteringTypes.thresholdCalculator]
        except KeyError:
            #Default value if no threshold calculator block was defined
            return ThresholdCalculatorHeaviside()

        try:
            type = thresholdCalculatorBlock[blockNames.ThresholdCalculator.type]
        except KeyError: 
            sys.exit("Threshold calculator must have a type")

        if type == blockNames.ThresholdCalculator.constant:
            try:
                value = thresholdCalculatorBlock[blockNames.ThresholdCalculator.value]
                return ThresholdCalculatorConstant(value)
            except KeyError:
                print "Using default parameters for constant threshold calculator"
                return ThresholdCalculatorConstant()
        elif type == blockNames.ThresholdCalculator.heaviside:
            try:
                values = thresholdCalculatorBlock[blockNames.ThresholdCalculator.values]
                conditions = thresholdCalculatorBlock[blockNames.ThresholdCalculator.conditions]
                return ThresholdCalculatorHeaviside(conditions, values)
            except KeyError:
                print "Using default parameters for Heaviside threshold calculator"
                return ThresholdCalculatorHeaviside()
        else:
            sys.exit("Unknown threshold calculator type! Choices are: " + str(thresholdcalculatortypes.THRESHOLD_CALCULATOR_TYPE_TO_STRING_DICTIONARY.values()))


from abc import ABCMeta, abstractmethod
class ThresholdCalculator():
    def __init__(self):
        self.type = "BaseClass" #change for abstract attribute
        pass

    @abstractmethod
    def calculate(self, contacts):
        pass

class ThresholdCalculatorConstant(ThresholdCalculator):
    def __init__(self, value = 2):
        self.type = thresholdcalculatortypes.THRESHOLD_CALCULATOR_TYPES.constant
        self.value = value

    def caclulate(self, contacts):
        return self.value

class ThresholdCalculatorHeaviside(ThresholdCalculator):
    def __init__(self, conditions=[15,10], values=[2,3,4]):
        self.type = thresholdcalculatortypes.THRESHOLD_CALCULATOR_TYPES.heaviside
    
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

