import blockNames

class DensityBuilder():
    def build(self, spawningBlock):
        try:
            densityBlock = spawningBlock[blockNames.SpawningParams.density]
        except KeyError:
            print "Using null density calculator (no preference for any cluster)"
            return NullDensityCalculator()

class DensityCalculator():
    #Mostly duplicated code with threshold calculator
    def __init__(self, conditions=[], values=[1.]):

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
    def calculate(seld, contacts):
        return 1.
