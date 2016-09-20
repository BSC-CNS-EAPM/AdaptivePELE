from spawning import densitycalculator
from spawning import densitycalculatortypes
import unittest

class densityCalculatorTest(unittest.TestCase):
    def testDensityCalculatorUndefinedBlock(self):
        spawningBlock = {}
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)
        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null)

    def testDensityCalculatorNullCalculator(self):
        spawningBlock = { 
            type: "irrelevant",
            "density" : {
                "type" : "null"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)
        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null)

    def testDensityCalculatorHeavisideNoParams(self):
        spawningBlock = { 
            type: "irrelevant",
            "density" : {
                "type" : "heaviside"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)


        goldenValues = [1.]
        goldenConditions = []

        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside)
        self.assertEqual(densityCalculator.values, goldenValues)
        self.assertEqual(densityCalculator.conditions, goldenConditions)

    def testDensityCalculatorHeavisideParams(self):
        spawningBlock = { 
            type: "irrelevant",
            "density" : {
                "type" : "heaviside",
                "conditions" : [1,2],
                "values" : [1,2]
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)


        goldenValues = [1., 2.]
        goldenConditions = [1., 2.]

        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside)
        self.assertEqual(densityCalculator.values, goldenValues)
        self.assertEqual(densityCalculator.conditions, goldenConditions)
