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
                "params" : {
                    "conditions" : [1,2],
                    "values" : [1,2]
                }
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)


        goldenValues = [1., 2.]
        goldenConditions = [1., 2.]

        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside)
        self.assertEqual(densityCalculator.values, goldenValues)
        self.assertEqual(densityCalculator.conditions, goldenConditions)

    def testDensityCalculatorContinuousParams(self):
        spawningBlock = {
            type: "irrelevant",
            "density" : {
                "type" : "continuous"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)


        self.assertAlmostEqual(densityCalculator.calculate(0.5, 8), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 8), 8)
        self.assertAlmostEqual(densityCalculator.calculate(0.2, 6), 1)
        self.assertAlmostEqual(densityCalculator.calculate(0.5, 6), 8)
        self.assertAlmostEqual(densityCalculator.calculate(0.05, 4), 1)
        self.assertAlmostEqual(densityCalculator.calculate(0.2, 4), 8)
