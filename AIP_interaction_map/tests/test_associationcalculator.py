"""
Script for testing of the association calculator.
@author: kate
"""

import unittest
import AIP_interaction_map.asc as associationcalculator

class TestSolvationCalculator(unittest.TestCase):
    """ test solvent association, solvation and binding energy calculation for an AIP 1/-1 pair """

    def setUp(self):
        self.aip1 = 1.0
        self.aip2 = -1.0
    
    def test_association(self):
        actual_value = associationcalculator.calculate_association_energy(self.aip1, self.aip2, solvent="chloroform")
        expected_value = -0.021166278657145154
        self.assertAlmostEqual(expected_value, actual_value, places=5)
    
    def test_solvation(self):
        actual_value = associationcalculator.calculate_solvation_energy(self.aip1, self.aip2, solvent="chloroform")
        expected_value = -2.97382219993
        self.assertAlmostEqual(expected_value, actual_value, places=5)

    def test_binding(self):
        actual_value = associationcalculator.calculate_binding_energy(self.aip1, self.aip2, solvent="chloroform")
        expected_value = -2.9949884785871452
        self.assertAlmostEqual(expected_value, actual_value, places=5)

    def test_solvents(self):
        water_value = associationcalculator.calculate_association_energy(self.aip1, self.aip2, solvent="water")
        hexadecane_value = associationcalculator.calculate_association_energy(self.aip1, self.aip2, solvent="n-hexadecane")
        cyclohexane_value = associationcalculator.calculate_association_energy(self.aip1, self.aip2, solvent="cyclohexane")
        acetonitrile_value = associationcalculator.calculate_association_energy(self.aip1, self.aip2, solvent="acetonitrile")
        toluene_value = associationcalculator.calculate_association_energy(self.aip1, self.aip2, solvent="toluene")

        water_expected = -3.3581811727062827
        hexadecane_expected = -0.28095227954893476
        cyclohexane_expected = -0.2611183164320676
        acetonitrile_expected = -0.21468583110729855
        toluene_expected = 0.01800831705787442

        self.assertAlmostEqual(water_expected, water_value, places=5)
        self.assertAlmostEqual(hexadecane_expected, hexadecane_value, places=5)
        self.assertAlmostEqual(cyclohexane_expected, cyclohexane_value, places=5)
        self.assertAlmostEqual(acetonitrile_expected, acetonitrile_value, places=5)
        self.assertAlmostEqual(toluene_expected, toluene_value, places=5)        
