"""
Script for testing of the Scoring class.
@author: kate
"""
import os
import unittest
from AIP_interaction_map.getScoring import Scoring
from AIP_interaction_map.AtomAIPClasses import Atom, AIP
import mdtraj as md
from simtk.openmm import app
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestScoring(unittest.TestCase):
    """ test scoring """

    def setUp(self):
        #guest-host 1
        self.pdb1 = f"{FIXTURE_DIR}/AB/A/compA.pdb"
        self.pdb2 = f"{FIXTURE_DIR}/AB/B/compB.pdb"
        self.aip1 = f"{FIXTURE_DIR}/AB/A/ssipA.xml"
        self.aip2 = f"{FIXTURE_DIR}/AB/B/ssipB.xml"
        #protein-ligand
        self.pdb5 = f"{FIXTURE_DIR}/4llx/ligand.pdb"
        self.pdb6 = f"{FIXTURE_DIR}/4llx/binding_site.pdb"
        self.aip5 = f"{FIXTURE_DIR}/4llx/ssip.xml"

    def test_atom_aip_classes(self):
        s = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2)

        s.Atom_L[0].set_sasa35(1)
        s.Atom_L[0].set_sasa70(2)
        correct_values = [1,2]
        self.assertEqual(correct_values[0], s.Atom_L[0].sasa35)
        self.assertEqual(correct_values[1], s.Atom_L[0].sasa70)

    def test_getHbonding(self):
        s = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2)
        correct_Hbond = [0.05869999900460243, 0.06750000268220901, 0.07050000131130219, 0.07339999824762344]
        self.assertCountEqual(correct_Hbond, list(s.H_bond_df.AIP_Distance.round(4)))
        ## test system of N/HN, O/HN/HN, O/OH, at correct proximity
        ## and at larger, snd also O/O for failure
        ## and N/OH/OH so that it has to pick one

    def test_RRHbonding(self):
        pass
        # s = Scoring(self.pdb5, self.pdb6, self.aip5)
        # correct_remaining_RR_H_bond_Atom = []
        # correct_remaining_RR_H_bond_AIP = []
        # self.assertEqual(correct_remaining_RR_H_bond_Atom, s.AipPairs)
        # self.assertEqual(correct_remaining_RR_H_bond_AIP, s.AipPairs)        
        ## for simple binding site, check the correct AIPs remain after Baker-Hubbard

    def test_getScoring(self):
        s = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2, bypass_at1=False, max_aip_dist=0.18,)
        correct_dist = [0.0582, 0.0932, 0.1682, 0.1728, 0.0694, 0.1421, 0.1447, 0.0735, 0.0926, 0.13, 0.1796, 0.0743, 0.0938, 0.0801, 0.1236, 0.0587, 0.0675, 0.0705, 0.0734]
        self.assertAlmostEqual(correct_dist[0], list(s.interaction_df.AIP_Distance.round(4))[0])

        s = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2, bypass_at1=True, max_aip_dist=0.18,)
        correct_dist = [0.1266, 0.1380, 0.1682, 0.1728, 0.1796, 0.05824, 0.0694, 0.0734, 0.0743, 0.0800, 0.0587, 0.0675, 0.0705, 0.0734]
        self.assertAlmostEqual(correct_dist[0], list(s.interaction_df.AIP_Distance.round(4))[0])
        ## test effect of bypass - for system where it would and wouldn't make a difference

    def test_SASA_calculations(self):
        s = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2, max_aip_dist=0.18,)
        v = [v.sasa35.round(5) for v in s.Atom_L.values()]
        correct_sasa = [0.05757, 0.04922, 0.05217, 0.04809, 0.07182, 0.03433, 0.07288, 
                        0.03697, 0.02174, 0.02325, 0.04921, 0.08816, 0.05857, 0.07668]
        self.assertAlmostEqual(correct_sasa[0], v[0])
        ## calculate SASA for guest, host, free and bound

    def test_AtomAipPairs(self):
        s = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2, max_aip_dist=0.18,)
        correct_L = [32.0, 33.0, 33.0, 35.0, 25.0, 24.0, 31.0, 22.0, 23.0, 30.0, 34.0, 15.0, 18.0, 19.0, 14.0]
        correct_R = [246.0, 254.0, 252.0, 212.0, 196.0, 203.0, 194.0, 262.0, 267.0, 265.0, 216.0, 273.0, 302.0, 272.0, 303.0]
        self.assertCountEqual(correct_L, list(s.interaction_df.L_AIP))    
        self.assertCountEqual(correct_R, list(s.interaction_df.R_AIP))    
        ## test get_atom_aip_pairing_lists for AipPairs against known set