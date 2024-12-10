"""
Script for testing of the ScoringTraj class.
@author: kate
"""
import os
import numpy as np
import unittest
from AIP_interaction_map.scoring_traj_hg import ScoringTraj
import mdtraj as md
from simtk.openmm import app
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestScoringTraj(unittest.TestCase):
    """ test scoring traj read in PDB / PL for simple cases """

    def setUp(self):
        #guest-host
        self.pdb1 = f"{FIXTURE_DIR}/AB/A/compA.pdb"
        self.pdb2 = f"{FIXTURE_DIR}/AB/B/compB.pdb"
        self.aip1 = f"{FIXTURE_DIR}/AB/A/ssipA.xml"
        self.aip2 = f"{FIXTURE_DIR}/AB/B/ssipB.xml"
        #protein-ligand
        self.pdb3 = f"{FIXTURE_DIR}4llx/ligand.pdb"
        self.pdb4 = f"{FIXTURE_DIR}4llx/binding_site.pdb"
        self.aip3 = f"{FIXTURE_DIR}4llx/ssip.xml"

    def test_read_in_PDB(self):
        """test scoring traj, read in pdb"""
        st = ScoringTraj(self.pdb1, self.pdb2, self.aip1, self.aip2)

        ## check for ligand name - and for ligand case, whether it's MOL or you should read it off the file
        self.assertEqual(st.ligand_name, "cpA")

        ## read in paths for pdb and ssip files consistently - especially when not explicitly stating ssip.xml path
        self.assertEqual(f"{self.pdb1}", f"{FIXTURE_DIR}/AB/A/compA.pdb")
        self.assertEqual(f"{self.pdb2}", f"{FIXTURE_DIR}/AB/B/compB.pdb")

        ## create mdtraj instance on loading in pdb
        # self.assertIsInstance(ScoringTraj, st)
        # self.assertIsInstance(md, st.mdtraj)

        ## AIP values read in correctly - check some AIP read off examples #e.g. 1st, 10th, 100th AIP values against aip.xml
        correct_AIP_values = -7.15
        self.assertEqual(correct_AIP_values, st.ssip_dict[19])
        
        ## same for fractions, isosurfaces and owners - both for guest and host
        correct_fraction_values = 0.5
        self.assertEqual(correct_fraction_values, st.fraction_dict[20])
        correct_isosurface_values = 0.002
        self.assertEqual(correct_isosurface_values, st.isosurface_dict[21])
        correct_owner_values = 2
        self.assertEqual(correct_owner_values, st.owner_dict[22])

    def test_read_in_PL(self):
        """test scoring traj, read in protein ligand"""
        st = ScoringTraj(self.pdb3, self.pdb4, self.aip3, False, protein_host=True)

        ## AIP values read in correctly - check some AIP read off examples #e.g. 1st, 10th, 100th AIP values against aip.xml
        correct_AIP_values = -1.81
        self.assertEqual(correct_AIP_values, st.ssip_dict[19])
        
        ## same for fractions, isosurfaces and owners - both for guest and host
        correct_fraction_values = 0.5
        self.assertEqual(correct_fraction_values, st.fraction_dict[20])
        correct_isosurface_values = 0.002
        self.assertEqual(correct_isosurface_values, st.isosurface_dict[21])
        correct_owner_values = 3
        self.assertEqual(correct_owner_values, st.owner_dict[22])
