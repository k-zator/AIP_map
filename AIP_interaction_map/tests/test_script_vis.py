"""
Script for testing of the association calculator.
@author: kate
"""
import os
import unittest
from AIP_interaction_map.score import create_parser
from AIP_interaction_map.getScoring import Scoring
from AIP_interaction_map.create_jmol_vis import create_ip_vis
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestScoringScript(unittest.TestCase):
    """ test scoring script """

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
    
    def test_script_paths(self):
        p = create_parser()
        args = p.parse_args(["--PDB_file1", self.pdb1,
                             "--PDB_file2", self.pdb2,
                             "--aip_file1", self.aip1,
                             "--aip_file2", self.aip2
                             ])
        self.assertEqual(args.PDB_file1, self.pdb1)
        self.assertEqual(args.PDB_file2, self.pdb2)
        self.assertEqual(args.aip_file1, self.aip1)
        self.assertEqual(args.aip_file2, self.aip2)
        ## check conversion of pdb/aip paths dependent on input - with/out explicit aip path 

    def test_script_paths(self):
        p = create_parser()
        args = p.parse_args(["--PDB_file1", self.pdb1,
                             "--PDB_file2", self.pdb2,
                             "--aip_file1", self.aip1,
                             "--aip_file2", self.aip2,
                             "--write", "A_B"
                             ])
        self.assertEqual(args.PDB_file1, self.pdb1)
        self.assertEqual(args.PDB_file2, self.pdb2)
        self.assertEqual(args.aip_file1, self.aip1)
        self.assertEqual(args.aip_file2, self.aip2)
        ## check write variation to paths

    """def test_script_variable(self):
        p = create_parser()
        args = p.parse_args(["--PDB_file1", self.pdb1,
                             "--PDB_file2", self.pdb2,
                             "--aip_file1", self.aip1,
                             "--aip_file2", self.aip2,
                             "--max_aip_dist", 0.18,
                             "--solvent", "water",
                             "--mbpa", False,
                             "--protein", False,
                             "--bypass_at1", False
                            ])
        self.assertEqual(args.max_aip_dist, 0.18)
        self.assertEqual(args.solvent, "chloroform")
        self.assertEqual(args.mbpa, False)
        self.assertEqual(args.protein, False)
        self.assertEqual(args.aip2, False)
        self.assertEqual(args.bypass_at1, False)        
        ## check, max_aip_dist, solvent, mbpa, protein, bypass,"""
        
class TestVisualisation(unittest.TestCase):
    """ test scoring script """

    def setUp(self):
        self.pdb3 = f"{FIXTURE_DIR}4llx/ligand.pdb"
        self.pdb4 = f"{FIXTURE_DIR}4llx/binding_site.pdb"
        self.aip3 = f"{FIXTURE_DIR}4llx/ssip.xml"
        self.correct_jmol = f"{FIXTURE_DIR}4llx/4llx_script"
    
    def test_visualisation(self):
        k = Scoring(self.pdb3, self.pdb4, self.aip3, "water", protein_host=True, max_aip_dist=0.18)
        create_ip_vis(k, "new_jmol", f"{FIXTURE_DIR}")
        self.assertListEqual(list(open(f"{FIXTURE_DIR}/new_jmol_script")), list(open(self.correct_jmol)))
        ## line by line