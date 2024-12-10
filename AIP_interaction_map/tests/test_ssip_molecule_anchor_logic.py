import os
import unittest
from AIP_interaction_map.ssip_molecule import SsipMolecule
import logging

logging.basicConfig()
logging.disable(logging.CRITICAL)

FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestAnchorLogicMethane(unittest.TestCase):
    def setUp(self):
        self.assertRaises(ValueError, SsipMolecule(FIXTURE_DIR + "/methane_ssip.xml"))
        #This should give an error because I have set it up so that an SSIPs 
        #need two heavy atoms for the anchor definition 
    
class TestAnchorLogicEthane(unittest.TestCase):
    def setUp(self):
        self.ssip_molecule = SsipMolecule(FIXTURE_DIR + "/ethane_ssip.xml")
    
    def test_anchors(self):
        self.ssip_ssip = self.ssip_molecule.get_vs_list(False)[0]
        self.assertEqual(self.ssip_ssip.anchor1.name, "a5")
        self.assertEqual(self.ssip_ssip.anchor2.name, "a1")
        self.assertEqual(self.ssip_ssip.anchor3.name, "a2")

class TestAnchorLogicHCN(unittest.TestCase):
    def setUp(self):
        self.assertRaises(ValueError, SsipMolecule(FIXTURE_DIR + "/hcn_ssip.xml"))
        #This should give a value error because three atoms in a line
        #cannot describe SSIPs in 3D space (no cross product) 

class TestAnchorLogicMeCN(unittest.TestCase):
    def setUp(self):
        self.ssip_molecule = SsipMolecule(FIXTURE_DIR + "/mecn_ssip.xml")
    
    def test_anchors(self):
        self.ssip_ssip = self.ssip_molecule.get_vs_list(False)[0]
        self.assertEqual(self.ssip_ssip.anchor1.name, "a3")
        self.assertEqual(self.ssip_ssip.anchor2.name, "a2")
        self.assertEqual(self.ssip_ssip.anchor3.elem, "H")
