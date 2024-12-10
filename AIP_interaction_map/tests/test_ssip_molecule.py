import os
import unittest
from AIP_interaction_map.ssip_molecule import SsipMolecule
import networkx as nx
import networkx.algorithms.isomorphism as iso
import logging

logging.basicConfig()
logging.disable(logging.CRITICAL)
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestSsipMolecule(unittest.TestCase):
    def setUp(self):
        self.ssip_molecule = SsipMolecule(FIXTURE_DIR + "/4llx/ssip.xml")
        
        self.test_ssip_molecule_network = nx.Graph()
        self.test_ssip_molecule_network.add_node('a1', elementType='C')
        self.test_ssip_molecule_network.add_node('a2', elementType='C')
        self.test_ssip_molecule_network.add_node('a3', elementType='C')
        self.test_ssip_molecule_network.add_node('a4', elementType='C')
        self.test_ssip_molecule_network.add_node('a5', elementType='C')
        self.test_ssip_molecule_network.add_node('a6', elementType='N')
        self.test_ssip_molecule_network.add_node('a7', elementType='C')
        self.test_ssip_molecule_network.add_node('a8', elementType='N')
        self.test_ssip_molecule_network.add_node('a9', elementType='N')
        self.test_ssip_molecule_network.add_node('a10', elementType='H')
        self.test_ssip_molecule_network.add_node('a11', elementType='H')
        self.test_ssip_molecule_network.add_node('a12', elementType='H')
        self.test_ssip_molecule_network.add_node('a13', elementType='H')
        self.test_ssip_molecule_network.add_node('a14', elementType='H')
        self.test_ssip_molecule_network.add_node('a15', elementType='H')
        self.test_ssip_molecule_network.add_node('a16', elementType='H')
        self.test_ssip_molecule_network.add_node('a17', elementType='H')
        self.test_ssip_molecule_network.add_node('a18', elementType='H')
        self.test_ssip_molecule_network.add_edge("a16", "a5")
        self.test_ssip_molecule_network.add_edge("a15", "a5")
        self.test_ssip_molecule_network.add_edge("a14", "a5")
        self.test_ssip_molecule_network.add_edge("a5", "a4")
        self.test_ssip_molecule_network.add_edge("a13", "a3")
        self.test_ssip_molecule_network.add_edge("a4", "a3")
        self.test_ssip_molecule_network.add_edge("a4", "a6")
        self.test_ssip_molecule_network.add_edge("a3", "a2")
        self.test_ssip_molecule_network.add_edge("a6", "a7")
        self.test_ssip_molecule_network.add_edge("a10", "a1")
        self.test_ssip_molecule_network.add_edge("a2", "a1")
        self.test_ssip_molecule_network.add_edge("a2", "a9")
        self.test_ssip_molecule_network.add_edge("a7", "a9")
        self.test_ssip_molecule_network.add_edge("a7", "a8")
        self.test_ssip_molecule_network.add_edge("a1", "a11")
        self.test_ssip_molecule_network.add_edge("a1", "a12")
        self.test_ssip_molecule_network.add_edge("a17", "a8")
        self.test_ssip_molecule_network.add_edge("a8", "a18")
    
    def test_network(self):
        nm = iso.categorical_node_match('elementType', 'cap')
        self.assertEqual(nx.is_isomorphic(self.test_ssip_molecule_network, self.ssip_molecule.network, node_match=nm), True)


class TestSsipAtom(unittest.TestCase):
    def setUp(self):
        self.ssip_molecule = SsipMolecule(FIXTURE_DIR + "/4llx/ssip.xml")
        self.ssip_atom = self.ssip_molecule["a3"]
        expected_value = "C.ar"
        self.assertAlmostEqual(expected_value, self.ssip_atom["type"]) 


class TestSsipSsip(unittest.TestCase):
    def setUp(self):
        self.ssip_molecule = SsipMolecule(FIXTURE_DIR + "/4llx/ssip.xml")
        self.ssip_ssip = self.ssip_molecule.get_vs_list(False)[0]

    def test_ssip_value(self):
        expected_value = -1.82
        self.assertAlmostEqual(expected_value, self.ssip_ssip.value)

    def test_neigh(self):
        expected_neigh = "a2"
        self.assertEqual(expected_neigh, self.ssip_ssip.neigh)

    def test_anchors(self):
        expected_anchor1 = "a2"
        expected_anchor2 = "a1"
        expected_anchor3 = "a3"
        self.assertEqual(expected_anchor1, self.ssip_ssip.anchor1.name)
        self.assertEqual(expected_anchor2, self.ssip_ssip.anchor2.name)
        self.assertEqual(expected_anchor3, self.ssip_ssip.anchor3.name)

    def test_weights(self):
        expected_w1 = -0.04076340602582107
        expected_w2 = -0.053912771726110864
        expected_w3 = 9.504270168114397
        self.assertAlmostEqual(expected_w1, self.ssip_ssip.weights[0])
        self.assertAlmostEqual(expected_w2, self.ssip_ssip.weights[1])
        self.assertAlmostEqual(expected_w3, self.ssip_ssip.weights[2])

