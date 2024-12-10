import os
import unittest
from AIP_interaction_map.binding_site import BindingSite
import numpy as np

FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestSsipMoleculeOneChain(unittest.TestCase):
    def setUp(self):
        pock = FIXTURE_DIR + "/4llx/4llx_pocket.pdb"
        prot = FIXTURE_DIR  + "/4llx/4llx_protein.pdb"
        self.binding_site = BindingSite(prot, pock)

    def test_one_chain(self):
        self.assertEqual(self.binding_site.traj_prot.topology.n_chains, 1)

    def test_pock_one_chain(self):
        list_res = []
        for residue in self.binding_site.traj.topology.residues:
            list_res.append(residue)
        
        self.assertEqual(len(list_res), 1)
        residue = list_res[0]
        
        self.assertEqual(residue.name, "TYR")
        
        self.assertEqual(residue.n_atoms, 21)
        
        self.assertEqual(residue.atom(12).name, "HA")
        
        expected_xyz = 0.1*np.array([13.504, 14.514, 48.183])
        self.assertEqual(np.allclose(self.binding_site.traj.xyz[0][12], expected_xyz), True)


class TestSsipMoleculeMultiChain(unittest.TestCase):
    def setUp(self):
        pock = FIXTURE_DIR + "/pocket_multi_chain.pdb"
        prot = FIXTURE_DIR  + "/protein_multi_chain.pdb"
        self.binding_site = BindingSite(prot, pock)

    def test_multi_chain(self):
        self.assertGreater(self.binding_site.traj_prot.topology.n_chains, 1)

    def test_pock_multi_chain(self):
        list_res = []
        for residue in self.binding_site.traj.topology.residues:
            list_res.append(residue)
        
        self.assertEqual(len(list_res), 1)
        residue = list_res[0]
        
        self.assertEqual(residue.name, "VAL")

        self.assertEqual(residue.n_atoms, 16)

        self.assertEqual(residue.atom(7).name, "HA")

        expected_xyz = 0.1*np.array([9.632, 10.545, 13.891])
        self.assertEqual(np.allclose(self.binding_site.traj.xyz[0][7], expected_xyz), True)
