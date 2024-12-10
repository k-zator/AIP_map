from simtk.openmm import app
import os
import unittest
import numpy as np
import tempfile
import shutil
from AIP_interaction_map.create_ff import create_parser, create_ff
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')


ligand_pdb = f"{FIXTURE_DIR}/4llx/ligand.pdb"
ligand_ssip = f"{FIXTURE_DIR}/4llx/ssip.xml"

class TestParser(unittest.TestCase):
    def setUp(self):
        self.parser = create_parser()

    def test_parsing_vs(self):
        parsed = self.parser.parse_args()
        self.assertEqual(parsed.vs, False)
        parsed = self.parser.parse_args(['--vs'])
        self.assertEqual(parsed.vs, True)

class TestCreateFfVs(unittest.TestCase):
    def setUp(self):
        self.parser = create_parser()
        self.tempfile = tempfile.mkdtemp()
        self.ligand_ff = f"{self.tempfile}/ligand_vs_ff.xml"
        self.custom_ff = f"{self.tempfile}/ligand_vs_custom.xml"
        parsed = self.parser.parse_args([
                '--out', f"{self.tempfile}",\
                '--ssip', ligand_ssip,
                '--vs', '-l'
                ]
            )

        create_ff(parsed)
        self.ligand = app.PDBFile(ligand_pdb)
        self.modeller = app.Modeller(self.ligand.topology, self.ligand.getPositions())
        self.forcefield = app.ForceField(self.ligand_ff, self.custom_ff)
        self.modeller.addExtraParticles(self.forcefield, ignoreExternalBonds=True)
        self.system = self.forcefield.createSystem(self.modeller.topology, ignoreExternalBonds=True)

    
    def test_file_creation(self):
        """
        Tests whether the ff files created by create_ff
        can be used to add the Virtual Sites onto the PDB file.
        """
        self.assertTrue(os.path.isfile(self.ligand_ff))
        self.assertTrue(os.path.isfile(self.custom_ff))

    def test_vs(self):
        self.assertEqual(self.system.getNumParticles(), 43)
        self.assertFalse(self.system.isVirtualSite(11))
        self.assertTrue(self.system.isVirtualSite(32))
    
    def test_SA_info(self):
        self.assertEqual(self.system.getForce(0).getGlobalParameterName(0), "SA")
        self.assertEqual(round(self.system.getForce(0).getGlobalParameterDefaultValue(0)), 165)
        self.assertEqual(self.system.getForce(0).getGlobalParameterName(1), "SApos")
        self.assertEqual(round(self.system.getForce(0).getGlobalParameterDefaultValue(1)), 93)
        self.assertEqual(self.system.getForce(0).getGlobalParameterName(2), "SAneg")
        self.assertEqual(round(self.system.getForce(0).getGlobalParameterDefaultValue(2)), 72)

    def test_zeros_info(self):
        self.assertEqual(self.system.getForce(0).getGlobalParameterName(3), "n_zeros")
        self.assertEqual(round(self.system.getForce(0).getGlobalParameterDefaultValue(3)), 0)
    
    def cleanUp(self):
        shutil.rmtree(self.tempfile)
