from simtk.openmm import app
import os
import unittest
import numpy as np
import tempfile
import shutil
from AIP_interaction_map.create_bs import create_parser, create_bs
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')


ligand = f"{FIXTURE_DIR}1s38/ligand_modified.pdb"
protein = f"{FIXTURE_DIR}1s38/1s38_all_flare.pdb"
pocket = f"{FIXTURE_DIR}1s38/1s38_pocket.pdb"
ff = f"{FIXTURE_DIR}just_protein_ff.xml"


class TestCreateFfBsTemplate(unittest.TestCase):
    def setUp(self):
        self.tempfile = tempfile.mkdtemp()
        self.bs_file = f"{self.tempfile}/binding_site.pdb"
        arguments = create_parser().parse_args([
            '-p', protein, 
            '-o', self.tempfile,
            'template',
            '-t', pocket,
            ])
        create_bs(arguments)
        self.assertTrue(os.path.isfile(self.bs_file))
        self.bs = app.PDBFile(self.bs_file)

    def test_binding_site(self):
        self.modeller = app.Modeller(self.bs.topology, self.bs.getPositions())
        self.forcefield = app.ForceField(ff)
        self.modeller.addExtraParticles(self.forcefield, ignoreExternalBonds=True)

    def cleanUp(self):
        shutil.rmtree(self.tempfile)

class TestCreateFfBsDistance(unittest.TestCase):
    def setUp(self):
        self.tempfile = tempfile.mkdtemp()
        self.bs_file = f"{self.tempfile}/binding_site.pdb"
        arguments = create_parser().parse_args([
            '-p', protein, 
            '-o', self.tempfile,
            'distance',
            '-l', ligand 
            ])
        create_bs(arguments)
        self.assertTrue(os.path.isfile(self.bs_file))
        self.bs = app.PDBFile(self.bs_file)

    def test_binding_site(self):
        self.modeller = app.Modeller(self.bs.topology, self.bs.getPositions())
        self.forcefield = app.ForceField(ff)
        self.modeller.addExtraParticles(self.forcefield, ignoreExternalBonds=True)

    def cleanUp(self):
        shutil.rmtree(self.tempfile)
