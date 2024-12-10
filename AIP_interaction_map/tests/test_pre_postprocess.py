"""
Script for testing of the preprocessing scripts: PDB_preprocess and create_combined.
@author: kate
"""
import os
import unittest
from  AIP_interaction_map.create_combined import create_combined
from  AIP_interaction_map.PDB_preprocess import PDB_preprocess
from AIP_interaction_map.getScoring import Scoring
FIXTURE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_files', '')

class TestCombined(unittest.TestCase):
    """ combined and DataDicts created the same as already present example files """

    def setUp(self):
        self.correct_combined =  f"{FIXTURE_DIR}AB/A/correct_combined.pdb"
        self.correct_DataDicts =  f"{FIXTURE_DIR}/AB/A/correct_DataDicts.txt"
        self.pdb1 = f"{FIXTURE_DIR}/AB/A/compA.pdb"
        self.aip1 = f"{FIXTURE_DIR}/AB/A/ssipA.xml"

    def test_combined(self):
        create_combined(self.pdb1, self.aip1)
        self.assertListEqual(list(open(self.correct_combined)), list(open(f"{FIXTURE_DIR}/AB/A/combined.pdb")))
        self.assertListEqual(list(open(self.correct_DataDicts)), list(open(f"{FIXTURE_DIR}AB/A/DataDicts.txt")))

class TestPDBPreproces(unittest.TestCase):
    """ test against correct formatting (as it overwrites it, copy the file first) """

    def setUp(self):
        self.correct_format =  f"{FIXTURE_DIR}/AB/A/compA.pdb"
        # original_format = copy to new name
        self.PDB_file =  f"{FIXTURE_DIR}/AB/A/compA_pre_compA.pdb"

    def test_preprocess(self):
        PDB_preprocess(self.PDB_file)
        self.assertListEqual(list(open(self.correct_format)), list(open(self.PDB_file)))

class TestCSVOutput(unittest.TestCase):
    """ test against correct formatting (as it overwrites it, copy the file first) """

    def setUp(self):
        self.correct_format =  f"{FIXTURE_DIR}/AB/correct_output.csv"
        self.pdb1 = f"{FIXTURE_DIR}/AB/A/compA.pdb"
        self.pdb2 = f"{FIXTURE_DIR}/AB/B/compB.pdb"
        self.aip1 = f"{FIXTURE_DIR}/AB/A/ssipA.xml"
        self.aip2 = f"{FIXTURE_DIR}/AB/B/ssipB.xml"

    def test_csvoutput(self):
        k = Scoring(self.pdb1, self.pdb2, self.aip1, self.aip2, max_aip_dist=0.18, solvent="chloroform")
        output_df = k.interaction_df[["L","L_type","L_value", "R", "R_type", "R_value", "Frac", "ddG_sc"]]
        output_df = output_df.round(1)
        output_df.L += 1
        output_df.R -= len(k.state.ligand_vs_indices) - 2
        output_df = output_df.rename(columns={"L": f"A", "R": f"B",
                                "L_type": f"A_type", "R_type": f"B_type",
                                "L_value": f"A_value", "R_value": f"B_value",
                                "L_frac": "f", "ddG_sc": "ddG (kJ/mol)"})
        output_df.to_csv(f"{FIXTURE_DIR}/AB/AB.csv")
        self.assertListEqual(list(open(self.correct_format)), list(open(f"{FIXTURE_DIR}/AB/AB.csv")))