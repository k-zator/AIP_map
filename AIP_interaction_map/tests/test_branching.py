"""
Script for testing of the branching algorithms.
@author: kate
"""

import unittest
from AIP_interaction_map.branch_pairing import branching, find_choices, drop_interactions, add_unique
from AIP_interaction_map.get_bipartite_pairing import branching_with_bipartite
import pandas as pd
import numpy as np
import pandas.testing as pd_testing

class TestScoringTraj(unittest.TestCase):
    """ test scoring traj read in PDB / PL for simple cases """

    def test_branch_pairing_2_3(self):
        
        AipPairs = np.array([[1, 0, 0.2, 11, 10, 0.1, 1.0, 1.0],
                            [1, 2, 0.2, 11, 12, 0.1, 1.0, 1.0],
                            [3, 2, 0.2, 13, 12, 0.1, 1.0, 1.0],
                            [3, 4, 0.2, 13, 14, 0.1, 1.0, 1.0]])
        final = []
        final_df = branching(AipPairs, final)
        correct_df = pd.DataFrame([[1.0, 0.0, 0.20, 11.0, 10.0, 0.10, 1.0],
                                   [3.0, 2.0, 0.20, 13.0, 12.0, 0.10, 1.0]],
            columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "Frac"])
        pd_testing.assert_frame_equal(correct_df, final_df)
        ## test whole branching against finding the correct answer for 2/3 and 3/3

    def test_branch_pairing_3_3(self):
        AipPairs = np.array([[1, 0, 0.2, 11, 10, 0.1, 1.0, 1.0],
                            [1, 2, 0.2, 11, 12, 0.1, 1.0, 1.0],
                            [3, 2, 0.2, 13, 12, 0.1, 1.0, 1.0],
                            [3, 4, 0.2, 13, 14, 0.1, 1.0, 1.0],
                            [5, 4, 0.2, 15, 14, 0.1, 1.0, 1.0]])
        final = []
        final_df = branching(AipPairs, final)
        correct_df = pd.DataFrame([[1.0, 0.0, 0.20, 11.0, 10.0, 0.10, 1.0],
                                   [3.0, 2.0, 0.20, 13.0, 12.0, 0.10, 1.0],
                                   [5.0, 4.0, 0.20, 15.0, 14.0, 0.10, 1.0]],
            columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "Frac"])
        pd_testing.assert_frame_equal(correct_df, final_df)

    def test_branch_uneven(self):
        AipPairs = np.array([[1, 0, 0.2, 11, 10, 0.1, 1.0, 1.0],
                            [1, 2, 0.21, 11, 12, 0.11, 1.0, 1.0],
                            [3, 2, 0.2, 13, 12, 0.1, 1.0, 1.0],
                            [3, 4, 0.21, 13, 14, 0.11, 1.0, 1.0],
                            [5, 4, 0.2, 15, 14, 0.1, 1.0, 1.0]])
        final = []
        final_df = branching(AipPairs, final)
        correct_sc = pd.DataFrame([[1.0, 0.0, 0.2, 11.0, 10.0, 0.1, 1.0],
                                   [3.0, 2.0, 0.2, 13.0, 12.0, 0.1, 1.0],
                                   [5.0, 4.0, 0.2, 15.0, 14.0, 0.1, 1.0]],
            columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "Frac"])
        pd_testing.assert_frame_equal(correct_sc, final_df)

    def test_find_choices(self):
        options1 = np.array([[1.0, 0.0, 0.2,  11.0, 10.0, 0.1,  1.0, 1.0],
                             [1.0, 2.0, 0.21, 11.0, 12.0, 0.11, 1.0, 1.0]])
        choices1 = find_choices(options1)
        correct1 = pd.DataFrame([[1.0, 0.0, 0.2,  11.0, 10.0, 0.1,  1.0, 1.0],
                                 [1.0, 2.0, 0.21, 11.0, 12.0, 0.11, 1.0, 1.0],
                                 [1.0, 0.0, 0.2,  11.0, 10.0, 0.1,  1.0, 0.5],
                                 [1.0, 2.0, 0.21, 11.0, 12.0, 0.11, 1.0, 0.5],
                                 [1.0, 0.0, 0.2 , 11.0, 10.0, 0.1 , 0.5, 1.0],
                                 [1.0, 2.0, 0.21, 11.0, 12.0, 0.11, 0.5, 1.0]])
        pd_testing.assert_frame_equal(correct1, pd.DataFrame(choices1[:-1]))
        
        options2 = np.array([[1.0, 0.0, 0.2,  11.0, 10.0, 0.1,  0.5, 1.0],
                             [1.0, 2.0, 0.21, 11.0, 12.0, 0.11, 1.0, 1.0]])
        choices2 = find_choices(options2)
        correct2 = pd.DataFrame([[1.0, 0.0, 0.2,  11.0, 10.0, 0.1,  0.5, 1.0],
                                 [1.0, 2.0, 0.21, 11.0, 12.0, 0.11, 1.0, 1.0]])
        pd_testing.assert_frame_equal(correct2, pd.DataFrame(choices2[:-1]))
    
    def test_drop_interactions(self):
        choice = np.array([1.0, 0.0, 0.2,  11.0, 10.0, 0.1,  0.5, 1.0])
        aippairs = np.array([[1, 0, 0.2, 11, 10, 0.1, 1.0, 1.0],
                             [1, 2, 0.2, 11, 12, 0.1, 1.0, 1.0],
                             [3, 2, 0.2, 13, 13, 0.1, 1.0, 1.0],
                             [3, 4, 0.2, 13, 14, 0.1, 1.0, 1.0],
                             [5, 4, 0.2, 15, 14, 0.1, 1.0, 1.0]])
        aippairs = drop_interactions(choice, aippairs)
        correct_aa = pd.DataFrame([[3.0, 2.0, 0.2, 13.0, 13.0, 0.1, 1.0, 1.0],
                                   [3.0, 4.0, 0.2, 13.0, 14.0, 0.1, 1.0, 1.0],
                                   [5.0, 4.0, 0.2, 15.0, 14.0, 0.1, 1.0, 1.0]])
        pd_testing.assert_frame_equal(correct_aa, pd.DataFrame(aippairs))
        ## drop interactions given 3 choices as given by elifs

    def test_add_unique(self):
        scen = []
        aippairs = np.array([[1, 0, 0.2, 11, 10, 0.1, 1.0, 0.5],
                             [3, 2, 0.2, 13, 13, 0.1, 1.0, 1.0],
                             [5, 4, 0.2, 15, 14, 0.1, 0.5, 0.5]])
        ar = add_unique(scen, aippairs)
        correct_aa = pd.DataFrame([[1.0, 0.0, 0.2, 11.0, 10.0, 0.1, 1.0, 0.5],
                                   [3.0, 2.0, 0.2, 13.0, 13.0, 0.1, 1.0, 1.0],
                                   [5.0, 4.0, 0.2, 15.0, 14.0, 0.1, 0.5, 0.5]])
        pd_testing.assert_frame_equal(correct_aa, pd.DataFrame(ar))
        ## add unique: A(1), B(2) 

    def test_add_unique_uneven(self):
        scen = []
        aippairs = np.array([[1, 0, 0.2, 11, 10, 0.1, 1.0, 0.5],
                             [1, 2, 0.2, 11, 12, 0.1, 1.0, 0.5],
                             [3, 4, 0.2, 13, 14, 0.1, 0.5, 1.0],
                             [5, 4, 0.2, 15, 14, 0.1, 0.5, 1.0]])
        ar = add_unique(scen, aippairs)
        correct_aa = pd.DataFrame([[1.0, 0.0, 0.2, 11.0, 10.0, 0.1, 1.0, 0.5],
                                   [3.0, 4.0, 0.2, 13.0, 14.0, 0.1, 0.5, 1.0],
                                   [5.0, 4.0, 0.2, 15.0, 14.0, 0.1, 0.5, 1.0]])
        pd_testing.assert_frame_equal(correct_aa, pd.DataFrame(ar))
        ## add unique: A 1.0 (1,2) 0.5 / A,B 0.5 (1) 1.0

    """def test_exception(self):
        "test edge cases"
        AipPairs = []
        final = []
        final_df = branching(AipPairs, final)
        correct_df = []
        self.assertEqual(correct_df, final_df)
        ## "that one exception" when I overwrite 1.0 to be 0.5 because the choice it paired with was 0.5 - check it shows up as 0.5 for a simple 1.0 - 0.5/0.5 pairing
        ## and reverse: 0.5/0.5 and 1.0 pairing"""
        ## when choice split network in two, A(1,2), B(1), C(2), depending on A choice, B or C becomes 2nd branch"""

    def test_bipartitie(self):
        AipPairs = np.array([[1, 0, 0.2, 11, 100, 0.1, 1.0, 1.0],
                             [1, 2, 0.2, 11, 120, 0.1, 1.0, 1.0],
                             [3, 2, 0.2, 13, 120, 0.1, 1.0, 1.0],
                             [3, 4, 0.2, 13, 140, 0.1, 1.0, 1.0],
                             [5, 4, 0.2, 15, 140, 0.1, 1.0, 1.0]])
        final = []
        final_df = branching_with_bipartite(AipPairs)
        correct_df = pd.DataFrame([[1.0, 0.0, 0.20, 11.0, 100.0, 0.10, 1.0],
                                   [5.0, 4.0, 0.20, 15.0, 140.0, 0.10, 1.0],
                                   [3.0, 2.0, 0.20, 13.0, 120.0, 0.10, 1.0]],
            columns=["L", "R", "Atom_Distance", "L_AIP", "R_AIP", "AIP_Distance", "Frac"])
        pd_testing.assert_frame_equal(correct_df, final_df)
        ## test whole branching against finding the correct answer for 2/3 and 3/3
