import numpy as np
import mdtraj as md
import pandas as pd
from AIP_interaction_map.AtomAIPClasses import Atom, AIP
from AIP_interaction_map.scoring_traj_hg import ScoringTraj
from AIP_interaction_map.get_desolvation import getDesolvation
from AIP_interaction_map.get_sasa import solvent_contact_search, exclude_solvent_contact
from AIP_interaction_map.get_Hbond import getHbonding, getRRHbonding
from AIP_interaction_map.get_non_polar import getNonPolar
from AIP_interaction_map.get_atom_aip_lists import getAtomAIPclasslists, getAtomAipDict

class Scoring():
    """Overall class for describing guest-host (ligand-binding site) interactions,
       considering both atoms and AIPs in a hybrid fashion.
       It takes two obligatory paths: guest_name and host_name - the names of molecules in the complex.
       It also requires paths to aip files for both but those are assumed to be located in the same
       directory and named ssip.xml (unless specified).
       Other possible arguments include maximum aip-aip distance, SASA probe radius, solvent,
       bypass at 0.8 A (which treats those closer contacts are preferential, hence shortening
       the time of calculation), SASA fraction required to desolvate, extension in aip distance upon
       desolvation, choice of branching algorithm searching for AIP pairings, and finally protein_host
       which, if True, uses precompiled aip values for protein host (hence uses OpenMM)."""

    def __init__(self, guest_name, host_name, aip_guest=False, aip_host=False,
                 protein_host=False, branch=True,
                 max_aip_dist=0.15, bypass_at1=False, solvent="chloroform",
                 probe_radius=0.2, frac_to_desolv=0.5, ext_upon_desolv=0.03):

        self.max_aip_dist = max_aip_dist
        self.solvent = solvent
        self.probe_radius = probe_radius
        self.frac_to_desolv = frac_to_desolv
        self.bypass_at1 = bypass_at1
        self.branch = branch
        self.ext_upon_desolv = ext_upon_desolv

        #Scoring Traj creates mdtraj trajectory from input files
        self.state = ScoringTraj(guest_name, host_name, aip_guest, aip_host, protein_host)
        if self.state.ligand_name == "error":
            return
        else:
            self.mdTrajectory = self.state.mdtraj

        #as Atoms and AIPs are linked, dictionary to convert between them
        getAtomAIPclasslists(self, protein_host)
        getAtomAipDict(self)
        #solvent accessible surface area measurements 
        solvent_contact_search(self)
        #exclude_solvent_contact(self)
        #first step is H bonding as their AIP contacts concern specific atom types
        getHbonding(self, trust_H_positions=True)
        #and similar search in protein to avoid double counting their interactions
        getRRHbonding(self)
        #now onto non-polar interactions
        exclude_solvent_contact(self)
        self.non_polar_df = getNonPolar(self)

        if len(self.non_polar_df) > 0:
            self.interaction_df = pd.concat(
                [self.non_polar_df, self.H_bond_df],  ignore_index=True)
            self.dG = sum(self.interaction_df.ddG_sc)
        else:
            self.interaction_df = self.H_bond_df
            self.dG = sum(self.interaction_df.ddG)

        #if there are not enough close distance AIP contacts, they could be found by
        #looking at desolvated AIPs as the source
        # getDesolvation(self)

    @staticmethod
    def split_df_row(df, target_col):
        '''origin: https://gist.github.com/jlln/338b4b0b55bd6984f883 @kleinias
           It splits row that has multiple entries in a cell to correspond to
           multiple rows with one entry in that cell'''

        row_accumulator = []
        def split_list_to_rows(row):
            split_row = row[target_col]
            if isinstance(split_row, list):
                for s in split_row:
                    new_row = row.to_dict()
                    new_row[target_col] = s
                    row_accumulator.append(new_row)
                if split_row == []:
                    new_row = row.to_dict()
                    new_row[target_col] = None
                    row_accumulator.append(new_row)
            else:
                new_row = row.to_dict()
                new_row[target_col] = split_row
                row_accumulator.append(new_row)
        df.apply(split_list_to_rows, axis=1)
        new_df = pd.DataFrame(row_accumulator)
        return new_df
