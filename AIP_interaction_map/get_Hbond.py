import numpy as np
import mdtraj as md
import pandas as pd
from AIP_interaction_map.asc import calculate_association_energy as ac

from .mdtraj.hbond import baker_hubbard

def getHbonding(self, trust_H_positions):
    """Determines hydrogen bonds between the guest and host molecules. The interacting atoms
       need to be of polar type and be O/N and H bound to O/N. It is possible to trust H positions
       and then the max D-A distance is 0.22 nm. Otherwise, the distance is made up by N/O radii."""
    # determine the aip aip separation distance
    if trust_H_positions:
        aip_aip_max_sep = 0.22  # unit: nm
    else:
        aip_aip_sep = 0.09
        d_H_0104 = 0.110  # estimate of vdW radius
        d_NO_03 = 0.155  # estimate of vdW radius
        d_NO_H = 0.1  # estimate of N/O-H bond length
        d_N_O = 0.3  # estimate of N-O max distance
        aip_aip_max_sep = d_N_O - d_NO_03 + d_NO_H + d_H_0104 + aip_aip_sep

    #meshgrids of distances between guest and host Atoms and AIPs
    L = np.array([l for l in list(self.Atom_L.keys()) if l not in self.solvent_contact])
    R = np.array([r for r in list(self.Atom_R.keys()) if r not in self.solvent_contact])
    atom_meshgrid = np.array(np.meshgrid(L, R)).T.reshape(-1, 2)
    atom_dist = md.compute_distances(self.mdTrajectory, atom_meshgrid, periodic=False)
    atom_array = atom_dist.T.reshape(len(atom_dist.T),)

    L_vs = [a.index for a in self.AIP_L.values() if a.polar]
    R_vs = [a.index for a in self.AIP_R.values() if a.polar]
    vs_meshgrid = np.array(np.meshgrid(L_vs, R_vs)).T.reshape(-1, 2)
    vs_dist = md.compute_distances(self.mdTrajectory, vs_meshgrid, periodic=False)
    vs_array = vs_dist[0]
    #select only AIP-AIPs pairings that are close enough 
    incontact = (vs_array <= aip_aip_max_sep)
    aip_L = vs_meshgrid.T[0][incontact]
    aip_R = vs_meshgrid.T[1][incontact]
    vs_dist = vs_array[incontact]
    type_L = np.array([self.AIP_L[L].type for L in aip_L])
    type_R = np.array([self.AIP_R[R].type for R in aip_R])

    #select only relevant atom types
    H_bond_correct = []
    A_A_contact = []
    for i, (L, R) in enumerate(zip(type_L, type_R)):
        if (L[0] in ["O", "N"] and R[0] in ["H"]):
            H_bond_correct.append(i)
        elif (L[0] in ["H"] and R[0] in ["O", "N"]):
            H_bond_correct.append(i)
        elif (L[0] in ["O", "N"] and R[0] in ["O", "N"]):
            A_A_contact.append(i)

    #exclude A_A contacts
    self.A_A_contact_list = list(np.unique(aip_L[A_A_contact]))+list(np.unique(aip_R[A_A_contact]))

    #hence correct H bond
    aip_L = aip_L[H_bond_correct]
    aip_R = aip_R[H_bond_correct]
    atom_L = [self.AIP_L[L].atom for L in aip_L]
    atom_R = [self.AIP_R[R].atom for R in aip_R]
    value_L = [self.AIP_L[L].value for L in aip_L]
    value_R = [self.AIP_R[R].value for R in aip_R]
    LR_dist = [atom_array[np.where((atom_meshgrid == (L, R)).all(axis=1))][0] 
               if len(np.where((atom_meshgrid == (L,R)).all(axis=1))[0]) != 0 
               else np.nan for L, R in zip(atom_L, atom_R)]
    ddG_LR = [ac(L, R, self.solvent) for L, R in zip(value_L, value_R)]

    s35_L = [self.Atom_L[L].sasa35 for L in atom_L]
    s35_R = [self.Atom_R[R].sasa35 for R in atom_R]
    s70_L = [self.Atom_L[L].sasa70 for L in atom_L]
    s70_R = [self.Atom_R[R].sasa70 for R in atom_R]
    
    #compile H_bond_df from all relevant AIP-AIP pairings
    H_bond_df = pd.DataFrame({
        "L_AIP": aip_L,
        "R_AIP": aip_R,
        "AIP_Distance": vs_dist[H_bond_correct],
        "Atom_Distance": LR_dist,
        "L": atom_L,
        "R": atom_R,
        "L_type": type_L[H_bond_correct],
        "R_type": type_R[H_bond_correct],
        "L_value": value_L,
        "R_value": value_R,
        "ddG": ddG_LR,
        "s35_L": s35_L,
        "s35_R": s35_R,
        "s70_L": s70_L,
        "s70_R": s70_R})
    H_bond_df.dropna(subset=['Atom_Distance'], inplace=True)
    H_bond_df = H_bond_df.sort_values(
        ["AIP_Distance", "Atom_Distance"], ignore_index=True)
    
    #make sure no AIPs are used for more than one AIP-AIP pairing
    already_paired_L = []
    already_paired_R = []
    for R in H_bond_df.R_AIP:
        if H_bond_df.R_AIP.value_counts()[R] == 1:
            already_paired_R.append(R)
    for i in already_paired_R:
        counter_L = H_bond_df.loc[H_bond_df.R_AIP == i].L_AIP.values
        for l in counter_L:
            H_bond_df.drop(
                H_bond_df.loc[H_bond_df.L_AIP == l].loc[H_bond_df.R_AIP != i].index, inplace=True)
    H_bond_df = H_bond_df.reset_index(drop=True)
    for L in H_bond_df.L_AIP:
        if H_bond_df.L_AIP.value_counts()[L] == 1:
            already_paired_L.append(L)
    for i in already_paired_L:
        counter_R = H_bond_df.loc[H_bond_df.L_AIP == i].R_AIP.values
        for r in counter_R:
            H_bond_df.drop(
                H_bond_df.loc[H_bond_df.R_AIP == r].loc[H_bond_df.L_AIP != i].index, inplace=True)
    H_bond_df = H_bond_df.reset_index(drop=True)

    #decide to keep based on shortest AIP-AIP dist
    for R in H_bond_df.R_AIP:
        if H_bond_df.R_AIP.value_counts()[R] > 1:
            dists = H_bond_df.loc[H_bond_df.R_AIP ==
                                  R].AIP_Distance.to_numpy()
            min_dist = min(dists)
            H_bond_df.drop(H_bond_df.loc[(H_bond_df.AIP_Distance != min_dist) & (
                H_bond_df.R_AIP == R)].index, inplace=True)
            counter_L = H_bond_df.loc[H_bond_df.R_AIP == R].L_AIP.values
            for l in counter_L:
                H_bond_df.drop(
                    H_bond_df.loc[H_bond_df.L_AIP == l].loc[H_bond_df.R_AIP != R].index, inplace=True)
    for L in H_bond_df.L_AIP:
        if H_bond_df.L_AIP.value_counts()[L] > 1:
            dists = H_bond_df.loc[H_bond_df.L_AIP ==
                                  L].AIP_Distance.to_numpy()
            min_dist = min(dists)
            H_bond_df.drop(H_bond_df.loc[(H_bond_df.AIP_Distance != min_dist)
                                         & (H_bond_df.L_AIP == L)].index, inplace=True)
            counter_R = H_bond_df.loc[H_bond_df.L_AIP == L].R_AIP.values
            for r in counter_R:
                H_bond_df.drop(
                    H_bond_df.loc[H_bond_df.R_AIP == r].loc[H_bond_df.L_AIP != L].index, inplace=True)
                
    self.H_bond_df = H_bond_df.reset_index(drop=True)
    self.H_bond_df["Frac"] = 1
    #sc = scaled by AIP fraction
    self.H_bond_df["ddG_sc"] = self.H_bond_df["ddG"]

def getRRHbonding(self):
    """If the host is much bigger, it could have intramolecular hydrogen bonding
        and those could interfere with non-polar interactions if are not excluded."""
    atom_traj = self.mdTrajectory.atom_slice(self.state.residue_all_indices)
    hbond_array = baker_hubbard(
        atom_traj, periodic=False, distance_cutoff=0.33, angle_cutoff=120, heavy_atom_distance=True)
    
    #now, the indices are shifted by len(self.state.ligand_all_indices)
    H_indices = hbond_array.T[1] + len(self.state.ligand_all_indices)
    NO_indices = hbond_array.T[2] + len(self.state.ligand_all_indices)
    atom_indices = np.concatenate((H_indices, NO_indices))
    vs_indices = [self.AtomAipDict[i] for i in atom_indices if i in self.AtomAipDict.keys()]
    vs_indices = [item for sublist in vs_indices for item in sublist]

    #check the indices are indeed polar
    vs_indices_polar = [i for i in vs_indices if (
        self.AIP_R[i].polar == True) or (self.AIP_R[i].polar == True)]
    self.RR_H_bond_Atom = np.array(atom_indices)
    self.RR_H_bond_AIP = np.array(vs_indices_polar)

    #and it's presence in these lists that will exclude them from interacting downstream
    #define_atom_aip_subset, line 75
