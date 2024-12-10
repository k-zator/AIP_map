import numpy as np
import mdtraj as md
import pandas as pd
from AIP_interaction_map.AtomAIPClasses import Atom, AIP

def getAtomAipDict(self):
    """Create dictionary that matches Atoms and AIPs in both guest and host compounds."""
    self.AtomAipDict = {}
    for k, v in self.AIP_L.items():
        if v.atom in self.AtomAipDict.keys():
            self.AtomAipDict[v.atom].append(k)
        else:
            self.AtomAipDict[v.atom] = [k]
    for k, v in self.AIP_R.items():
        if v.atom in self.AtomAipDict.keys():
            self.AtomAipDict[v.atom].append(k)
        else:
            self.AtomAipDict[v.atom] = [k]

def getAtomAIPclasslists(self, protein_host):
    """Create dictionaries of Atom and AIP with index as key and class objects as values."""
    AIP_L = {}
    for i in self.state.ligand_vs_indices:
        if protein_host == False:
            atom_type = self.state.atom_type_dict[i]
        else:
            atom_type = self.mdTrajectory.topology.atom(i).name.split("_", 1)[1]
        if self.state.dual:
            aip = AIP(i, self.state.ssip_dict[i], self.state.owner_dict[i], atom_type,
                    self.mdTrajectory.xyz[0][i], self.state.isosurface_dict[i], self.state.fraction_dict[i], 
                    self.state.dual_dict[i])
        else:
            aip = AIP(i, self.state.ssip_dict[i], self.state.owner_dict[i], atom_type,
                    self.mdTrajectory.xyz[0][i], self.state.isosurface_dict[i], self.state.fraction_dict[i])

        AIP_L[i] = aip
    self.AIP_L = AIP_L

    AIP_R = {}
    for i in self.state.residue_vs_indices:
        if protein_host == False:
            atom_type = self.state.atom_type_dict[i]
        else:
            atom_type = self.mdTrajectory.topology.atom(i).name.split("_", 1)[1]
        if self.state.dual:
            aip = AIP(i, self.state.ssip_dict[i], self.state.owner_dict[i], atom_type,
                    self.mdTrajectory.xyz[0][i], self.state.isosurface_dict[i], self.state.fraction_dict[i], 
                    self.state.dual_dict[i])
        else:                
            aip = AIP(i, self.state.ssip_dict[i], self.state.owner_dict[i], atom_type,
                    self.mdTrajectory.xyz[0][i], self.state.isosurface_dict[i], self.state.fraction_dict[i])
        AIP_R[i] = aip
    self.AIP_R = AIP_R

    Atom_L = {}
    for i in self.state.ligand_atom_indices:
        try:
            atom_type = [j.type for j in self.AIP_L.values() if j.atom == i][0]
        except:
            atom_type = "X"
        atom = Atom(i, atom_type, self.mdTrajectory.xyz[0][i])
        Atom_L[i] = atom
    self.Atom_L = Atom_L

    Atom_R = {}
    for i in self.state.residue_atom_indices:
        try:
            atom_type = [j.type for j in self.AIP_R.values() if j.atom == i][0]
        except:
            atom_type = "X"
        atom = Atom(i, atom_type, self.mdTrajectory.xyz[0][i])
        Atom_R[i] = atom
    self.Atom_R = Atom_R

def define_atom_aip_subset(self):
    """Remove already interacting species (by H bonding) from lists of AIPs and Atoms."""
    L_OH_in_H_bond_df = [l for l in self.H_bond_df.L if self.Atom_L[l].type[0] != "N"]
    R_OH_in_H_bond_df = [l for l in self.H_bond_df.R if self.Atom_R[l].type[0] != "N"]
    L_Atom = [l for l in list(self.Atom_L.keys()) if l not in 
              (L_OH_in_H_bond_df+self.solvent_contact)]
    L_AIP = [l for l in list(self.AIP_L.keys()) if l not in 
             (list(self.H_bond_df.L_AIP)+self.A_A_contact_list+self.solvent_contact)]
    R_Atom = [r for r in list(self.Atom_R.keys()) if r not in 
              (R_OH_in_H_bond_df+self.solvent_contact)]
    if hasattr(self, "RR_H_bond_AIP"):
        R_AIP = [r for r in list(self.AIP_R.keys()) if r not in 
                 (list(self.H_bond_df.R_AIP)+list(self.RR_H_bond_AIP) \
                  +self.A_A_contact_list+self.solvent_contact)]
    else:
        R_AIP = [r for r in list(self.AIP_R.keys()) if r not in list(self.H_bond_df.R_AIP)+self.solvent_contact]
    return L_Atom, L_AIP, R_Atom, R_AIP

def get_atom_aip_pairing_lists(self, SASA_pass, L_Atom, L_AIP, R_Atom, R_AIP):
    """Create dataframe with AIP-AIP (atom respective Atom) pairings to be used 
       by the branching algorithm, subject to the max AIP-AIP distance criterion."""
    #getDesolvation SASA pass with a larger atom_atom_distance
    if SASA_pass == False:
        L_Atom, L_AIP, R_Atom, R_AIP = define_atom_aip_subset(self)
        atom_atom_dist = 0.09 # estimate of solvent diameter in nm
    else:
        atom_atom_dist = 0.09 + 0.015

    #Atom-Atom distances are first pass, pairing need to be at least van der Waals
    #radii + solvent diameter close, otherwise, they could be captured by SASA pass
    atom_meshgrid = np.array(np.meshgrid(L_Atom, R_Atom)).T.reshape(-1, 2)
    atom_dist = md.compute_distances(self.mdTrajectory, atom_meshgrid, periodic=False)
    atom_array = atom_dist[0]

    vs_meshgrid = np.array(np.meshgrid(L_AIP, R_AIP)).T.reshape(-1, 2)
    vs_dist = md.compute_distances(self.mdTrajectory, vs_meshgrid, periodic=False)
    vs_array = vs_dist[0]

    atom_contact = []
    for L, R in atom_meshgrid:
        L_r = self.mdTrajectory.topology.atom(L).element.radius
        R_r = self.mdTrajectory.topology.atom(R).element.radius
        atom_contact.append(L_r + R_r + atom_atom_dist)
    incontact_atom = (atom_array < np.array(atom_contact)) 
    atom_dists = atom_array[incontact_atom]
    atom_pairs = atom_meshgrid[list(incontact_atom)]    

    incontact_vs = (vs_array <= self.max_aip_dist)  # unit: nm
    vs_dists = vs_array[incontact_vs]
    vs_pairs = vs_meshgrid[list(incontact_vs)]

    #create pandas DataFrame with Atom/AIP contacts
    AllAtomPairs = pd.DataFrame(
        {"L": atom_pairs.T[0], "R": atom_pairs.T[1], "Atom_Distance": atom_dists})
    AllAtomPairs["L_AIP"] = [self.AtomAipDict[l] if l in self.AtomAipDict.keys()
                             else np.nan for l in AllAtomPairs.L]
    AllAtomPairs["R_AIP"] = [self.AtomAipDict[r] if r in self.AtomAipDict.keys()
                             else np.nan for r in AllAtomPairs.R]
    #if atoms have no corresponding AIPs in contact 
    AllAtomPairs = AllAtomPairs.dropna()
    #split AIP cell into row with one AIP per row
    AllAtomPairs = self.split_df_row(AllAtomPairs, "L_AIP")
    AllAtomPairs = self.split_df_row(AllAtomPairs, "R_AIP")
    #add distance information to DataFrame
    AllAtomPairs["AIP_Distance"] = np.nan
    for i, row in AllAtomPairs.iterrows():
        aips = np.array(row[["L_AIP", "R_AIP"]])
        if np.any(np.all(vs_pairs == aips, axis=1)):
            i_pair = np.where(np.all(vs_pairs == aips, axis=1))
            AllAtomPairs.at[i, 'AIP_Distance'] = vs_dists[i_pair]
    AtomPairs = AllAtomPairs.dropna().reset_index(drop=True)

    #C.ar-C.ar interaction at right angle excision
    #select C.ar-C.ar interactions only
    AtomPairs["L_type"] = [self.AIP_L[int(L)].type for L in AtomPairs["L_AIP"]]
    AtomPairs["R_type"] = [self.AIP_R[int(R)].type for R in AtomPairs["R_AIP"]]
    CC_AtomPairs = AtomPairs[(AtomPairs["L_type"]== "C.ar") & (AtomPairs["R_type"]== "C.ar")]
    AtomPairs.drop("L_type",axis=1,inplace=True)
    AtomPairs.drop("R_type",axis=1,inplace=True)
    #calculate C.normal - C.normal angle
    L_xyz = [self.Atom_L[int(L)].xyz for L in CC_AtomPairs["L"]]
    R_xyz = [self.Atom_R[int(R)].xyz for R in CC_AtomPairs["R"]]
    vector_L = [-self.AIP_L[int(L)].xyz + L_xyz[i] for i, L in enumerate(CC_AtomPairs["L_AIP"])]
    vector_R = [-self.AIP_R[int(R)].xyz + R_xyz[i] for i, R in enumerate(CC_AtomPairs["R_AIP"])]
    nv_L = np.array([L/np.linalg.norm(L) for L in vector_L])
    nv_R = np.array([R/np.linalg.norm(R) for R in vector_R])
    angle = np.array([np.arccos(np.dot(L, R)) for L, R in zip(nv_L, nv_R)]) * 180 / np.pi
    CC_AtomPairs.loc[:, "angle"] = angle
    #if angle smaller than 130, drop interaction
    CC_drop = CC_AtomPairs[CC_AtomPairs.angle < 130]
    [AtomPairs.drop(i, inplace=True) for i in CC_drop.index]
    AtomPairs = AtomPairs.reset_index(drop=True)

    #add AIP fraction information 
    if len(AtomPairs > 0):
        AtomPairs["L_frac"] = [self.AIP_L[l].fraction for l in AtomPairs["L_AIP"]]
        AtomPairs["R_frac"] = [self.AIP_R[r].fraction for r in AtomPairs["R_AIP"]]
        #AIP-AIP distance takes precendence in determining contacts
        AtomPairs = AtomPairs.sort_values(["AIP_Distance", "Atom_Distance"], ignore_index=True)
        AipPairs = np.array([np.array(row) for _, row in AtomPairs.iterrows()])
        return AipPairs
    else:
        return []