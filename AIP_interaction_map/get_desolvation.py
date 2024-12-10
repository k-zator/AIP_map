import numpy as np
import mdtraj as md
from AIP_interaction_map.get_sasa import get_sasa_per_AIP

polynomials = {
"water_alpha": [0.838237, -0.51907000, -7.714317E-02, -8.1380530E-01, 3.524312E-01, -7.0570040E-02, 7.6254820E-03, -4.2989670E-04, 9.9333680E-06],
"water_beta": [0.7751790, 1.258973006E-02, -3.7236050E-01, 4.1098470E-02, 1.7334470E-02, 2.097013000E-03, 1.2569490E-04, 3.79972E-06, 4.62349E-08],
"hexadecane_alpha": [-1.257667, -1.343249E-02, -1.0086660E-03, 2.992371E-06, 9.5324600E-05, -2.9965860E-05, 4.0427420E-06, -2.6054980E-07, -6.573138E-09],
"hexadecane_beta": [-1.256501, 1.30147200E-01, -3.56959E-03, 2.6307740E-04, 3.59959700E-05, 2.771214000E-06, 1.3473240E-07, 3.622555E-09, 4.0709470E-11],
"chloroform_alpha": [-1.220126, -1.662607E-01, -9.94286E-03, -4.672309E-04, -2.678255E-06, 1.475957E-05, -2.411077E-06, 1.74012E-07, -4.822481E-09],
"chloroform_beta": [-1.187073, 4.507678E-01, 1.42246E-01, 9.4427E-02, 1.417198E-02, 1.072633E-03, 4.489395E-05, 9.877984E-07, 8.891901E-09]}
theta = {
"theta_water": 0.73848929,
"theta_hexadecane": 0.41933333,
"theta_chloroform": 0.45466667}

def getDesolvation(self):
    """Calculate SASA changes for all atoms and AIPs and figure out if there is a discrepancy
       between paired and desolvated AIPs. If so, the latter ones are allowed to pair again 
       with a larger distance criterion (+ extension upon desolvation). If they are not found,
       perhaps they are just transferred to a non-polar medium, hence (phase) transfer energies
       are calculated."""
    
    self.max_aip_dist = self.max_aip_dist + self.ext_upon_desolv

    self.sasa_free, self.frac_desolv, desolv_energy_1st, desolv_energy_2nd, transfer_energy_1st = \
        get_desolvated_AIPs(self)
    
    self.desolv_energy = sum(desolv_energy_1st) + sum(desolv_energy_2nd)
    self.transfer_energy = sum(transfer_energy_1st)

    all_atoms = np.array(list(self.state.ligand_atom_indices)+list(self.state.residue_atom_indices))
    desolv_nonint_atoms = all_atoms[transfer_energy_1st != 0]
    for _, row in self.interaction_df.iterrows():
        if row.L_frac == 1 and row.R_frac == 0.5 and sum(self.interaction_df.L_AIP == row.L_AIP) == 1:
            self.AIP_L[row.L_AIP].fraction = 0.5
            desolv_nonint_atoms = np.append(desolv_nonint_atoms, int(row.L))
        elif row.R_frac == 1 and row.L_frac == 0.5 and sum(self.interaction_df.R_AIP == row.R_AIP) == 1:
            self.AIP_R[row.R_AIP].fraction = 0.5
            desolv_nonint_atoms = np.append(desolv_nonint_atoms, int(row.R))

    L_Atom = [a for a in desolv_nonint_atoms if a in self.Atom_L]
    R_Atom = [a for a in desolv_nonint_atoms if a in self.Atom_R]
    L_AIP = [k for l in L_Atom for k in self.AtomAipDict[l]]
    R_AIP = [k for l in R_Atom for k in self.AtomAipDict[l]]

    #second pass of the branch algorithm, now just for desolvated sites
    self.non_polar_add = self.getNonPolar(True, L_Atom, L_AIP, R_Atom, R_AIP)
    

def get_free_energy_of_solvation(aip, solvent):
    """Calculate free energy of solvation for a given AIP using the AIP value and solvent"""

    if float(aip) >= 0:
        sign = "alpha"
    else:
        sign = "beta"

    #choice of solvent
    if solvent not in ["water", "chloroform", "hexadecane"]:
        return "Please pick an available solvent"
    else:
        poly = polynomials[f"{solvent}_{sign}"]
    poly_indices = np.arange(9)
    dG = sum(np.power(aip, poly_indices) * poly)
    return dG

def get_free_energy_of_phase_transfer(aip, s_from, s_to):
    """Calculate free energy of phase transfer for a given AIP using the AIP value and 
       solvent from and solvent to which the AIP is transferred."""

    dG_solv_from = get_free_energy_of_solvation(aip, s_from)
    dG_solv_to = get_free_energy_of_solvation(aip, s_to)
    return dG_solv_to - dG_solv_from

def get_PrimInterAipDict(self):
    """Primarily Interacting AIP dictionary:
       dictionary ordered in sequence of likely interacting and hence to be desolvated AIPs
       deprecated use soon?"""
    
    L_AIP = [l for l in list(self.AIP_L.keys())]
    R_AIP = [r for r in list(self.AIP_R.keys())]
    # obtaining meshgrid for determining pairings
    vs_meshgrid = np.array(np.meshgrid(L_AIP, R_AIP)).T.reshape(-1, 2)
    vs_dist = md.compute_distances(self.mdTrajectory, vs_meshgrid, periodic=False)
    vs_array = vs_dist[0]

    PrimInterAipDict = {}
    for k, aips in self.AtomAipDict.items():
        if len(aips) == 1:
            #for single atom-AIP, just that one
            PrimInterAipDict[k] = aips
        else:
            #for multiple AIPs, find min distance for cross-AIPs for each
            if aips[0] in self.AIP_L.keys():
                distances = [vs_array[vs_meshgrid.T[0] == i].min() for i in aips]
            #0 index only works for L, then need to switch to R, so index 1 and look for the same value set
            else:
                distances = [vs_array[vs_meshgrid.T[1] == i].min() for i in aips]
            #order them from closest to furthest away
            PrimInterAipDict[k] = [val for (_, val) in sorted(zip(distances, aips), 
                                   key=lambda x: x[0], reverse=False)]
    return PrimInterAipDict #remember they'll be ordered lists

def get_desolvated_AIPs(self):
    """Find AIPs whose percentage change from free species to bound ones is at least 50% (default). 
       SASA is calculated using the shrake-rupley algorithm at a given probe radius (default 0.2)."""
    
    #find trajectories of free and bound species 
    all_atoms = list(self.state.ligand_atom_indices)+list(self.state.residue_atom_indices)
    all_area = [4*np.pi*(self.mdTrajectory.topology.atom(mol).element.radius + self.probe_radius)**2 for mol in all_atoms]
    traj = self.mdTrajectory.atom_slice(all_atoms)
    traj_lig = self.mdTrajectory.atom_slice(self.state.ligand_atom_indices)
    traj_res = self.mdTrajectory.atom_slice(self.state.residue_atom_indices)

    #calculate sasa_per_AIP
    sasa_bound = get_sasa_per_AIP(self, traj, self.probe_radius, bound=True, L=True)[0] / all_area
    _ = get_sasa_per_AIP(self, traj, self.probe_radius, bound=True, L=False)
    sasa_lig = get_sasa_per_AIP(self, traj_lig, self.probe_radius, bound=False, L=True)[0]
    sasa_res = get_sasa_per_AIP(self, traj_res, self.probe_radius, bound=False, L="Rinit")[0]
    sasa_free = np.concatenate((sasa_lig, sasa_res)) / all_area

    #and calculate fraction of SASA change
    sasa_diff = (sasa_free - sasa_bound) 
    frac_desolv = np.divide(sasa_diff, sasa_free, out=np.zeros_like(sasa_diff), where=sasa_free!=0)

    #find corresponding AIPs which are desolvated as the effect and energy of desolvation and their phase transfer
    ####NEEDS UPDATING
    PrimInterAipDict = get_PrimInterAipDict(self)
    AIP_dict = self.AIP_L.copy()
    AIP_dict.update(self.AIP_R)
    desolv_energy_1st = [-get_free_energy_of_solvation(AIP_dict[PrimInterAipDict[a][0]].value, self.solvent) 
                         if frac_desolv[i] > self.frac_to_desolv and a in PrimInterAipDict.keys() 
                         else 0 for i, a in enumerate(all_atoms)]
    desolv_energy_2nd = [-get_free_energy_of_solvation(AIP_dict[PrimInterAipDict[a][1]].value, self.solvent) 
                         if frac_desolv[i] > 0.75 and a in PrimInterAipDict.keys() and len(PrimInterAipDict[a]) > 1 
                         else 0 for i, a in enumerate(all_atoms)]
    
    int_aips =  list(self.interaction_df.R_AIP)+list(self.interaction_df.L_AIP)
    transfer_energy_1st = [get_free_energy_of_phase_transfer(AIP_dict[PrimInterAipDict[a][0]].value, self.solvent, "hexadecane") 
                         if frac_desolv[i] > self.frac_to_desolv and a in PrimInterAipDict.keys() 
                                                            and PrimInterAipDict[a][0] not in int_aips
                         else 0 for i, a in enumerate(all_atoms)]
    
    #and scaling by their respective AIP fractions
    fractions = np.array([AIP_dict[PrimInterAipDict[a][0]].fraction if a in PrimInterAipDict.keys() else 0 for a in all_atoms])
    desolv_energy_1st = np.array(desolv_energy_1st) * fractions
    desolv_energy_2nd = np.array(desolv_energy_2nd) * fractions
    transfer_energy_1st = np.array(transfer_energy_1st) * fractions

    return sasa_free, frac_desolv, desolv_energy_1st, desolv_energy_2nd, transfer_energy_1st 