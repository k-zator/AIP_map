import numpy as np
import mdtraj as md

def solvent_contact_search(self):
    """Calculates SASA for all atoms and AIPs present in the complex."""
    limit035 = 0.098
    limit070 = 0.08
    all_atoms = list(self.state.ligand_atom_indices) + list(self.state.residue_atom_indices)
    free_ligand_traj = self.mdTrajectory.atom_slice(list(self.state.ligand_atom_indices))
    free_residue_traj = self.mdTrajectory.atom_slice(list(self.state.residue_atom_indices))    
    traj = self.mdTrajectory.atom_slice(all_atoms)

    try:
        sa_free, _ = md.shrake_rupley(
            free_ligand_traj, probe_radius=0.0, n_sphere_points=1000, get_mapping=False)
        sasa35f, _ = md.shrake_rupley(
            free_ligand_traj, probe_radius=0.035, n_sphere_points=1000, get_mapping=False)
        sasa70f, _ = md.shrake_rupley(
            free_ligand_traj, probe_radius=0.070, n_sphere_points=1000, get_mapping=False)
        sasa35Rf, _ = md.shrake_rupley(
            free_residue_traj, probe_radius=0.035, n_sphere_points=1000, get_mapping=False)
        sasa70Rf, _ = md.shrake_rupley(
            free_residue_traj, probe_radius=0.070, n_sphere_points=1000, get_mapping=False)    
        sasa35 = get_sasa_per_AIP(self, traj, 0.035)
        sasa70 = get_sasa_per_AIP(self, traj, 0.070)
    except:
        sa_free = md.shrake_rupley(
            free_ligand_traj, probe_radius=0.0, n_sphere_points=1000, get_mapping=False)
        sasa35f = md.shrake_rupley(
            free_ligand_traj, probe_radius=0.035, n_sphere_points=1000, get_mapping=False)
        sasa70f = md.shrake_rupley(
            free_ligand_traj, probe_radius=0.070, n_sphere_points=1000, get_mapping=False)
        sasa35Rf = md.shrake_rupley(
            free_residue_traj, probe_radius=0.035, n_sphere_points=1000, get_mapping=False)
        sasa70Rf = md.shrake_rupley(
            free_residue_traj, probe_radius=0.070, n_sphere_points=1000, get_mapping=False)
        sasa35 = md.shrake_rupley(
            traj, probe_radius=0.035, n_sphere_points=1000, get_mapping=False)
        sasa70 = md.shrake_rupley(
            traj, probe_radius=0.070, n_sphere_points=1000, get_mapping=False)
    # rid of AIPs with too large a SASA values per AIP:
    self.solvent_contact = []
    for i, a in enumerate(all_atoms):
        if a in self.Atom_L.keys():
            self.Atom_L[a].set_sasa35(sasa35[0][i])
            self.Atom_L[a].set_sasa70(sasa70[0][i])
            self.Atom_L[a].set_sasa35f(sasa35f[0][i])
            self.Atom_L[a].set_sasa70f(sasa70f[0][i])
            self.Atom_L[a].set_sa(sa_free[0][i])
        else:
            self.Atom_R[a].set_sasa35(sasa35[0][i])
            self.Atom_R[a].set_sasa70(sasa70[0][i])
            self.Atom_R[a].set_sasa35f(sasa35Rf[0][i - len(self.state.ligand_atom_indices)])
            self.Atom_R[a].set_sasa70f(sasa70Rf[0][i - len(self.state.ligand_atom_indices)])           
        if sasa35[0][i] > limit035 and sasa70[0][i] > limit070:
            if a < len(self.state.ligand_atom_indices):
                if self.Atom_L[a].polar == True: #due to 4k18
                    self.solvent_contact.append(a)
            else:
                if self.Atom_R[a].polar == True: #due to 4k18
                    self.solvent_contact.append(a)

def exclude_solvent_contact(self):
    """Based off limit criteria, decides if the AIPs are exposed to solvent. If so, they get deleted
        from the AIP dictionaries and therefore completely excluded from further consideration"""
    limit035 = 0.098
    limit070 = 0.034 #due to 4twp
    limit035s = 0.04
    for a in self.AtomAipDict.keys():
        if a < len(self.state.ligand_atom_indices):
            if self.Atom_L[a].sasa35 > limit035*2 and self.Atom_L[a].sasa70 > limit070*2 \
                and self.Atom_L[a].type in ["Cl", "Br", "I", "S.3", "S.2.phene"]:
                self.solvent_contact.append(a)
            elif self.Atom_L[a].sasa35 > limit035 and self.Atom_L[a].sasa70 > limit070*2 \
                and self.Atom_L[a].type in ["F"]:
                self.solvent_contact.append(a)
            elif self.Atom_L[a].sasa35 > limit035 and self.Atom_L[a].sasa70 > limit070 \
                and self.Atom_L[a].type not in ["F", "Cl", "Br", "I", "S.3", "S.2.phene"]:
                self.solvent_contact.append(a)
            try:
                for ai in self.AtomAipDict[a]:
                    aip_sasa = self.AIP_L[ai].sasa_b35
                    if aip_sasa > limit035s and self.Atom_L[a].polar == False and self.AIP_L[ai].type[0] in ["N", "C", "O"]:
                        self.solvent_contact.append(ai)
                    elif aip_sasa > limit035s*2 and self.AIP_L[ai].type in ["Cl", "Br", "I", "S.3", "S.2.phene"]:
                        self.solvent_contact.append(a)
            except:
                pass
        else:
            if self.Atom_R[a].sasa35 > limit035*3 and self.Atom_R[a].sasa70 > limit070*3 \
                and self.Atom_R[a].type in ["S.3", "S.2.phene"]:
                self.solvent_contact.append(a)
            elif self.Atom_R[a].sasa35 > limit035 and self.Atom_R[a].sasa70 > limit070 \
                and self.Atom_R[a].type not in ["S.3", "S.2.phene"]:
                self.solvent_contact.append(a)
            try:
                for ai in self.AtomAipDict[a]:
                    aip_sasa = self.AIP_R[ai].sasa_b35
                    if aip_sasa > limit035s and self.Atom_R[a].polar == False and self.AIP_R[ai].type[0] in ["N", "C", "O"]:
                        self.solvent_contact.append(ai)
                    elif aip_sasa > limit035s*2 and self.AIP_R[ai].type in ["S.3", "S.2.phene"]:
                        self.solvent_contact.append(a)
            except:
                pass

def get_sasa_per_AIP(self, traj, probe_radius, bound=True):
    """wrapper for shrake_rupley to also quantify sasa changed per AIP"""
    Atom = self.Atom_L
    AIP = self.AIP_L
    no_atoms = traj.n_atoms

    y, a = md.shrake_rupley(traj, probe_radius=probe_radius, n_sphere_points=1000, get_mapping=False)
    b = a.reshape(no_atoms,1000,3)
        
    for atom_i in Atom.keys():   
        aip_ind = np.array([i.index for i in AIP.values() if (i.atom == atom_i)])
        if len(aip_ind) == 0:
            continue
 
        area_atom_i = y[0][atom_i]

        if area_atom_i > 0:
            aip_xyz = np.array([i.xyz for i in AIP.values() if (i.atom == atom_i)])
            data = b[atom_i]
            data = data[~np.all(data == 0, axis=1)]
            JJ = np.array([np.linalg.norm(data-a, axis=1) for a in aip_xyz]).T
            jj = aip_ind[JJ.argmin(axis=1)]
            no_jj = len(jj) / area_atom_i
            for k in aip_ind:
                if bound == True and probe_radius == 0.035:
                    AIP[k].set_sasa_b35(np.count_nonzero(jj == k) / no_jj)
                elif bound == True and probe_radius == 0.070:
                    AIP[k].set_sasa_b70(np.count_nonzero(jj == k) / no_jj)
                else:
                    AIP[k].set_sasa_f(np.count_nonzero(jj == k) / no_jj)
        else:       
            for k in aip_ind:
                if bound == True and probe_radius == 0.035:
                    AIP[k].set_sasa_b35(0.0)
                elif bound == True and probe_radius == 0.070:
                    AIP[k].set_sasa_b70(0.0)
                else:
                    AIP[k].set_sasa_f(0.0)

    for jind, atom_j in enumerate(self.Atom_R.keys()):   
        aip_ind = np.array([j.index for j in self.AIP_R.values() if (j.atom == atom_j)])
        if len(aip_ind) == 0:
            continue
            
        atom_jj = atom_j
        atom_j = jind + len(Atom)
        area_atom_j = y[0][atom_j]

        if area_atom_j > 0:
            aip_xyz = np.array([j.xyz for j in self.AIP_R.values() if (j.atom == atom_jj)])
            data = b[atom_j]
            data = data[~np.all(data == 0, axis=1)]
            JJ = np.array([np.linalg.norm(data-a, axis=1) for a in aip_xyz]).T
            jj = aip_ind[JJ.argmin(axis=1)]
            no_jj = len(jj) / area_atom_j
            for k in aip_ind:
                if bound == True and probe_radius == 0.035:
                    self.AIP_R[k].set_sasa_b35(np.count_nonzero(jj == k) / no_jj)
                elif bound == True and probe_radius == 0.070:
                    self.AIP_R[k].set_sasa_b70(np.count_nonzero(jj == k) / no_jj)
                else:
                    self.AIP_R[k].set_sasa_f(np.count_nonzero(jj == k) / no_jj)
        else:       
            for k in aip_ind:
                if bound == True and probe_radius == 0.035:
                    self.AIP_R[k].set_sasa_b35(0.0)
                elif bound == True and probe_radius == 0.070:
                    self.AIP_R[k].set_sasa_b70(0.0)
                else:
                    self.AIP_R[k].set_sasa_f(0.0)
    return y