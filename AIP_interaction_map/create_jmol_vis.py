import numpy as np
from AIP_interaction_map.Jmol_backbone_IP import Jmol_backbone_IP as backbone_IP
from AIP_interaction_map.Jmol_backbone import Jmol_backbone as backbone


def create_jmol_vis(self, filename, path="."):
    """Creates the visualisation using a backbone Jmol script to display all AIPs and colour 
       the interacting AIPs distinctly. Excludes H bonds"""
    self.mdTrajectory.save_pdb(f"{path}/{filename}.pdb")

    #with open(f"./ssipscore/Jmol_backbone", 'r') as f:
    #    backbone = f.readlines()

    ff = filename.split("/")[-1].split(".")[0]
    file = [s.replace('PDB_COMPLEX', f"{path}/{ff}.pdb")
            for s in backbone]
    guest_list = [int(x) for x in self.interaction_df.L_AIP]
    file = [s.replace('GUEST_LIST', " ".join(str(x)
            for x in sorted(guest_list))) for s in file]
    host_list = [int(x) for x in self.interaction_df.R_AIP]
    file = [s.replace('HOST_LIST',  " ".join(str(x)
            for x in sorted(host_list))) for s in file]

    #Adds distance information
    measure_line = '  measure ({GUEST_DIS}) ({HOST_DIS});\n'
    for i in range(len(guest_list)):
        guest = str(guest_list[i])
        host = str(host_list[i])
        insertion_index = file.index(measure_line)
        file = [s.replace('GUEST_DIS', guest) for s in file]
        file = [s.replace('HOST_DIS', host) for s in file]
        file.insert(insertion_index, measure_line)
    file.remove(measure_line)

    with open(f"{path}/{filename}_script", 'w') as fw:
        fw.writelines(file)


def create_ip_vis(self, filename, path=".", lower=0.5, strong=-5):
    """Creates the visualisation using a backbone Jmol_IP script to display interacting points (IP)
       IPs are mid-points between the interacting points, colour-coded by energy of the interaction.
       If one AIP interacts with multiple others, their energies are averaged and points superimposed."""
    
    vis = self.state.mdtraj.atom_slice(
        list(self.state.ligand_atom_indices)+list(self.state.residue_atom_indices))
    vis.save_pdb(f"{path}/{filename}.pdb")
    ff = filename.split("/")[-1].split(".")[0]
    file = [s.replace('PDB_COMPLEX', f"{path}/{ff}.pdb")
            for s in backbone_IP]

    #IP POSITIONS
    IP_positions = (np.array([self.AIP_L[i].xyz for i in self.interaction_df.L_AIP]) +
                    np.array([self.AIP_R[i].xyz for i in self.interaction_df.R_AIP]))/2
    #Will be amended for if there are many contacts to one AIP
    dG_value = np.array(self.interaction_df.ddG_sc)

    #If single AIP interacting with multiple others, they get superimposed in space
    nu = dict(self.interaction_df.L_AIP.value_counts())
    for k, v in nu.items():
        if v > 1:
            #twice interacting AIP
            tiAIP = (
                list(self.interaction_df.L_AIP[self.interaction_df.L_AIP == k].index))
            new_position = (IP_positions[tiAIP[0]] + IP_positions[tiAIP[1]]) / 2
            new_value = dG_value[tiAIP[0]] + dG_value[tiAIP[1]]
            IP_positions[tiAIP[0]] = new_position
            IP_positions[tiAIP[1]] = new_position
            dG_value[tiAIP[0]] = new_value
            dG_value[tiAIP[1]] = new_value
    nu2 = dict(self.interaction_df.R_AIP.value_counts())

    for k, v in nu2.items():
        if v > 1:
            #twice interacting AIP
            tiAIP = (
                list(self.interaction_df.R_AIP[self.interaction_df.R_AIP == k].index))
            new_position = (IP_positions[tiAIP[0]] + IP_positions[tiAIP[1]]) / 2
            new_value = dG_value[tiAIP[0]] + dG_value[tiAIP[1]]
            IP_positions[tiAIP[0]] = new_position
            IP_positions[tiAIP[1]] = new_position
            dG_value[tiAIP[0]] = new_value
            dG_value[tiAIP[1]] = new_value

    IP_line = '  data "append example"|1|IP|VS IP_POS|end "append example";show data\n'
    for i in IP_positions:
        position = '{} {} {}'.format(
            np.round(i[0]*10, 2), np.round(i[1]*10, 2), np.round(i[2]*10, 2))
        insertion_index = file.index(IP_line)
        file = [s.replace('IP_POS', position) for s in file]
        file.insert(insertion_index, IP_line)
    file.remove(IP_line)

    #IP COLOUR-CODING
    aphb = vis.n_atoms
    IP_indices = np.arange(aphb, aphb+len(self.interaction_df))
    IP_indices = np.flip(IP_indices)

    rep = IP_indices[dG_value >= lower]
    lar = IP_indices[dG_value < strong]
    sma = IP_indices[np.where((dG_value > strong) & (dG_value < -lower))]
    zero = IP_indices[np.where((dG_value > -lower) & (dG_value < lower))]

    file = [s.replace('LAR', " ".join(str(x)
                      for x in sorted(lar))) for s in file]
    file = [s.replace('SMA', " ".join(str(x)
                      for x in sorted(sma))) for s in file]
    file = [s.replace('REP', " ".join(str(x)
                      for x in sorted(rep))) for s in file]
    file = [s.replace('ZERO', " ".join(str(x)
                      for x in sorted(zero))) for s in file]

    #Also converting self.interaction_df.R indices to those in only protein_atom_indices
    atom_interspersed_dict = {}
    no_lig_atoms = len(self.Atom_L)
    for i, k in enumerate(self.Atom_R.keys()):
        atom_interspersed_dict[k] = i+no_lig_atoms

    guest_list = [int(i) for i in self.interaction_df.L]
    host_list = [atom_interspersed_dict[i] for i in self.interaction_df.R]

    #Adds distance information
    measure_line = '  draw indx ({GUEST_DIS}) ({HOST_DIS}) diameter 0.02 translucent 200 [x000000];\n'
    for i in range(len(guest_list)):
        guest = str(guest_list[i])
        host = str(IP_indices[i])
        insertion_index = file.index(measure_line)
        file = [s.replace('indx', 'obj'+str(i)) for s in file]
        file = [s.replace('GUEST_DIS', guest) for s in file]
        file = [s.replace('HOST_DIS', host) for s in file]
        file.insert(insertion_index, measure_line)
    for i in range(len(host_list)):
        guest = str(host_list[i])
        host = str(IP_indices[i])
        insertion_index = file.index(measure_line)
        file = [s.replace('indx', 'objh'+str(i)) for s in file]
        file = [s.replace('GUEST_DIS', guest) for s in file]
        file = [s.replace('HOST_DIS', host) for s in file]
        file.insert(insertion_index, measure_line)
    file.remove(measure_line)

    with open(f"{path}/{filename}_script", 'w') as fw:
        fw.writelines(file)
