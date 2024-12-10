import os
import numpy as np
from AIP_interaction_map.ssip_molecule import SsipMolecule


def create_combined(PDB_file, aip_path=False, dual=False, set_custom_anchors=False):
    """Creates cpmbined PDB file with the AIPs as virtual sites. Other AIP
       information is saved as a DataDicts.txt file in the same directory."""
    
    #read in AIP information
    path_to_folder = os.path.dirname(PDB_file)
    if aip_path == False:
        mol = SsipMolecule(path_to_folder+"/ssip.xml", set_custom_anchors)
    else:
        mol = SsipMolecule(aip_path, set_custom_anchors)
    #read in PDB, it contains connect record
    with open(PDB_file) as f:
        pdb_atoms = f.readlines()
    f.close()
    
    atom_pos_dict = {}
    for line in pdb_atoms:
        if "ATOM" in line:
            a = line.split()
            atom_pos_dict[int(a[1])] = np.array(
                    [float(a[6]), float(a[7]), float(a[8])])
    atom_pos_dict
    no_atoms = len(atom_pos_dict)
    comp_name = pdb_atoms[10].split()[3]

    ssip_dict = {}
    fraction_dict = {}
    isosurface_dict = {}
    atom_type_dict = {}
    owner_dict = {}
    if dual:
        dual_dict = {}

    for k, a in mol.atom_dict.items():
        #because dict need to start from 0 but the names of atoms do not
        i = int(k[1:]) - 1
        ssip_dict[i] = 0
        fraction_dict[i] = 0
        isosurface_dict[i] = 0
        atom_type_dict[i] = a.type
        owner_dict[i] = i
        if dual:
            dual_dict[i] = False

    #find positions of AIPs
    vs_pos_dict = {}
    for i, s in enumerate(mol.get_vs_list(set_custom_anchors)):
        rel_atoms = [s.anchor1.name[1:], s.anchor2.name[1:], s.anchor3.name[1:]]
        weights = s.weights
        r1 = atom_pos_dict[int(rel_atoms[1])] - atom_pos_dict[int(rel_atoms[0])]
        r2 = atom_pos_dict[int(rel_atoms[2])] - atom_pos_dict[int(rel_atoms[0])]
        rcross = np.cross(r1, r2)
        vs_pos_dict[i] = atom_pos_dict[int(rel_atoms[0])] \
                + r1 * float(weights[0]) \
                + r2 * float(weights[1]) \
                + rcross * float(weights[2])/10
        
        ind = no_atoms + i
        ssip_dict[ind] = s.value
        fraction_dict[ind] = s.fraction
        isosurface_dict[ind] = s.isosurface
        atom_type_dict[ind] = s.type
        owner_dict[ind] = int(s.neigh[1:]) - 1
        if dual:
            dual_dict[ind] = s.dual

    #add AIPs as virtual sites
    vs_indices = list(vs_pos_dict.keys())
    vs_indices.reverse()
    
    for k in vs_indices:
        i = k + no_atoms + 1
        if i < 100:
            space1 = " "
        else:
            space1 = ""
        if int(k) < 10:
            space2 = "   "
        elif int(k) >= 100:
            space2 = " "
        else:
            space2 = "  "
        x = round(vs_pos_dict[k][0], 3)
        if float(x) >= 0 and float(x) < 10:
            space3 = "  "
        elif float(x) <= -10:
            space3 = ""
        else:
            space3 = " "
        y = round(vs_pos_dict[k][1], 3)
        if float(y) >= 0 and float(y) < 10:
            space4 = "   "
        elif float(y) <= -10:
            space4 = " "
        else:
            space4 = "  "
        z = round(vs_pos_dict[k][2], 3)
        if float(z) >= 0 and float(z) < 10:
            space5 = "   "
        elif float(z) <= -10:
            space5 = " "
        else:
            space5 = "  "

        PDB_line = "ATOM    "+space1+str(i)+space2+"M"+str(k)+" "+comp_name+" A   1     "+space3 + \
            '%.3f' % x+space4+'%.3f' % y+space5+'%.3f' % z+"  1.00  0.00          VS\n"
        pdb_atoms.insert(no_atoms+2, PDB_line)

    #save both files
    with open(path_to_folder+"/"+"combined.pdb", 'w') as fw:
        fw.writelines(pdb_atoms)

    if dual:
        data_dicts = [ssip_dict, isosurface_dict, fraction_dict, atom_type_dict, owner_dict, dual_dict]
    else:
        data_dicts = [ssip_dict, isosurface_dict, fraction_dict, atom_type_dict, owner_dict]
    outputFile = open(path_to_folder+"/"+"DataDicts.txt", "w")
    outputFile.write(str(data_dicts))
    outputFile.flush()
    outputFile.close()
