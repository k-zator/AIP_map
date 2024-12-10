#/usr/env/pytho
import copy
import mdtraj as md
import numpy as np
import logging
import re
from AIP_interaction_map.constants import AA_LIST
logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

def copy_pdb (ligand_file, protein_file):
    """By inserting the pdb code of the desired structure, 
       this function returns a deep copy of the pdb."""
    traj1 = md.load(ligand_file)
    traj2 = md.load(protein_file)
    traj = traj1.stack(traj2) 
    return copy.deepcopy(traj)

def find_ligand_atoms (traj):
    """Finds the indices of all the atoms of the ligand present in the pdb."""
    ligand_indices = []
    for atom in traj.topology.residue(0).atoms:
        ligand_indices.append(atom.index)
    return ligand_indices

def find_protein_atoms (traj):
    protein_indices = []
    for atom in traj.topology.atoms:
        if atom.residue.name.upper() in AA_LIST :
             protein_indices.append(atom.index)
    return protein_indices

def find_neighbours (traj, cutoff):
    """Finds the residue sequences of the residues that have at least
       one atom within a cutoff distance of the ligand."""
    protein_res_less_than_cutoff = set ()
    table, bond = traj.topology.to_dataframe()
    ligand_atom_index = find_ligand_atoms (traj)
    list_seq = []
    dict_atoms = {}
    for i in ligand_atom_index:
        a = np.array ([int (i)])
        for index in md.compute_neighbors(traj, cutoff, a)[0]:
            res = traj.topology.atom(index).residue
            if res.name in AA_LIST:
                if res.index not in list_seq:
                    list_seq.append(res.index)
    return list_seq

def indices_of_all_the_neigh_res_atoms(traj, cutoff):
    """Finds the indices of all the atoms of the amino acid residues that are neighbouring the ligand."""
    indices = []
    protein_res_less_than_cutoff = find_neighbours (traj, cutoff)
    for seq in protein_res_less_than_cutoff:
        res_id= traj.topology.select('resid {}'.format(seq))
        indices.extend(res_id)
    indices.sort()
    return indices

def protein_neighbours_pdb(ligand_file, protein_file, output_file, cutoff):
    """Creates a pdb file which has the protein atoms of the neighbouring residues only/"""
    traj = copy_pdb (ligand_file, protein_file)
    indices = indices_of_all_the_neigh_res_atoms(traj, cutoff)
    traj.atom_slice(indices, inplace=True)
    traj.save(output_file)
    LOGGER.info("created binding site")
    return traj
