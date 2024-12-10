import mdtraj as md
from AIP_interaction_map.constants import AA_LIST, LIG
import numpy as np
import logging

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

class BindingSite():
    def __init__(self, file_prot, file_pock):
        """This function takes the trajectory of the protein (which contains
        the protein with the hydrogens, and the trajectory with the pocket (which
        does not contain hydrogens) and prints out a pdb file with the binding
        pocket with the hydrogens present. out_file assigns the name of the pdb file"""
        self.__file_prot__=file_prot
        self.__file_pock__=file_pock
        self.traj_prot = md.load_pdb(self.__file_prot__)
        self.traj_pock = md.load_pdb(self.__file_pock__)
        self.traj = None

        if self.traj_prot.topology.n_chains >= 1:
            self.traj = self.pock_multi_chain()

        else:
            LOGGER.error("There are no residues in the file")

    def pock_multi_chain(self):
        """This function considers multiple chains when writing out the binding site
           The reason why this function is a bit more complicated is that if there are
           subunits of the same amino acid chains, then the names of the residues repeat 
           themselves."""

        dict_res = {}
        list_res = [] 
        res_not = []
        prot_pock_index =[]
        for i, j in enumerate(self.traj_prot.xyz[0]):
            dict_res[tuple(j)]=i
        for k in self.traj_pock.xyz[0]:
            try:
                i=dict_res[tuple(k)]
                tp = self.traj_prot.topology.atom(i).residue
                if tp not in list_res and tp.name in AA_LIST:
                    list_res.append(tp)
                    list_new = self.traj_prot.topology.select("resid == {} and chainid == {}".format(
                                tp.index, tp.chain.index))
                    prot_pock_index.extend(list_new)
                elif tp.name not in AA_LIST and tp.name not in res_not and tp.name != "HOH":
                    LOGGER.error(" Residue {} is being ignored".format(tp.name))
                    res_not.append(tp.name)
            except:
                pass
        traj = self.traj_prot.atom_slice(prot_pock_index)
        return traj

    def write_file(self, out_file):
        self.traj.save(out_file)
        LOGGER.info("written binding site")
