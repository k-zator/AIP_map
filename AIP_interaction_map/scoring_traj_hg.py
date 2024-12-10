import os
import ast
import shutil
import logging
import tempfile
import numpy as np
import mdtraj as md
from simtk import unit
from simtk.openmm import app
from simtk.openmm import VerletIntegrator
from AIP_interaction_map.create_ff import create_ff_int
from AIP_interaction_map.create_combined import create_combined
from AIP_interaction_map.constants import FF_VS_PATH, FF_VS_PATH_DUAL
logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class ScoringTraj():
    """Creates an MDTraj Trajectory and dictionaries for the AIP information necessary for the scoring.
       Takes input of guest and host PDB and AIP file paths. For proteins, the AIPs have been pre-compiled."""

    def __init__(self, guest_name, host_name, aip_guest, aip_host, protein_host=False):
        """
        Parameters
        ----------
        guest_name : string
            Path to the guest ligand files of interest
        host_name : string
            Path to the host ligand files of interest
        protein_host : bool
            if True uses precompiled xml files of residue AIPs
            to be assigned to a protein PDB instead of host
        """
        self.ligand_name = self.get_ligand_name(guest_name, protein_host)

        if aip_guest == False:
            path_to_guest_folder = os.path.dirname(host_name)
            aip_guest = f"{path_to_guest_folder}/ssip.xml"
        self.dual = self.is_dual(aip_guest)

        if protein_host == False:
            status = self.read_in_PDB(guest_name, host_name, aip_guest, aip_host)
            if status != "error":
                self.set_indices()
            else:
                self.ligand_name = status
        else:
            self.read_in_PL(guest_name, host_name, aip_guest)
    
    @staticmethod
    def is_dual(aip, dual_before=False):
        if dual_before:
            return True
        with open(aip, "r") as fp: 
            for line in fp: 
                if 'valueMax' in line: 
                    return True 
            else: 
                return False 

    @staticmethod
    def get_ligand_name(ligand_name, protein_host):
        if protein_host:
            name = "MOL"
        else:
            comp = ligand_name.split("/")[-1].split(".")[0]
            name = comp[0]+comp[-2:]
        return name

    @staticmethod
    def get_mdtraj(modeller, system):
        integrator = VerletIntegrator(1*unit.femtosecond)
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        temp_folder = tempfile.mkdtemp()+"/"
        simulation.saveState(temp_folder+"state.xml")
        mtop = md.Topology.from_openmm(modeller.topology)
        traj = md.load_xml(temp_folder + "state.xml", mtop)
        shutil.rmtree(temp_folder)
        return traj

    @staticmethod
    def get_dict_anchors(system):
        dict_anchors = {}
        for i in range(system.getNumParticles()):
            if system.isVirtualSite(i) == True:
                dict_anchors[i] = np.array([system.getVirtualSite(i).getParticle(0),
                                            system.getVirtualSite(
                                                i).getParticle(1),
                                            system.getVirtualSite(i).getParticle(2)])
        return dict_anchors

    @staticmethod
    def get_ssip_dict(system, dict_anchors, dual):
        ssip_dict = {}
        fraction_dict = {}
        isosurface_dict = {}
        owner_dict = {}
        dual_dict = {}
        for i in range(system.getForce(0).getNumParticles()):
            ssip_dict[i] = system.getForce(0).getParticleParameters(i)[0]
            fraction_dict[i] = system.getForce(0).getParticleParameters(i)[1]
            isosurface_dict[i] = system.getForce(0).getParticleParameters(i)[2]
            if dual:
                dual_dict[i] = system.getForce(0).getParticleParameters(i)[3]
            else:
                dual_dict[i] = False
            if i in dict_anchors.keys():
                owner_dict[i] = dict_anchors[i][0]
            else:
                owner_dict[i] = i
        if dual:
            return ssip_dict, fraction_dict, isosurface_dict, owner_dict, dual_dict
        else:
            return ssip_dict, fraction_dict, isosurface_dict, owner_dict
        
    def write_traj(self, out_file):
        self.out_file = out_file
        md.Trajectory.save(self.mdtraj, out_file)

    def set_indices(self):
        """This function creates all the guest and host/ ligand and residue indices necessary."""
        self.ligand_all_indices = self.mdtraj.topology.select(
            "(resname == {})".format(self.ligand_name))
        self.residue_all_indices = self.mdtraj.topology.select(
            "(resname != {})".format(self.ligand_name))
        self.ligand_vs_indices = self.mdtraj.topology.select(
            "(resname == {}) and (name =~ 'M.*')".format(self.ligand_name))
        self.residue_vs_indices = self.mdtraj.topology.select(
            "(resname != {}) and (name =~ 'M.*')".format(self.ligand_name))
        self.ligand_atom_indices = self.mdtraj.topology.select(
            "(resname == {}) and not (name =~  'M.*')".format(self.ligand_name))
        self.residue_atom_indices = self.mdtraj.topology.select(
            "(resname != {}) and not (name =~ 'M.*')".format(self.ligand_name))

    @staticmethod
    def get_sa(system):
        return system.getForce(0).getGlobalParameterDefaultValue(0), \
            system.getForce(0).getGlobalParameterDefaultValue(1), \
            system.getForce(0).getGlobalParameterDefaultValue(2), \
            system.getForce(0).getGlobalParameterDefaultValue(3)

    def read_in_PDB(self, guest, host, aip_guest, aip_host):
        """Default reading in of PDB and ssip.xml files to create combined.pdb (+vs AIPs)
           and AIP information in DataDicts.txt (saved to save time on repeated calculations).
           Combined Trajectory handed over along with dictonaries to the main scoring function."""
        path_to_guest_folder = os.path.dirname(guest)
        if not os.path.isfile(f"{path_to_guest_folder}/combined.pdb"):
            create_combined(guest, aip_guest, self.dual)
        GuestDict = open(f"{path_to_guest_folder}/DataDicts.txt", "r")
        lines = GuestDict.readlines()
        objects = []
        for line in lines:
            objects.append(ast.literal_eval(line))
            if self.dual:
                ssip_dict, isosurface_dict, fraction_dict, atom_type_dict, owner_dict, dual_dict = objects[0]
            else:
                ssip_dict, isosurface_dict, fraction_dict, atom_type_dict, owner_dict = objects[0]

        if host is not None:
            self.is_dual(aip_host, dual_before=self.dual)            
            path_to_host_folder = os.path.dirname(host)
            if not os.path.isfile(f"{path_to_host_folder}/combined.pdb"):
                create_combined(host, aip_host, self.dual)
            HostDict = open(f"{path_to_host_folder}/DataDicts.txt", "r")
            lines = HostDict.readlines()
            objects = []
            for line in lines:
                objects.append(ast.literal_eval(line))
            if self.dual:
                ssip_dict2, isosurface_dict2, fraction_dict2, atom_type_dict2, owner_dict2, dual_dict2 = objects[0]
            else:
                ssip_dict2, isosurface_dict2, fraction_dict2, atom_type_dict2, owner_dict2 = objects[0]

            guest_len = len(ssip_dict.keys())
            for i in ssip_dict2.keys():
                ssip_dict[guest_len+i] = ssip_dict2[i]
                isosurface_dict[guest_len+i] = isosurface_dict2[i]
                fraction_dict[guest_len+i] = fraction_dict2[i]
                atom_type_dict[guest_len+i] = atom_type_dict2[i]
                owner_dict[guest_len+i] = owner_dict2[i] + guest_len
                if self.dual:
                    dual_dict[guest_len+i] = dual_dict2[i]


            guest_md = md.load_pdb(f"{path_to_guest_folder}/combined.pdb")
            self.mdtraj = guest_md.stack(md.load_pdb(
                f"{path_to_host_folder}/combined.pdb"))
        self.ssip_dict = ssip_dict
        self.isosurface_dict = isosurface_dict
        self.fraction_dict = fraction_dict
        self.atom_type_dict = atom_type_dict
        self.owner_dict = owner_dict
        if self.dual:
            self.dual_dict = dual_dict
        return True

    def read_in_PL(self, guest_name, host_name, aip_guest):
        """Protein-ligand variety: it requires OpenMM to read in protein forcefield (virtual site AIPs),
           hence also creates the forcefield files for the ligands. Hands over the trajectory and all AIP
           information to MDTraj and main scoring function also."""
        guest = app.PDBFile(guest_name)
        path_to_guest_folder = os.path.dirname(guest_name)
        if not os.path.isfile(f"{path_to_guest_folder}/ligand_vs_ff.xml"):
            create_ff_int(aip_guest, path_to_guest_folder)

        guest_ff = f"{path_to_guest_folder}/ligand_vs_ff.xml"
        guest_custom = f"{path_to_guest_folder}/ligand_vs_custom.xml"
        modeller = app.Modeller(guest.topology, guest.getPositions(frame=0))
        if host_name != None:
            host = app.PDBFile(host_name)
            if self.dual:
                host_ff = FF_VS_PATH_DUAL
            else:
                host_ff = FF_VS_PATH
            modeller.add(host.topology, host.positions)
            forcefield = app.ForceField(host_ff, guest_ff, guest_custom)
        else:
            forcefield = app.ForceField(guest_ff, guest_custom)

        modeller.addExtraParticles(forcefield, ignoreExternalBonds=True)
        system = forcefield.createSystem(modeller.topology, ignoreExternalBonds=True)
        self.mdtraj = self.get_mdtraj(modeller, system)
        self.dict_anchors = self.get_dict_anchors(system)
        self.set_indices()
        if self.dual:
            self.ssip_dict, self.fraction_dict, self.isosurface_dict, self.owner_dict, self.dual_dict = \
            self.get_ssip_dict(system, self.dict_anchors, self.dual)
        else:
            self.ssip_dict, self.fraction_dict, self.isosurface_dict, self.owner_dict = \
            self.get_ssip_dict(system, self.dict_anchors, self.dual)            
        self.sa_tot, self.sa_pos, self.sa_neg, self.n_zeros = self.get_sa(system)
