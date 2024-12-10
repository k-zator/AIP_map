import networkx as nx
from networkx.algorithms import shortest_path_length
import numpy as np
from lxml import etree
import logging
import copy
import time
import re
from AIP_interaction_map.constants import CML_NS, SSIP_NS
import networkx.algorithms.isomorphism as iso
from mendeleev import element

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class SsipMolecule:
    """
    A class used to represent information regarding atoms and ssips described in the ssip.xml files
    Attributes
    ----------
    atom_dict : dict
        A dictionary that describes the atoms present in the molecule. 
        The key is the cml atom name and the value is an SsipAtom object.
    ssip_list : list
       A list of SsipSsip objects describing all the ssips present in the molecule.
    bond_list : list
       A list of SsipBond objects describing the connectivity of the molecule.
    ssip_network : nx.Graph
       A nx.Graph representation of the molecule. 
       The NODES are the cml atom names of the molecule, and these have an attribute elementType which represents the element of the atom
       The EDGES represent which atoms are connected and these have an attribute bondOrder which represents bond order of the bond."""

    def __init__(self, ssip_file, remove_zeros=False,
                 averaged_atom_ssip=False, set_custom_anchors=False):
        self.__remove_zeros__ = remove_zeros
        self.tree = copy.deepcopy(etree.parse(ssip_file))
        self.n_removed_zeros = 0
        self.ligand_tree = None
        self.dual = False
        self.resname = "MOL"
        self.ssip_list = self.__get_ssip_list__()
        self.atom_dict = self.__get_atom_dict__()
        self.bond_list = self.__get_bond_list__()
        self.network = self.__get_network__()
        self.surface = self.SsipSurface(self.__get_surface_path__()[0])
        #self.dict_average = self.__get_dict_average__()
        self.atom_tree = self.set_atom_tree(
            averaged_atom_ssip=averaged_atom_ssip)
        self.get_ssip_name_dict()
        try:
            self.vs_tree = self.set_vs_tree(set_custom_anchors)
        except ValueError:
            ValueError("Could not find a VS force field")

    @staticmethod
    def write_ligand_xml(tree, out_file):
        with open(out_file, "wb") as writefile:
            writefile.write(etree.tostring(tree, pretty_print=True))
            LOGGER.info("written forcefield")
        return

    def write_custom_xml(self, tree, out_file, protein_xml=None):
        custom_tree = self.get_custom_xml(tree, protein_xml=protein_xml)
        with open(out_file, "wb") as writefile:
            writefile.write(etree.tostring(custom_tree, pretty_print=True))
        LOGGER.info("written custom")
        return

    def set_vs_tree(self, set_custom_anchors):
        vs_tree = copy.deepcopy(self.atom_tree)
        root = vs_tree.xpath("//AtomTypes")
        new_element = etree.SubElement(root[0], "Type")
        new_element.attrib["class"] = "ssip"
        new_element.attrib["name"] = "ssip"
        new_element.attrib["mass"] = "0.0"
        root = vs_tree.xpath(
            "//Residues/Residue[@name='{}']/Atom".format(self.resname))
        root_tag = root[0].tag
        root_parent = root[0].getparent()
        try:
            ssip_list = self.get_vs_list(set_custom_anchors)
            for i, ssip in enumerate(ssip_list):
                new_element = etree.SubElement(root_parent, root_tag)
                new_element.attrib["ssip_charge"] = "{}".format(ssip.value)
                name = "M{}_{}".format(i, ssip.type)
                new_element.attrib["fraction"] = str(ssip.fraction)
                new_element.attrib["isosurface"] = str(ssip.isosurface)
                if self.dual:
                    new_element.attrib["dual"] = str(ssip.dual)
                new_element.attrib["name"] = name
                new_element.attrib["type"] = "ssip"
                vs_element = etree.SubElement(root_parent, "VirtualSite")
                vs_element.attrib["type"] = "outOfPlane"
                vs_element.attrib["siteName"] = name
                vs_element.attrib["atomName1"] = ssip.anchor1.name
                vs_element.attrib["atomName2"] = ssip.anchor2.name
                vs_element.attrib["atomName3"] = ssip.anchor3.name
                vs_element.attrib["weight12"] = str(round(ssip.weights[0], 1))
                vs_element.attrib["weight13"] = str(round(ssip.weights[1], 1))
                vs_element.attrib["weightCross"] = str(round(ssip.weights[2]))
            return vs_tree
        except:
            raise ValueError

    def set_atom_tree(self, averaged_atom_ssip=True):
        atom_tree = copy.deepcopy(self.get_xml_tree())
        atom_dict = self.get_atom_value_dict(average=averaged_atom_ssip)
        for k, v in atom_dict.items():
            root = atom_tree.xpath(
                '//Residues/Residue[@name="{}"]/Atom[@name="{}"]'.format(self.resname, k))
            root[0].attrib["ssip_charge"] = str(v)
            root[0].attrib["fraction"] = str(0)
            root[0].attrib["isosurface"] = str(0)
            if self.dual:
                root[0].attrib["dual"] = str(0)
        return atom_tree

    def get_xml_tree(self):
        tree = etree.Element("ForceField")
        atype = etree.SubElement(tree, "AtomTypes")
        for k, v in self.atom_dict.items():
            types = etree.SubElement(atype, "Type")
            types.attrib["class"] = "{}".format(k)
            types.attrib["name"] = "ligand-{}".format(k)
            types.attrib["element"] = "{}".format(v.elem)
            types.attrib["mass"] = "{}".format(element(v.elem).atomic_weight)
        residues = etree.SubElement(tree, "Residues")
        residue = etree.SubElement(residues, "Residue")
        residue.attrib["name"] = self.resname
        for k, v in self.atom_dict.items():
            atom_elem = etree.SubElement(residue, "Atom")
            atom_elem.attrib["name"] = k
            atom_elem.attrib["type"] = "ligand-{}".format(k)
        for bond in self.bond_list:
            bond_elem = etree.SubElement(residue, "Bond")
            bond_elem.attrib["atomName1"] = bond.pair[0]
            bond_elem.attrib["atomName2"] = bond.pair[1]
        return tree

    def get_custom_xml(self, ligand_tree, protein_xml=None):
        """This function takes in the ligand force field file and the protein
           force field file and outputs a CustomNonBonded force field which
           is necessary to assign the ssip value to each atom present in the
           protein ligand system of interest."""

        tree_new = etree.Element("ForceField")
        custom_elem = etree.SubElement(tree_new, "CustomNonbondedForce")
        custom_elem.attrib["energy"] = "n_zeros*SA*SAneg*SApos*ssip_charge1*ssip_charge2"
        custom_elem.attrib["bondCutoff"] = "1"
        perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
        perpar_elem.attrib["name"] = "ssip_charge"
        perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
        perpar_elem.attrib["name"] = "fraction"
        perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
        perpar_elem.attrib["name"] = "isosurface"
        if self.dual:
            perpar_elem = etree.SubElement(custom_elem, "PerParticleParameter")
            perpar_elem.attrib["name"] = "dual"        
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "SA"
        global_elem.attrib["defaultValue"] = self.surface.tot
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "SApos"
        global_elem.attrib["defaultValue"] = self.surface.pos
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "SAneg"
        global_elem.attrib["defaultValue"] = self.surface.neg
        global_elem = etree.SubElement(custom_elem, "GlobalParameter")
        global_elem.attrib["name"] = "n_zeros"
        global_elem.attrib["defaultValue"] = str(self.n_removed_zeros)
        attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
        attfrom_elem.attrib["name"] = "ssip_charge"
        attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
        attfrom_elem.attrib["name"] = "fraction"
        attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
        attfrom_elem.attrib["name"] = "isosurface"
        if self.dual:
            attfrom_elem = etree.SubElement(custom_elem, "UseAttributeFromResidue")
            attfrom_elem.attrib["name"] = "dual"
        for i in ligand_tree.xpath("//AtomTypes/Type"):
            atom_elem = etree.SubElement(custom_elem, "Atom")
            atom_elem.attrib["type"] = i.attrib["name"]
        if protein_xml != None:
            tree_prot = copy.deepcopy(etree.parse(protein_xml))
            for i in tree_prot.xpath("//AtomTypes/Type"):
                atom_elem = etree.SubElement(custom_elem, "Atom")
                atom_elem.attrib["type"] = i.attrib["name"]
        return tree_new

    def get_atom_value_dict(self, average=True):
        if average == True:
            atom_value_dict = {}
            for k, v in self.atom_dict.items():
                atom_value_dict[k] = v.get_atom_ssip(
                    self.__get_dict_average__())
            return atom_value_dict
        elif average == False:
            atom_value_dict = {}
            for k, v in self.atom_dict.items():
                atom_value_dict[k] = v.get_atom_ssip_no_average(self.ssip_list)
            return atom_value_dict

    def get_ssip_name_dict(self, average=False):
        ssip_name_dict = {}
        for k, v in self.atom_dict.items():
            atom_index = int(re.sub('\D', '', k))
            ssip_name_dict[atom_index] = v.neigh_ssip_name
        self.ssip_name_dict = ssip_name_dict

    def get_vs_list(self, set_custom_anchors):
        virtual_site_ssip_list = []
        ssip_none = []
        for i in self.ssip_list:
            ssip = i.get_anchors_weights(self.network, self.atom_dict, set_custom_anchors)
            if ssip.anchor1 != None:
                virtual_site_ssip_list.append(ssip)
            else:
                ssip_none.append(ssip)
        if ssip_none == []:
            return virtual_site_ssip_list
        else:
            #LOGGER.error("It was not possible to find the virtual site descriptions of all the SSIPs")
            raise ValueError(
                "Not possible to find the virtual site descriptions of all the SSIPs")

    def __get_ssip_list__(self):
        ssip_list = []
        self.n_removed_zeros = 0
        if self.__remove_zeros__ == True:
            for i in self.tree.xpath("//ssip:SSIP", namespaces={"ssip": SSIP_NS, "cml": CML_NS}):
                ssip = self.SsipSsip(i)
                if ssip.value != 0.0:
                    ssip_list.append(ssip)
                else:
                    self.n_removed_zeros += 1
        elif self.__remove_zeros__ == False:
            for i in self.tree.xpath("//ssip:SSIP", namespaces={"ssip": SSIP_NS, "cml": CML_NS}):
                ssip = self.SsipSsip(i)
                ssip_list.append(ssip)
        else:
            LOGGER.error("remove_zeros is of type {}, but needs to be boolean".format(
                str(type(self.__remove_zeros__))))
        return ssip_list

    def __get_atom_dict__(self):
        atom_dict = {}
        for i in self.tree.xpath("//cml:atom", namespaces={"ssip": SSIP_NS, "cml": CML_NS}):
            atom = self.SsipAtom(i, self.ssip_list)
            atom_dict[atom.name] = atom
        return atom_dict

    def __get_bond_list__(self):
        bond_list = []
        for i in self.tree.xpath("//cml:bondArray/cml:bond", namespaces={"ssip": SSIP_NS, "cml": CML_NS}):
            bond_list.append(self.SsipBond(i))
        return bond_list

    def __get_network__(self):
        """A function to obtain the graph representation of the molecule of interest
        """
        Atom = nx.Graph()
        for i in self.atom_dict.values():
            Atom.add_node(i.name, elementType=i.elem)
        for i in self.bond_list:
            Atom.add_edge(i.pair[0], i.pair[1], bondOrder=i.order)
        return Atom

    def __get_surface_path__(self):
        surface_path = self.tree.xpath(
            "//ssip:SurfaceInformation/ssip:Surfaces/ssip:Surface", namespaces={"ssip": SSIP_NS, "cml": CML_NS})
        if len(surface_path) > 1:
            LOGGER.warn(
                "There are multiple sets of surface information in the file, only one is being read in")
        return surface_path

    def __get_equivalence_dict__(self):
        nm = iso.categorical_node_match('elementType', 'cap')
        GM = iso.GraphMatcher(self.network, self.network, node_match=nm)
        ds = []
        for subgraph in GM.subgraph_isomorphisms_iter():
            ds.append(subgraph)
        dn = {}
        list_k = ["a1"]
        for k in ds[0].keys():
            dn[k] = []
            for d in ds:
                if d[k] not in dn[k]:
                    dn[k].append(d[k])
        return dn

    def __get_equivalence_array_list__(self):
        list_array = []
        dict_list = self.__get_equivalence_dict__()
        for i in dict_list.values():
            if isinstance(i, list) == True:
                i.sort()
                a = tuple(i)
            elif isinstance(i, str) == True:
                a = tuple(i)
            else:
                a = None
            if any((a == x) for x in list_array):
                pass
            else:
                list_array.append(a)
        return list_array

    def __get_dict_average__(self):
        list_array = self.__get_equivalence_array_list__()
        dict_average = {}
        for j in list_array:
            dict_average[j] = []
            for i in j:
                ssip = self.atom_dict[i].neigh_ssip
                if isinstance(ssip, self.SsipSsip) == True:
                    dict_average[j].append(ssip.value)
            if len(dict_average[j]) != 0:
                dict_average[j] = np.average(dict_average[j])
            else:
                dict_average[j] = None
        dict_atom_average = {}
        for k, v in dict_average.items():
            for j in k:
                dict_atom_average[j] = v
        return dict_atom_average

    class SsipBond:
        """This is a class that helps to organise information regarding the connectivity described in the ssip.xml file
           Attributes
           ----------
           name : tuple 
               a tuple containing the name of the cml id of the atoms present in the bond 
           order : str
               the value of the bond order of the bond"""
        
        def __init__(self, bond_elem):
            bondpair = bond_elem.attrib["{{{}}}atomRefs2".format(CML_NS)]
            b1, b2 = bondpair.split(' ')
            self.pair = (b1, b2)
            self.order = bond_elem.attrib["{{{}}}order".format(CML_NS)]

    class SsipAtom:
        """This is a class that helps to organise information regarding an atom present in the ssip.xml file
           Attributes
           ----------
           name : str 
               the cml id of the atom
           elem : str 
               the chemical element of the atom  if a while loop restarts stop
           coord : np.array([float, float, float])
               a numpy array containing the x, y, z coordinates of the atom in nm."""

        def __init__(self, cml_elem, list_ssip):
            self.list_ssip = list_ssip
            self.name = cml_elem.attrib['{{{}}}id'.format(CML_NS)]
            self.elem = cml_elem.attrib['{{{}}}elementType'.format(CML_NS)]
            self.coord = np.array([0.1 * float(cml_elem.attrib['{{{}}}x3'.format(CML_NS)]),
                                   0.1 *
                                   float(
                                       cml_elem.attrib['{{{}}}y3'.format(CML_NS)]),
                                   0.1 * float(cml_elem.attrib['{{{}}}z3'.format(CML_NS)])])
            self.type = cml_elem.attrib['{{{}}}aipAtomType'.format(CML_NS)]
            self.neigh_ssip = self.get_neigh_ssip(self.list_ssip)
            self.neigh_ssip_name = self.get_neigh_ssip_name(self.list_ssip)

        def get_neigh_ssip(self, list_ssip):
            """A method used to represent information regarding an atom centered ssip
               Attributes
               ----------
               atom_ssip : SsipMolecule.SsipAtom
                   The same as SsipMolecule.SsipAtom with extra information about the value of the atom centered ssip."""
            ssip_value = 0
            count = 0
            atom_ssip = None
            for i in list_ssip:
                if i.neigh == self.name and abs(i.value) >= ssip_value:
                    atom_ssip = i
                    ssip_value = i.value
            return atom_ssip

        def get_atom_ssip(self, dict_average):
            if dict_average[self.name] == None:
                dis_min = 5
                for i in self.list_ssip:
                    dis = np.linalg.norm(self.coord - i.coord)
                    if dis < dis_min:
                        dis_min = dis
                        atom_ssip = i.value
            else:
                atom_ssip = dict_average[self.name]
            return atom_ssip

        def get_atom_ssip_no_average(self, list_ssip):
            atom_ssip = None
            if self.neigh_ssip == None:
                distance = 5
                for i in list_ssip:
                    dis = np.linalg.norm(self.coord - i.coord)
                    if dis < distance:
                        distance = dis
                        atom_ssip = i.value

            elif self.neigh_ssip != None:
                atom_ssip = self.neigh_ssip.value
            return atom_ssip

        def get_neigh_ssip_name(self, list_ssip):
            atom_ssip_name = None
            distance = 5
            for index, i in enumerate(list_ssip):
                dis = np.linalg.norm(self.coord - i.coord)
                if dis < distance:
                    distance = dis
                    atom_ssip_name = "M{}".format(index)
            return atom_ssip_name

    class SsipSsip:
        """This is a class that helps to organise information regarding an ssip present in the ssip.xml file
           Attributes
           ----------
           value : float 
               the value of the SSIP
           neigh : str 
               the cml id of the nearest atom to the ssip
           coord : np.array([float, float, float])
               a numpy array containing the x, y, z coordinates of the ssip in nm."""

        def __init__(self, ssip_elem):
            self.value = float(ssip_elem.attrib['{{{}}}value'.format(SSIP_NS)])
            self.fraction = float(
                ssip_elem.attrib['{{{}}}aipAreaFraction'.format(SSIP_NS)])
            self.isosurface = float(
                ssip_elem.attrib['{{{}}}isosurface'.format(SSIP_NS)])
            self.type = str(
                ssip_elem.attrib['{{{}}}aipAtomType'.format(SSIP_NS)])
            #try:    
            #    self.dual = str(ssip_elem.attrib['{{{}}}valueMax'.format(SSIP_NS)])
            #except:
            #    self.dual = 0
            self.neigh = ssip_elem.attrib["{{{}}}nearestAtomID".format(
                SSIP_NS)]
            self.coord = np.array([0.1 * float(ssip_elem.attrib['{{{}}}x3'.format(CML_NS)]),
                                   0.1 *
                                   float(
                                       ssip_elem.attrib['{{{}}}y3'.format(CML_NS)]),
                                   0.1 * float(ssip_elem.attrib['{{{}}}z3'.format(CML_NS)])])

        @staticmethod
        def compute_cross(c1, c2, c3):
            v1 = c2 - c1
            v2 = c3 - c1
            vcross = np.cross(v1, v2)
            mcross = np.linalg.norm(vcross)
            return mcross

        @staticmethod
        def compute_transformation_matrix(r12, r13, rcross):
            m = np.column_stack((r12, r13, rcross))
            try:
                inv = np.linalg.inv(m)
            except:
                inv = np.linalg.pinv(m)
            return inv

        def compute_weights(self, ssip_coord, r1, r12, r13, rcross):
            a = ssip_coord - r1
            transf_matrix = self.compute_transformation_matrix(
                r12, r13, rcross)
            return np.dot(transf_matrix, a)

        def get_anchors_weights(self, network, dict_atom, set_custom_anchors):
            self2 = copy.deepcopy(self)
            self2.anchor1 = None
            self2.anchor2 = None
            self2.anchor3 = None
            self2.weights = None

            anchor_list = []
            h_num = 0
            j = 0
            cml = dict_atom[self.neigh]
            if cml.elem == 'H':
                h_num += 1
            start = time.time()
            while len(anchor_list) < 2:
                j += 1
                for i in list(network.nodes):
                    if shortest_path_length(network, source=cml.name, target=i, weight=None, method='dijkstra') == j:
                        if dict_atom[i].elem == 'H':
                            h_num += 1
                        if h_num < 2 and dict_atom[i].elem == 'H':
                            anchor_list.append(i)
                        elif dict_atom[i].elem != 'H':
                            anchor_list.append(i)
                        if len(anchor_list) >= 2 and 0.0 < self.compute_cross(cml.coord, dict_atom[anchor_list[0]].coord, dict_atom[anchor_list[1]].coord) < 0.0001:
                            anchor_list.remove(anchor_list[1])
                end = time.time()
                if end-start > 1:
                    LOGGER.error(
                        "While loop timed out: Could not find the appropriate anchors")
                    return self2
            self2.anchor1 = dict_atom[self.neigh]
            self2.anchor2 = dict_atom[anchor_list[0]]
            self2.anchor3 = dict_atom[anchor_list[1]]
            if set_custom_anchors:
                self2.anchor1 = dict_atom[self.neigh]
                neigh = int(self.neigh[1:])
                if f"a{neigh+3}" in dict_atom.keys():
                    self2.anchor2 = dict_atom[f"a{neigh+3}"]
                else: 
                    overh = neigh+3 - len(dict_atom.keys())
                    self2.anchor2 = dict_atom[f"a{overh}"]
                if f"a{neigh+5}" in dict_atom.keys():
                    self2.anchor3 = dict_atom[f"a{neigh+5}"]
                else: 
                    overh = neigh+5 - len(dict_atom.keys())
                    self2.anchor3 = dict_atom[f"a{overh}"]

            r1 = self2.anchor1.coord
            r12 = self2.anchor2.coord-self2.anchor1.coord
            r13 = self2.anchor3.coord-self2.anchor1.coord
            rcross = np.cross(r12, r13)
            self2.weights = self.compute_weights(
                self.coord, r1, r12, r13, rcross)
                
            return self2

    class SsipSurface:
        """This is a class that helps to organise information regarding the surface of the molecule as described  in the ssip.xml file
           Attributes
           ----------
           tot : float 
               Total surface area
           pos : float 
               Positive surface area
           neg : float 
               Negative surface area
           elec_dens_iso : float 
               The electron denisty of the isosurface considered find the electrostatic potential
           n_meps : int."""

        def __init__(self, surface_elem):

            tot_path = surface_elem.xpath("ssip:TotalSurfaceArea", namespaces={
                                          "ssip": SSIP_NS, "cml": CML_NS})
            self.tot = tot_path[0].text
            if len(tot_path) > 1:
                LOGGER.warn(
                    "There are multiple sets of total surface information in the file, only one is being read in")

            pos_path = surface_elem.xpath("ssip:PositiveSurfaceArea", namespaces={
                                          "ssip": SSIP_NS, "cml": CML_NS})
            self.pos = pos_path[0].text
            if len(pos_path) > 1:
                LOGGER.warn(
                    "There are multiple sets of positive surface information in the file, only one is being read in")

            neg_path = surface_elem.xpath("ssip:NegativeSurfaceArea", namespaces={
                                          "ssip": SSIP_NS, "cml": CML_NS})
            self.neg = neg_path[0].text
            if len(neg_path) > 1:
                LOGGER.warn(
                    "There are multiple sets of negative surface information in the file, only one is being read in")
