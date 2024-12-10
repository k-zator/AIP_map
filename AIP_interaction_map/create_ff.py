import argparse
from AIP_interaction_map.ssip_molecule import SsipMolecule
import sys
from AIP_interaction_map.constants import FF_VS_PATH


def create_parser():
    parser = argparse.ArgumentParser()
    help_text = 'Create forcefield files for OpenMM from ssip.xml files'
    sign_off = 'Author: MCStorer <mcs92>'
    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)
    parser.add_argument(
        '--ssip',
        '-s',
        dest='ssip',
        type=str,
        default=None,
        action='store',
        help='SSIP file to be converted into force fields',
        metavar='a'
        )
    parser.add_argument(
        '--vs',
        '-v',
        dest='vs',
        default=False,
        action='store_true',
        help='Boolean. Add SSIPs as virtual sites.\
              If False, the description is only atom centered',
        )
    parser.add_argument(
        '--output',
        '-o',
        dest='out',
        type=str,
        default=".",
        action='store',
        help='filepath for ligand and custom ff files',
        metavar='a'
        )
    parser.add_argument(
        '--ligand_only',
        '-l',
        dest='l',
        default=False,
        action='store_true',
        help='if used, this flag indicates that you just want to look at a ligand,\
                and does not consider the protein in the forecefield creation',
        )
    parser.add_argument(
        '--not_averaged',
        '-n',
        dest='n',
        default=True,
        action='store_false',
        help='if used, this flag indicates that you don\'t \
                want the atom centered ff to be averaged.'
        )
    return parser 

def create_ff(arguments):
    """Create forcefield files from ssip file."""
    if type(arguments.ssip) != None:
        mol = SsipMolecule(arguments.ssip, averaged_atom_ssip = arguments.n)
        if arguments.vs==True and arguments.l==True:
            mol.write_ligand_xml(mol.vs_tree, f"{arguments.out}/ligand_vs_ff.xml")
            mol.write_custom_xml(mol.vs_tree, f"{arguments.out}/ligand_vs_custom.xml")
        elif arguments.vs==True and arguments.l==False:
            mol.write_ligand_xml(mol.vs_tree, f"{arguments.out}/ligand_vs_ff.xml")
            mol.write_custom_xml(mol.vs_tree, f"{arguments.out}/ligand_vs_custom.xml",\
            FF_VS_PATH)
        elif arguments.vs==False and arguments.l==True:
            mol.write_ligand_xml(mol.atom_tree, f"{arguments.out}/ligand_atom_ff.xml")
            mol.write_custom_xml(mol.atom_tree, f"{arguments.out}/ligand_atom_custom.xml")
        else:
            mol.write_ligand_xml(mol.atom_tree, f"{arguments.out}/ligand_atom_ff.xml")
            mol.write_custom_xml(mol.atom_tree, f"{arguments.out}/ligand_atom_custom.xml",\
            FF_VS_PATH)

def create_ff_int(ssip, out):
    """Internal create forcefield function use hence limited to use in hunter_scoring.
       By default, uses explicit non-averaged virtual sites for ligand and protein."""
    if type(ssip) != None:
        mol = SsipMolecule(ssip, averaged_atom_ssip = False)
        mol.write_ligand_xml(mol.vs_tree, f"{out}/ligand_vs_ff.xml")
        mol.write_custom_xml(mol.vs_tree, f"{out}/ligand_vs_custom.xml", FF_VS_PATH)

if __name__ == "__main__":
    # execute only if run as a script
    arguments = create_parser().parse_args(sys.argv[1:])
    create_ff(arguments)
