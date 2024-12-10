from AIP_interaction_map.binding_site_distance import protein_neighbours_pdb
from AIP_interaction_map.binding_site import BindingSite
import argparse
import sys

def create_parser():
    parser = argparse.ArgumentParser()
    help_text = 'Output the SUrface Area,  nearest atom ID, the value of the SSIP'
    sign_off = 'Author: MCStorer <mcs92>'
    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)
    parser.add_argument(
        '--protein',
        '-p',
        dest='protein',
        default=None, 
        type=str,
        action='store',
        help='Path to the protonated pdb file.', 
        metavar='a'
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
    subparsers = parser.add_subparsers(help='sub-command help', dest="command")
    parser_template = subparsers.add_parser('template', 
            help='Obtain the binding site from a non-protonated binding site file.')
    parser_distance = subparsers.add_parser('distance', 
            help='obtain the bidning site from a cutoff distance from the ligand')
    parser_template.add_argument(
                '--tf',
                '-t',
                dest='tf',
                default=None, 
                type=str, 
                help='Path to template file.'
                )
    parser_distance.add_argument(
                '--cutoff',
                '-c',
                dest='cutoff',
                default=0.5, 
                type=float, 
                help='Cutoff distance to define binding site.'
                )
    parser_distance.add_argument(
                '--ligand',
                '-l',
                dest='ligand',
                default=None, 
                type=str, 
                help='Path to ligand file.'
                )
    return parser

def create_bs(arguments):
    if arguments.command == 'template':
        binding_site = BindingSite(arguments.protein, arguments.tf)
        binding_site.write_file(f"{arguments.out}/binding_site.pdb")
    elif arguments.command == 'distance':
        protein_neighbours_pdb(
                arguments.ligand, 
                arguments.protein, 
                f"{arguments.out}/binding_site.pdb", 
                arguments.cutoff)
        
if __name__ == "__main__":
    # execute only if run as a script
    arguments = create_parser().parse_args(sys.argv[1:])
    create_bs(arguments)

