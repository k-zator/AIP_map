"""
Parses input for scoring from command line.
@author: Katarzyna Joanna Zator (kz265)
"""
import os
import pandas as pd
import logging
import argparse
from AIP_interaction_map.getScoring import Scoring
from AIP_interaction_map.create_jmol_vis import create_ip_vis
from AIP_interaction_map.PDB_preprocess import PDB_preprocess
logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


def create_parser():
    help_text = 'AIP pairing scoring algorithm to find interacting contacts made by AIPs within a parametric distance. \
                 Requires input of two molecule PDB files and two aip.xml files with complete pathes. \
                 For the Protein-ligand option, all protein AIPs are already encoded in  data/protein_aip_ff.xml. \
                 The function returns a csv table of AIP contacts and a jmol visualisation script named after the \
                 two PDB names: {PDB1}_{PDB2}_script.'
    sign_off = 'Author: Katarzyna Joanna Zator <kz265>'

    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)
    parser.add_argument(
        '--PDB_file1',
        '-p1',
        dest='PDB_file1',
        type=str,
        default=".",
        action='store',
        help='PDB file with first of the interacting structures',
        metavar='PDB_file1'
    )
    parser.add_argument(
        '--PDB_file2',
        '-p2',
        dest='PDB_file2',
        type=str,
        default=".",
        action='store',
        help='PDB file with second of the interacting structures',
        metavar='PDB_file2'
    )
    parser.add_argument(
        '--aip_file1',
        '-a1',
        dest='aip_file1',
        type=str,
        default=False,
        action='store',
        help='AIP file with first of the interacting structures',
        metavar='aip_file1'
    )
    parser.add_argument(
        '--aip_file2',
        '-a2',
        dest='aip_file2',
        type=str,
        default=False,
        action='store',
        help='AIP file with second of the interacting structures',
        metavar='aip_file2'
    )
    parser.add_argument(
        '--max_aip_dist',
        '-d',
        dest='max_aip_dist',
        type=float,
        default=0.15,
        action='store',
        help='default: 0.15, the maximum distance allowed for aip-aip contact to be considered an interaction'
    )
    parser.add_argument(
        '--solvent',
        '-s',
        dest='solvent',
        default="chloroform",
        action='store',
        help='default: chloroform, determines the solvent in which the interaction takes place,'
    )
    parser.add_argument(
        '--mbpa',
        dest='branch',
        default=True,
        action='store_false',
        help='if present, uses the maximum bipartite algorithm rather than the custom branching algorithm'
    )
    parser.add_argument(
        '--protein',
        dest='protein_host',
        default=False,
        action='store_true',
        help='if present, it treats the host, second structure, as protein and uses default AIP parametrisation'
    )
    parser.add_argument(
        '--bypass_at1',
        dest='bypass_at1',
        default=False,
        action='store_true',
        help='if present, it treats the aip-aip contacts closer than 1 A as preferential and does a  \
              first-pass of pairing for those alone, a speed-up for large systems'
    )
    parser.add_argument(
        '--write',
        '-w',
        dest='write',
        default=False,
        action='store',
        help='if used, determines name for the jmol and csv file instead of {PDB1}_{PDB2}'
    )
    parser.add_argument(
        '--lower',
        '-low',
        dest='lower',
        type=float,
        default=0.5,
        action='store',
        help='default: 0.5 kJ/mol, AIP pairs with an absolute value that is less than \
              this value are considered low and are plotted as small spheres'
    )
    parser.add_argument(
        '--strong',
        '-str',
        dest='strong',
        type=float,
        default=-5,
        action='store',
        help='default: -5 kJ/mol, AIP pairs more negative than this value \
              are considered strongly attractive and are plotted as dark green spheres'
    )    
    parser.set_defaults(func=scoring)
    return parser

def scoring(args):
    """Creates a folder structure for the pairing function: a main forlder and one for each interacting
       species so that name formatting is unique. Then runs scoring and creates a visualisation."""
    comp1 = args.PDB_file1.split("/")[-1].split(".")[0]
    comp2 = args.PDB_file2.split("/")[-1].split(".")[0]
    if args.aip_file1 != False:
        aip_name1 = args.aip_file1.split("/")[-1].split(".")[0]
    else:
        path_to_guest_folder = os.path.dirname(args.PDB_file1)
        args.aip_file1 = f"{path_to_guest_folder}/ssip.xml"
        aip_name1 = "ssip"

    if args.protein_host == False and args.aip_file2 != False:
        PDB_preprocess(args.PDB_file1, args.PDB_file2)
        aip_name2 = args.aip_file2.split("/")[-1].split(".")[0]
    elif args.protein_host == False:
        PDB_preprocess(args.PDB_file1, args.PDB_file2)
        path_to_host_folder = os.path.dirname(args.PDB_file2)
        args.aip_file2 = f"{path_to_host_folder}/ssip.xml"
        aip_name2 = "ssip"
    else:
        a2 = False

    if args.write != False:
        os.system(f"mkdir {args.write}")
        if args.protein_host != False:
            os.system(f"cp {args.aip_file1} {args.PDB_file1} {args.PDB_file2} {args.write}")
            p1 = f"{args.write}/{comp1}.pdb"
            p2 = f"{args.write}/{comp2}.pdb"
            a1 = f"{args.write}/{aip_name1}.xml"
        else:
            os.system(f"mkdir {args.write}/{comp1}")
            os.system(f"mkdir {args.write}/{comp2}")
            os.system(f"cp {args.aip_file1} {args.PDB_file1} {args.write}/{comp1}")
            os.system(f"cp {args.aip_file2} {args.PDB_file2} {args.write}/{comp2}")
            p1 = f"{args.write}/{comp1}/{comp1}.pdb"
            p2 = f"{args.write}/{comp2}/{comp2}.pdb"
            a1 = f"{args.write}/{comp1}/{aip_name1}.xml"
            a2 = f"{args.write}/{comp2}/{aip_name2}.xml"
    else:
        os.system(f"mkdir {comp1}_{comp2}")
        if args.protein_host != False:
            os.system(f"cp {args.aip_file1} {args.PDB_file1} {args.PDB_file2} {comp1}_{comp2}")
            p1 = f"{comp1}_{comp2}/{comp1}.pdb"
            p2 = f"{comp1}_{comp2}/{comp2}.pdb"
            a1 = f"{comp1}_{comp2}/{aip_name1}.xml"
        else:
            os.system(f"mkdir {comp1}_{comp2}/{comp1}")
            os.system(f"mkdir {comp1}_{comp2}/{comp2}")
            os.system(f"cp {args.aip_file1} {args.PDB_file1} {comp1}_{comp2}/{comp1}")
            os.system(f"cp {args.aip_file2} {args.PDB_file2} {comp1}_{comp2}/{comp2}")
            p1 = f"{comp1}_{comp2}/{comp1}/{comp1}.pdb"
            p2 = f"{comp1}_{comp2}/{comp2}/{comp2}.pdb"
            a1 = f"{comp1}_{comp2}/{comp1}/{aip_name1}.xml"
            a2 = f"{comp1}_{comp2}/{comp2}/{aip_name2}.xml"
    
    k = Scoring(p1, p2, aip_guest=a1, aip_host=a2, solvent=args.solvent,
                protein_host=args.protein_host, branch=args.branch,
                max_aip_dist=args.max_aip_dist, bypass_at1=args.bypass_at1)
    output_df = k.interaction_df[["L","L_type","L_value", "R", "R_type", "R_value", "Frac", "ddG_sc"]]
    output_df = output_df.round(1)
    output_df.L += 1
    if args.protein_host == False:
        output_df.R -= len(k.state.ligand_vs_indices) - 2
    else:
        R_cont_ind_dict = {}
        l = len(k.state.ligand_atom_indices) + 2
        for ind, i in enumerate(k.state.residue_atom_indices):
            R_cont_ind_dict[i] = ind + l
        output_df.R = pd.Series([R_cont_ind_dict[r] for r in output_df.R])
    output_df = output_df.rename(columns={"L": f"{comp1}", "R": f"{comp2}",
                              "L_type": f"{comp1}_type", "R_type": f"{comp2}_type",
                              "L_value": f"{comp1}_value", "R_value": f"{comp2}_value",
                              "L_frac": "f", "ddG_sc": "ddG (kJ/mol)"})

    if args.write != False:
        output_df.to_csv(f"{args.write}/{args.write}.csv")
        create_ip_vis(k, f"{args.write}/{args.write}", lower=args.lower, strong=args.strong)
    else:
        comp1 = args.PDB_file1.split("/")[-1].split(".")[0]
        comp2 = args.PDB_file2.split("/")[-1].split(".")[0]
        output_df.to_csv(f"{comp1}_{comp2}/{comp1}_{comp2}.csv")
        create_ip_vis(k, f"{comp1}_{comp2}/{comp1}_{comp2}",  lower=args.lower, strong=args.strong)

def main():
    """Main function. Returns logger outputs if encounters failures."""
    parser = create_parser()
    args = parser.parse_args()
    _ = args.func(args)
    LOGGER.info("Scoring completed. Written AIP pairings and visualisation.")

if __name__ == "__main__":
    main()
