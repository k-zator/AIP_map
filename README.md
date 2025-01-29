# AIP_map

Module for pairing AIPs in dimeric guest-host and protein-ligand complexes.

Author: Katarzyna Zator

#### Environment setup.

Clone the repository:

    git clone https://github.com/k-zator/AIP_map.git
    cd AIP_map

Then install this environment instead (based on python 3.6):

        conda env create -f environment.yml
        conda activate aip_map

#### Calculation.

The code is as a module so that the entire calculation is available as CLI.
The calculation requires two molecular structure input files and two AIP files. Molecular files should be in the PDB format, and AIP files can be calculated using the AIP module (https://github.com/k-zator/AIP). To calculate:

        python -m AIP_interaction_map.score -p1 compA.pdb -p2 compB.pdb -a1 aipA.xml -a2 aipB.xml -w A_B 

The module produces a folder (A_B) with the PDB output (A_B.pdb), a jmol script (A_B_script) to visualise the AIP-AIP interactions, and a CSV file (A_B.csv) that contains the list of AIP-AIP interactions, the relevant AIPs, and their 
contributions to the binding free energy.

#### Further options.

There are a series of additional parameters which can further refine the search:

  -h, --help            show this help message and exit
  --PDB_file1 PDB_file1, -p1 PDB_file1
                        PDB file with first of the interacting structures
  --PDB_file2 PDB_file2, -p2 PDB_file2
                        PDB file with second of the interacting structures
  --aip_file1 aip_file1, -a1 aip_file1
                        AIP file with first of the interacting structures
  --aip_file2 aip_file2, -a2 aip_file2
                        AIP file with second of the interacting structures
  --max_aip_dist MAX_AIP_DIST, -d MAX_AIP_DIST
                        default: 0.15, the maximum distance allowed for aip-aip contact to be considered an interaction
  --solvent SOLVENT, -s SOLVENT
                        default: chloroform, determines the solvent in which the interaction takes place, see full list in the paper SI
  --mbpa                
                        if present, uses the maximum bipartite algorithm rather than the custom branching algorithm
  --protein             
                        if present, it treats the host, second structure, as a protein and uses default AIP parametrisation for protein
  --bypass_at1          
                        if present, it treats the AIP-AIP contacts closer than 1 angstrom as preferential and does a first-pass of pairing for those alone, a speed-up for large systems
  --write WRITE, -w WRITE
                        if used, determines name for the jmol and csv file instead of {PDB1}_{PDB2}
  --lower LOWER, -low LOWER
                        default: 0.5 kJ/mol, AIP pairs with an absolute value that is less than this value are considered low and are plotted as small spheres
  --strong STRONG, -str STRONG
                        default: -5 kJ/mol, AIP pairs more negative than this value are considered strongly attractive and are plotted as dark green spheres

They are saved so that they do not need to be compiled every time but if deleted, but otherwise will just be written out again. The process length is proportional to molecular size. The AIP-AIP pairings can be visualised with:

        jmol A_B_script

#### Protein-ligand complexes.

It is possible to investigate protein-ligand complexes with the AIP_interaction_map module as the protein can be 
treated as an embedded series of amino acids whose AIPs have been precalculated. It requires producing a custom 
description of the protein binding site (primarily to minimise the calculation cost and ease visualisation) and
forcefield files for mapping the AIPs onto the structure with:

        python -m AIP_interaction_map.create_bs distance 2 -p protein.pdb -o .
        python -m AIP_interaction_map.create_ff.py -s ligand.xml -v -o .

This produces a binding_site.pdb, and ligand_vs.ff and ligand_custom.xml which contain the description of the protein
and which will be implicitly read by the AIP_interaction_map.score module. The command to produce the AIP-AIP pairings
is now:

        python -m AIP_interaction_map.score -p1 ligand.pdb -p2 binding_site.pdb -a1 ligand.xml --protein -w A_B 

Similarly, it produces the protein_ligand.pdb, jmol script, and a csv table of results.


Who do I talk to?
Any queries please contact Katarzyna Zator, kz265@cam.ac.uk

License
© Katarzyna Zator, Maria Chiara Storer, Christopher Hunter at the University of Cambridge
This is released under an AGPLv3 license for academic use.
