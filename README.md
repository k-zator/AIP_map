# hunter_scoring

Pre-release code for pairing AIPs in guest-host complexes.
Author: Katarzyna Zator, kz265

#### Environment setup.

Clone the repository:

    git clone git@gitlab.developers.cam.ac.uk:ch/hunter/hunter_scoring.git
    cd hunter_scoring

Then install this environment instead (based on python 3.6):

        conda create --name openmm --file spec-file.txt
        conda activate openmm

#### Code running.

The scoring algorithm is currently set up to work in an IDE: jupyter notebook.
To open it:

        jupyter-notebook

There are a series of available notebooks, or you can create your own.
In any case, the key command is:

        from ssipscore.getScoring import Scoring
        self = Scoring(PDB_file1, PDB_file2)

where PDB files need to have their full paths specified and they cannot be in the same folder.
It is an expectation that the ssip.xml file (named so) is also present in that folder and will
be used as the AIP source. The order technically matters as the first mol becomes ligand (L), and the latter residue (R).

#### Further options.

There are a series of default-specified parameters:

1. max_aip_dist=0.10 - max distance between AIPs to allow a contact, it has a very precaious scaling hence do not go beyond 0.20,

2. from_PDB=True - both guest and host are specified by PDB and ssip.xml explicitly. Otherwise, it assumes protein-host and requires only the PDB (prepped with binding_site.py to match atom types) and uses pre-calculated AIPs,

3. bypass_at1=False - an alternative algorithm for searching for AIP contacts with primary search at 0.08 nm, so that those have the algoritm performs three concentric searches. First at 0.08, then at max_aip_dist, then for remainder at max_aip_dist+ext_upon_desolv,

4. solvent="chloroform" - solvent in which the interactions are considered, currently the choice is between chloroform and water

5. probe_radius=0.2 - distance beyond the vdW radius at which SASAs are calculated,

6. frac_to_desolv=0.5 - fraction of desolvated SASA at which the atom is considered desolvated, this is per atom so far hence 2nd
   site corresponds to a constant of 0.7 (when available) and it does not reconise full desolvation of particular non-contiguous faces,

7. ext_upon_desolv=0.03 - extension to the max_aip_dist at which other contacts within the desolvation cavity are searched for.

The code does create a DataDicts.txt and combined.pdb files in each respective folder for guest_host pairing, and ligand_vs_ff.xml and ligand_vs_custom.xml for guest_protein option.
They are saved so that they do not need to be compiled every time but if deleted, but otherwise will just be written out again. The process length is proportional to molecular size.

#### Useful commands.

Of notable interest will be particular objects:

        #Dataframe of present interactions
        self.interaction_df

which also can be saved: self.interaction_df.to_csv(f"{name}.csv").
It contains information about the atom-atom and AIP-AIP contacts, namely the distances, indices (L, R for atoms, or \_AIP), Atom types, AIP fractions, AIP values, the contribution to overall binding: ddG (and ddG_sc: scaled by AIP fractions)

The indices are maintained througout the algorithm but the R molecule will have its shifted to ensure they are all unique. This might make it slightly harder to see the contacts, hence for visualisation:

        from ssipscore.create_jmol_vis import create_jmol_vis, create_ip_vis
        #visualisation of all AIPs without interacting ones bein colour-coded:
        create_jmol_vis
        #visualisation of the AIP contacts (as mid-point between interacting AIPs)
        create_ip_vis(self, f"{name}")
        #this will create a .pdb and a jmol _script. To see the contacts, run
        jmol {name}_script

Additionally, because for the guest-host contacts the interacting surface is frequently concave, more than just the interacting AIPs should be desolvated. The total desolvation energy (of both interacting and hidden but not interacting AIPs), identity of the non-interacting but hidden AIPs, and any contacts that could be made with those (i.e. at the extended aip-aip distance) are available.

        #Desolvation energy for all sites with sufficient SASA change
        self.desolv_energy
        #AIPs that are desolvated upon complexation but did not have a contact found
        self.desolvated_AIPs
        #Additional interactions present beyond cut-off but within the desolvated cavity
        self.non_polar_add
