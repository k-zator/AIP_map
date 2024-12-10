
def prepare_PDB(PDB_file, mol_name=None):
    """Read in PDB file and set a distinct molecule name and format atom string to fit the pattern:
       'ATOM      1  C1  CPD A   1 <xyz>  C"""
    with open(PDB_file) as f:
        data = f.readlines()
        f.close()

    if sum([line.find("HETATM") for line in data]) == -len(data):
        #if HETATM is not present, hence ATOM is correct
        first_atom_index = [line.find("ATOM") for line in data].index(0)
    else:
        data = [line.replace("HETATM", "ATOM  ") if (
            line.find("HETATM") >= 0) else line for line in data]
        first_atom_index = [line.find("ATOM") for line in data].index(0)

    line = data[first_atom_index]
    if line[:13] == "ATOM      1  " and line[20:26] == " A   1":
        #correct line
        pdb_mol_name = data[2].split()[3]
    elif line[:13] == "ATOM      1  ":
        data = [line[:20] + " A   1" + line[26:]
                if (line.find("ATOM") >= 0) else line for line in data]
        pdb_mol_name = data[2].split()[3]
    else:
        print("Watch out, there's an issue with the PDB:", PDB_file)
        return None

    if mol_name == None:  #else keep custom defined one
        full_name = PDB_file.split("/")[-1].split(".")[0]
        mol_name = full_name[0]+full_name[-2:]
    data = [line.replace(pdb_mol_name, mol_name) if (
        line.find("ATOM") >= 0) else line for line in data]

    with open(PDB_file, 'w') as fw:
        fw.writelines(data)

    return mol_name


def PDB_preprocess(PDB_file1, PDB_file2=False):
    """Read in file and make sure the names of molecules contained in them are distinct,
       by making them a function of their PDB file name. Secondly, prepare combined.pdb
       for each and their respective ssip information xml files. It requies a file
       ligand_vs_ff.xml in the same folder as the PDB."""
    mol1_name = prepare_PDB(PDB_file1)
    if PDB_file2 is not False:
        mol2_name = prepare_PDB(PDB_file2)
        if mol1_name == mol2_name:
            prepare_PDB(PDB_file2, mol_name="CPD")

        if mol1_name == None or mol2_name == None:
            print("PDB preprocessing failed. Please check the input PDB files")
