############################################################################
#script extracts nonstandard residues from a pdb file and calculated the\
#molecular weight with RDkit. These residues can later be replaced with standard\
#residues with close MW.
#Standard Residues Embedded: The list of standard residues is embedded\
#in the script under the main function. updare as required.
#Amirhossein Taghavi
#UF Scripps
#07/27/2024
############################################################################

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def calculate_molecular_weight(pdb_file, residue_name):
    # Load the PDB file
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False)
    if not mol:
        raise ValueError(f"Could not read PDB file: {pdb_file}")

    # Find the residue with the given name
    residue_atoms = [atom for atom in mol.GetAtoms() if atom.GetPDBResidueInfo().GetResidueName().strip() == residue_name]

    if not residue_atoms:
        raise ValueError(f"Residue {residue_name} not found in the PDB file.")

    # Create a new molecule for the residue
    residue_mol = Chem.RWMol()
    atom_mapping = {}

    for atom in residue_atoms:
        new_atom = Chem.Atom(atom.GetAtomicNum())
        new_idx = residue_mol.AddAtom(new_atom)
        atom_mapping[atom.GetIdx()] = new_idx

    # Add bonds to the residue molecule
    added_bonds = set()
    for atom in residue_atoms:
        atom_idx = atom.GetIdx()
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in atom_mapping:
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                bond_pair = tuple(sorted((atom_mapping[atom_idx], atom_mapping[neighbor_idx])))
                if bond_pair not in added_bonds:
                    residue_mol.AddBond(atom_mapping[atom_idx], atom_mapping[neighbor_idx], bond.GetBondType())
                    added_bonds.add(bond_pair)

    # Finalize the molecule
    residue_mol = residue_mol.GetMol()

    # Sanitize the molecule and calculate implicit valence
    Chem.SanitizeMol(residue_mol)
    AllChem.Compute2DCoords(residue_mol)

    # Calculate the molecular weight
    mol_weight = AllChem.CalcExactMolWt(residue_mol)
    return mol_weight

def identify_non_standard_residues(pdb_file, standard_residues):
    # Load the PDB file
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False)
    if not mol:
        raise ValueError(f"Could not read PDB file: {pdb_file}")

    non_standard_residues = set()
    for atom in mol.GetAtoms():
        residue_name = atom.GetPDBResidueInfo().GetResidueName().strip()
        if residue_name not in standard_residues:
            non_standard_residues.add(residue_name)

    return list(non_standard_residues)

def main(args):
    standard_residues = ["A", "G", "C", "U", "HOH", "GLN", "GLU", "GLY", "HIS"]
    
    if args.action == "calculate":
        mol_weight = calculate_molecular_weight(args.pdb_file, args.residue_name)
        print(f"The molecular weight of the residue {args.residue_name} is: {mol_weight:.4f}")
    elif args.action == "identify":
        non_standard_residues = identify_non_standard_residues(args.pdb_file, standard_residues)
        print("Non-standard residues found in the PDB file:")
        for res in non_standard_residues:
            print(res)
            mol_weight = calculate_molecular_weight(args.pdb_file, res)
            print(f"The molecular weight of the non-standard residue {res} is: {mol_weight:.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files to calculate molecular weights or identify non-standard residues.")
    subparsers = parser.add_subparsers(dest="action", help="Action to perform")

    # Subparser for calculating molecular weight
    parser_calculate = subparsers.add_parser("calculate", help="Calculate the molecular weight of a residue.")
    parser_calculate.add_argument("pdb_file", type=str, help="Path to the PDB file.")
    parser_calculate.add_argument("residue_name", type=str, help="Name of the residue to calculate the molecular weight for.")

    # Subparser for identifying non-standard residues
    parser_identify = subparsers.add_parser("identify", help="Identify non-standard residues in a PDB file.")
    parser_identify.add_argument("pdb_file", type=str, help="Path to the PDB file.")

    args = parser.parse_args()
    main(args)

