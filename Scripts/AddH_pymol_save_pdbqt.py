import os
import glob
from pymol import cmd, finish_launching, CmdException

def add_hydrogens_to_pdbs(pdb_dir, output_dir):
    # Initialize PyMOL in quiet mode (headless)
    finish_launching(['pymol', '-qc'])

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get a list of all PDB files in the directory
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

    # Loop over each PDB file
    for pdb_file in pdb_files:
        try:
            # Load the PDB file
            cmd.load(pdb_file)
            
            # Get the base name of the PDB file (without the directory and extension)
            base_name = os.path.basename(pdb_file).split('.')[0]
            
            # Add missing hydrogens
            cmd.h_add()
            
            # Save the modified structure as a PDBQT file
            output_file = os.path.join(output_dir, f"{base_name}.pdbqt")
            cmd.save(output_file, selection="all")  # Save the entire structure in PDBQT format
            
            # Clear the structure to prepare for the next file
            cmd.delete("all")
        except CmdException as e:
            print(f"CmdException encountered while processing {pdb_file}: {e}")
        except Exception as e:
            print(f"Unexpected error encountered while processing {pdb_file}: {e}")

    print("Hydrogen addition and conversion to PDBQT completed for all PDB files.")
    cmd.quit()

# Define the directory containing the PDB files
pdb_dir = "/path/to/your/pdb_files"
output_dir = "/path/to/save/processed_files"

add_hydrogens_to_pdbs(pdb_dir, output_dir)

