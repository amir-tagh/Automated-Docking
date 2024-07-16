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
            
            # Check for existing hydrogens
            hydrogen_count = cmd.count_atoms("elem H")
            if hydrogen_count > 0:
                print(f"File {pdb_file} already contains hydrogens. Skipping.")
                cmd.delete("all")
                continue
            
            # Get the base name of the PDB file (without the directory and extension)
            base_name = os.path.basename(pdb_file).split('.')[0]
            
            # Add missing hydrogens
            cmd.h_add()
            
            # Save the modified PDB file to the output directory
            output_file = os.path.join(output_dir, f"{base_name}_with_hydrogens.pdb")
            cmd.save(output_file, selection="all")  # Use selection="all" to save the entire structure
            
            # Clear the structure to prepare for the next file
            cmd.delete("all")
        except CmdException as e:
            print(f"CmdException encountered while processing {pdb_file}: {e}")
        except Exception as e:
            print(f"Unexpected error encountered while processing {pdb_file}: {e}")

    print("Hydrogen addition completed for all PDB files.")
    cmd.quit()


# Define the directory containing the PDB files
pdb_dir = "/home/amir/Projects/DOCKING-OF-HARIBOSS-AI-COLLABORATION/HARIBOSS-PDB-NEW/hariboss_pdb/rna"
output_dir = "/home/amir/Projects/DOCKING-OF-HARIBOSS-AI-COLLABORATION/HARIBOSS-PDB-NEW/hariboss_pdb/rna_added_H"

add_hydrogens_to_pdbs(pdb_dir, output_dir)
