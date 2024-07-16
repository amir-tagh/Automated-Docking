#####################################
#AmirhosseinTaghavi              
#UF Scripps
#07/16/2023
#####################################
#add the missinh hydrogens with pymol
#convert the pdb to pdbqt with obabel
import os
import subprocess
from pymol import cmd, finish_launching, CmdException

def add_hydrogens_to_pdb(input_file, output_dir):
    # Initialize PyMOL in quiet mode (headless)
    finish_launching(['pymol', '-qc'])

    try:
        # Load the PDB file
        cmd.load(input_file)
        
        # Get the base name of the PDB file (without the directory and extension)
        base_name = os.path.basename(input_file).split('.')[0]
        
        # Add missing hydrogens
        cmd.h_add()
        
        # Save the modified structure as a PDB file
        temp_pdb_file = os.path.join(output_dir, f"{base_name}.pdb")
        cmd.save(temp_pdb_file, selection="all")  # Save the entire structure
        
        # Convert PDB to PDBQT using OpenBabel
        output_pdbqt_file = os.path.join(output_dir, f"{base_name}.pdbqt")
        obabel_convert(temp_pdb_file, output_pdbqt_file)
        
        # Clear the structure to prepare for the next file
        cmd.delete("all")
    except CmdException as e:
        print(f"CmdException encountered while processing {input_file}: {e}")
    except Exception as e:
        print(f"Unexpected error encountered while processing {input_file}: {e}")
    finally:
        cmd.quit()

def obabel_convert(input_file, output_file):
    # Command to convert PDB to PDBQT using OpenBabel
    obabel_cmd = f"obabel {input_file} -O {output_file} -opdbqt"

    # Execute the command using subprocess
    try:
        subprocess.run(obabel_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while executing OpenBabel: {e}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python add_hydrogens_and_convert.py <input_file> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    add_hydrogens_to_pdb(input_file, output_dir)

