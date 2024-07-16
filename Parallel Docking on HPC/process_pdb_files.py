#####################################
#AmirhosseinTaghavi              
#UF Scripps
#07/16/2023
#manage the parallel processing of PDB\
files using the multiprocessing module
#####################################
import os
import glob
import sys
import multiprocessing
from subprocess import run, CalledProcessError

def process_file(pdb_file, output_dir):
    try:
        run(['python', 'add_hydrogens_and_convert.py', pdb_file, output_dir], check=True)
    except CalledProcessError as e:
        print(f"Error while processing {pdb_file}: {e}")

def main(input_dir, output_dir, num_workers):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get a list of all PDB files in the directory
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))

    # Set up multiprocessing pool
    pool = multiprocessing.Pool(num_workers)
    
    # Process files in parallel
    for pdb_file in pdb_files:
        pool.apply_async(process_file, args=(pdb_file, output_dir))
    
    pool.close()
    pool.join()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_pdb_files.py <input_dir> <output_dir> <num_workers>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    num_workers = int(sys.argv[3])

    main(input_dir, output_dir, num_workers)

