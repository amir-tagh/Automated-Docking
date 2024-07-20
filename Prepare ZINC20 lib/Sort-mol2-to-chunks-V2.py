import os
import glob

def count_mol2_structures(file_path):
    structure_count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == '@<TRIPOS>MOLECULE':
                structure_count += 1
    return structure_count

def combine_mol2_files(input_dir, chunk_size=30000):
    # Get list of all .mol2 files in the directory
    mol2_files = sorted(glob.glob(os.path.join(input_dir, '*.mol2')))
    
    if not mol2_files:
        print("No .mol2 files found in the directory.")
        return

    chunk_counter = 0
    current_structure_count = 0
    output_file = None

    # Create directory to store combined files
    combined_dir = os.path.join(input_dir, 'combined')
    os.makedirs(combined_dir, exist_ok=True)

    for mol2_file in mol2_files:
        num_structures = count_mol2_structures(mol2_file)

        if current_structure_count + num_structures > chunk_size:
            if output_file:
                output_file.close()
            chunk_counter += 1
            current_structure_count = 0
            output_file_path = os.path.join(combined_dir, f'combined_{chunk_counter:03d}.mol2')
            output_file = open(output_file_path, 'w')

        with open(mol2_file, 'r') as f:
            output_file.write(f.read())
            current_structure_count += num_structures

    if output_file:
        output_file.close()

    print(f"Combined files are saved in: {combined_dir}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Combine .mol2 files in a directory into chunks of specified size.')
    parser.add_argument('input_dir', type=str, help='Path to the directory containing .mol2 files.')
    parser.add_argument('--chunk_size', type=int, default=30000, help='Number of structures per combined file.')

    args = parser.parse_args()
    combine_mol2_files(args.input_dir, args.chunk_size)

