###########################################################################################################################################
Amirhossein Taghavi
07/20/2024
#Convert SMILES to SDF: The script converts SMILES to SDF files in parallel using Open Babel and saves them in the sdf_files directory.
#Apply GAFF Force Field: It applies the GAFF force field to these SDF files.
#Calculate NPRs: The script calculates NPR1 and NPR2 for each molecule.
#Plot NPRs: It plots NPR1 vs NPR2 using Matplotlib.
###########################################################################################################################################

import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from multiprocessing import Pool, cpu_count

def run_obabel(smiles, index, output_dir):
    input_file = os.path.join(output_dir, f"temp_{index}.smi")
    sdf_file = os.path.join(output_dir, f"temp_{index}.sdf")
    gaff_sdf_file = os.path.join(output_dir, f"temp_{index}_gaff.sdf")

    with open(input_file, 'w') as f:
        f.write(smiles)

    subprocess.run(['obabel', '-i', 'smi', input_file, '-o', 'sdf', '-O', sdf_file, '--gen3d'])
    subprocess.run(['obabel', sdf_file, '-O', gaff_sdf_file, '--minimize', '--ff', 'GAFF'])

    return gaff_sdf_file

def convert_smiles_to_sdf_parallel(smiles_list, num_cpus, output_dir):
    with Pool(num_cpus) as pool:
        sdf_files = pool.starmap(run_obabel, [(smiles, i, output_dir) for i, smiles in enumerate(smiles_list)])
    return sdf_files

def compute_inertia_tensor(coords):
    com = np.mean(coords, axis=0)
    coords_centered = coords - com
    
    inertia_tensor = np.zeros((3, 3))
    for i in range(coords_centered.shape[0]):
        x, y, z = coords_centered[i]
        inertia_tensor[0, 0] += (y**2 + z**2)
        inertia_tensor[1, 1] += (x**2 + z**2)
        inertia_tensor[2, 2] += (x**2 + y**2)
        inertia_tensor[0, 1] -= x*y
        inertia_tensor[0, 2] -= x*z
        inertia_tensor[1, 2] -= y*z
    
    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]
    
    return inertia_tensor

def calculate_principal_moments_of_inertia(mol):
    mol = Chem.AddHs(mol)

    try:
        AllChem.EmbedMolecule(mol, randomSeed=42, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        AllChem.UFFOptimizeMolecule(mol, maxIters=200, vdwThresh=10.0, confId=-1, ignoreInterfragInteractions=False)
        
        conf = mol.GetConformer()
        coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
        coords = np.array(coords)
        
        inertia_tensor = compute_inertia_tensor(coords)
        moments_of_inertia = np.linalg.eigvals(inertia_tensor)
        sorted_moments = sorted(moments_of_inertia)
        
        if len(sorted_moments) >= 3 and sorted_moments[0] > 0:
            npr1 = sorted_moments[0] / sorted_moments[1]  # PMI1 / PMI2
            npr2 = sorted_moments[1] / sorted_moments[2]  # PMI2 / PMI3
            return npr1, npr2
        else:
            print("Invalid principal moments of inertia.")
            return None, None
    except Exception as e:
        print(f"Error calculating NPR for molecule: {e}")
        return None, None

def plot_npr(npr_ratios, output_plot):
    valid_ratios = [(npr1, npr2) for npr1, npr2 in npr_ratios if npr1 is not None and npr2 is not None]
    
    if not valid_ratios:
        print("No valid NPR ratios to plot.")
        return
    
    npr1, npr2 = zip(*valid_ratios)

    fig, ax = plt.subplots(figsize=(12, 8))

    # Draw upside-down triangle
    triangle = np.array([[0, 0], [1, 0], [0.5, -np.sqrt(3)/2], [0, 0]])
    ax.plot(triangle[:, 0], triangle[:, 1], 'k-')
    
    # Convert NPR ratios to triangle coordinates
    x = 0.5 * (2 * np.array(npr1) + np.array(npr2)) / (np.array(npr1) + np.array(npr2) + 1)
    y = -(np.sqrt(3) / 2) * np.array(npr2) / (np.array(npr1) + np.array(npr2) + 1)
    
    # Create a gradient color map based on NPR1 and NPR2 values
    colors = np.sqrt(np.array(npr1)**2 + np.array(npr2)**2)

    sc = ax.scatter(x, y, s=50, edgecolor='k', c=colors, cmap='viridis', alpha=0.7)  # Reduced dot size to 50
    plt.colorbar(sc, ax=ax, label='Gradient based on NPR values')
    
    # Annotate corners with slight offset to avoid overlapping
    ax.text(0, 0.05, 'Sphere', verticalalignment='bottom', horizontalalignment='left', fontsize=12)
    ax.text(1, 0.05, 'Rod', verticalalignment='bottom', horizontalalignment='right', fontsize=12)
    ax.text(0.5, -np.sqrt(3)/2 - 0.05, 'Disc', verticalalignment='top', horizontalalignment='center', fontsize=12)
    
    # Label axes
    ax.set_xlabel('NPR1 (PMI1 / PMI2)')
    ax.set_ylabel('NPR2 (PMI2 / PMI3)')

    ax.set_aspect('equal')
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-np.sqrt(3)/2 - 0.1, 0.1)
    
    plt.title('Normalized Principal Moments of Inertia Ratios (NPR1 vs NPR2)')
    plt.savefig(output_plot)
    #plt.show()

def main(input_file, num_cpus, output_plot):
    output_dir = 'sdf_files'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    smiles_list = []
    npr_ratios = []

    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 1:  # Ensure there is at least one column
                smiles = parts[0].strip()  # Read SMILES from the first column
                smiles_list.append(smiles)
            else:
                print("Line does not contain enough columns.")

    sdf_files = convert_smiles_to_sdf_parallel(smiles_list, num_cpus, output_dir)

    for sdf_file in sdf_files:
        suppl = Chem.SDMolSupplier(sdf_file)
        for mol in suppl:
            if mol is None: continue
            npr_ratio = calculate_principal_moments_of_inertia(mol)
            npr_ratios.append(npr_ratio)

    plot_npr(npr_ratios, output_plot)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and plot Normalized Principal Moments of Inertia Ratios (NPR1 vs NPR2) for molecules.')
    parser.add_argument('input_file', type=str, help='Path to the tab-separated file containing SMILES strings.')
    parser.add_argument('--num_cpus', type=int, default=cpu_count(), help='Number of CPU cores to use for parallel processing.')
    parser.add_argument('--output_plot', type=str, default='npr_plot.png', help='Path to save the NPR scatter plot.')
    args = parser.parse_args()
    main(args.input_file, args.num_cpus, args.output_plot)

