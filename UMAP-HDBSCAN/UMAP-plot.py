######################################################################
#Parse the SMILES strings.\
#Convert the SMILES strings to molecular descriptors or fingerprints.\
#Use UMAP to reduce the dimensionality of these descriptors.\
#Plot the UMAP embeddings.\
#A list of valid SMILES strings (valid_smiles) is maintained\
#to keep track of successfully processed SMILES.

#####################################################################
#Amirthossein Taghavi
#UF Scripps
#07/17/22
######################


import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import umap
import matplotlib.pyplot as plt


import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import umap
import matplotlib.pyplot as plt

def smiles_to_fingerprints(smiles_list):
    """Convert a list of SMILES to RDKit fingerprints."""
    fingerprints = []
    valid_smiles = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
                arr = np.zeros((0,), dtype=np.int8)
                DataStructs.ConvertToNumpyArray(fingerprint, arr)
                fingerprints.append(arr)
                valid_smiles.append(smiles)
            else:
                print(f"Invalid SMILES string: {smiles}")
        except Exception as e:
            print(f"Error processing SMILES '{smiles}': {e}")
    return np.array(fingerprints), valid_smiles

def plot_umap(fingerprints, output_file):
    """Perform UMAP on fingerprints and plot the result."""
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(fingerprints)

    plt.figure(figsize=(10, 8))
    plt.scatter(embedding[:, 0], embedding[:, 1], s=5, cmap='Spectral', alpha=0.5)
    plt.title('UMAP projection of SMILES data')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.colorbar(label='UMAP coordinates')
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python smiles_to_umap.py <input_smiles_file> <output_image_file>")
        sys.exit(1)

    input_smiles_file = sys.argv[1]
    output_image_file = sys.argv[2]

    # Read SMILES from tab-separated file
    df = pd.read_csv(input_smiles_file, sep='\t', header=None)
    smiles_list = df.iloc[:, 0].dropna().tolist()

    # Convert SMILES to fingerprints
    fingerprints, valid_smiles = smiles_to_fingerprints(smiles_list)

    if len(fingerprints) == 0:
        print("No valid SMILES strings found. Exiting.")
        sys.exit(1)

    # Plot UMAP
    plot_umap(fingerprints, output_image_file)

