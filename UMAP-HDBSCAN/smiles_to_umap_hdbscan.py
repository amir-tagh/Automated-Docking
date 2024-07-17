##########################################################################################
#HDBSCAN Clustering:
#After obtaining the UMAP embeddings, apply HDBSCAN clustering with hdbscan.\
#HDBSCAN(min_cluster_size=5). You can adjust the min_cluster_size parameter.\
#Plotting with Clusters:
#The plot_umap_with_clusters function generates a scatter plot of the UMAP-reduced data.
#The script takes a tab-separated file with SMILES strings in the first column.
#Run the script:
#python smiles_to_umap_hdbscan.py smiles.tsv umap_hdbscan_plot.png
#########################################################################################



import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import umap
import hdbscan
import matplotlib.pyplot as plt

def smiles_to_fingerprints(smiles_list):
    """Convert a list of SMILES to RDKit fingerprints."""
    fingerprints = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
                arr = np.zeros((0,), dtype=np.int8)
                DataStructs.ConvertToNumpyArray(fingerprint, arr)
                fingerprints.append(arr)
            else:
                print(f"Invalid SMILES string: {smiles}")
        except Exception as e:
            print(f"Error processing SMILES '{smiles}': {e}")
    return np.array(fingerprints)

def plot_umap_with_clusters(fingerprints, output_file):
    """Perform UMAP on fingerprints, apply HDBSCAN, and plot the result."""
    # Perform UMAP
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(fingerprints)

    # Apply HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(min_cluster_size=5)
    cluster_labels = clusterer.fit_predict(embedding)

    # Create a scatter plot of UMAP results colored by cluster labels
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(embedding[:, 0], embedding[:, 1], c=cluster_labels, cmap='Spectral', s=20, alpha=0.5)
    plt.title('UMAP projection with HDBSCAN Clustering')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.colorbar(scatter, label='Cluster Labels')
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python smiles_to_umap_hdbscan.py <input_smiles_file> <output_image_file>")
        sys.exit(1)

    input_smiles_file = sys.argv[1]
    output_image_file = sys.argv[2]

    # Read SMILES from tab-separated file
    df = pd.read_csv(input_smiles_file, sep='\t', header=None)
    smiles_list = df.iloc[:, 0].dropna().tolist()

    # Convert SMILES to fingerprints
    fingerprints = smiles_to_fingerprints(smiles_list)

    if len(fingerprints) == 0:
        print("No valid SMILES strings found. Exiting.")
        sys.exit(1)

    # Plot UMAP with HDBSCAN Clustering
    plot_umap_with_clusters(fingerprints, output_image_file)

