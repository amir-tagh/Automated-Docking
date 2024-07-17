###########################################
#The script reads SMILES strings from a tab-separated input file.
#The SMILES strings are converted to RDKit molecule objects
#UMAP is used to reduce the dimensionality of the molecular fingerprints.
#The UMAP embeddings are clustered using K-means
#Amirhossein Taghavi 07/17/2024
##########################################

#usage
#python umap_kmeans.py your_input_file.tsv [num of clusters]


import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import umap
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def smiles_to_molecules(smiles_list):
    molecules = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecules.append(mol)
        else:
            print(f"Invalid SMILES string: {smiles}")
    return molecules

def compute_umap_embedding(molecules):
    # Compute fingerprints
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024) for mol in molecules]
    arr = np.zeros((len(molecules), 1024))
    for i, fp in enumerate(fingerprints):
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr[i])
    # Compute UMAP embedding
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(arr)
    return embedding

def apply_kmeans(embedding, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(embedding)
    cluster_labels = kmeans.labels_
    return cluster_labels

def plot_umap_with_clusters(embedding, cluster_labels):
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(embedding[:, 0], embedding[:, 1], c=cluster_labels, cmap='Spectral', s=5, alpha=0.5)
    plt.title('UMAP projection with K-means Clustering')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.colorbar(scatter, label='Cluster Labels')
    plt.show()

if __name__ == "__main__":
    # Check if input file is provided as argument
    if len(sys.argv) != 3:
        print("Usage: python umap_kmeans.py <input_smiles_file> <n_clusters>")
        sys.exit(1)

    # Read input file path from command line argument
    input_smiles_file = sys.argv[1]
    n_clusters = int(sys.argv[2])

    # Read SMILES from file
    df = pd.read_csv(input_smiles_file, sep='\t', header=None)
    smiles_list = df.iloc[:, 0].dropna().tolist()

    # Convert SMILES to RDKit molecules
    molecules = smiles_to_molecules(smiles_list)

    # Compute UMAP embedding
    embedding = compute_umap_embedding(molecules)

    # Apply K-means clustering
    cluster_labels = apply_kmeans(embedding, n_clusters)

    # Plot UMAP with clusters
    plot_umap_with_clusters(embedding, cluster_labels)

