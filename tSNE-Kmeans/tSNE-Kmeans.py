import sys
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt

def compute_tsne_embedding(data, perplexity=30, n_iter=1000):
    tsne = TSNE(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=42)
    embedding = tsne.fit_transform(data)
    return embedding

def apply_kmeans(embedding, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    kmeans.fit(embedding)
    cluster_labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_
    return cluster_labels, cluster_centers

def save_smiles_of_cluster_centers(smiles_list, cluster_labels, cluster_centers, output_filename):
    cluster_center_smiles = []

    for i in range(len(cluster_centers)):
        cluster_idx = np.where(cluster_labels == i)[0]  # Find indices of data points in this cluster
        if cluster_idx.size > 0:
            center_idx = cluster_idx[np.argmin(np.linalg.norm(cluster_centers[i] - cluster_idx[:, None], axis=1))]
            cluster_center_smiles.append(smiles_list[center_idx])

    with open(output_filename, 'w') as f:
        for i, smiles in enumerate(cluster_center_smiles):
            f.write(f"Cluster {i}: {smiles}\n")

    print(f"Saved SMILES of cluster centers to {output_filename}")

def plot_tsne_with_clusters(embedding, cluster_labels, cluster_centers, output_image_file):
    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(embedding[:, 0], embedding[:, 1], c=cluster_labels, cmap='Spectral', s=5)
    plt.colorbar(scatter, boundaries=np.arange(len(np.unique(cluster_labels)) + 1) - 0.5).set_ticks(np.arange(len(np.unique(cluster_labels))))
    
    for i, center in enumerate(cluster_centers):
        plt.text(center[0], center[1], str(i), fontsize=12, ha='center', va='center', fontweight='bold')

    plt.title('t-SNE projection of the dataset', fontsize=24)
    plt.xlabel('t-SNE1')
    plt.ylabel('t-SNE2')
    plt.savefig(output_image_file)
    plt.close()
    print(f"Saved t-SNE plot with clusters to {output_image_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_data_file> <n_clusters> <output_image_file> <output_smiles_file>")
        sys.exit(1)

    input_data_file = sys.argv[1]
    n_clusters = int(sys.argv[2])
    output_image_file = sys.argv[3]
    output_smiles_file = sys.argv[4]

    # Load SMILES strings from input file
    try:
        df = pd.read_csv(input_data_file, sep='\t', header=None)
        smiles_list = df.iloc[:, 0].dropna().tolist()
    except FileNotFoundError:
        print(f"Error: File '{input_data_file}' not found.")
        sys.exit(1)

    # Generate Morgan fingerprints from SMILES strings
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024) for mol in molecules]
    data = np.zeros((len(fingerprints), 1024))
    for i, fp in enumerate(fingerprints):
        data[i] = np.array(fp)

    # Compute t-SNE embedding
    embedding = compute_tsne_embedding(data)

    # Apply K-means clustering
    cluster_labels, cluster_centers = apply_kmeans(embedding, n_clusters)

    # Save SMILES of cluster centers to a file
    save_smiles_of_cluster_centers(smiles_list, cluster_labels, cluster_centers, output_smiles_file)

    # Plot t-SNE with clusters
    plot_tsne_with_clusters(embedding, cluster_labels, cluster_centers, output_image_file)

    print("Process completed successfully.")

