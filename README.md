Explanations:\
a series of scripts to prepare the small molecule and RNA for docking\
Python Script: The Python script add_hydrogens_and_convert.py\
takes two command-line arguments: the input directory containing PDB files and the output directory where the processed files will be saved.\
It processes each PDB file, adds hydrogens, and converts it to the PDBQT format using OpenBabel.\
SLURM Batch Script:\
The batch script submit_job.sh sets up the job parameters (such as job name, output/error logs, resources required) and runs the Python script.\
It loads the necessary modules (pymol and openbabel) and specifies the input and output directories.\
Submission:\
The job is submitted to the HPC scheduler with sbatch submit_job.sh, which schedules and runs the job on the available compute nodes.\
smiles_to_umap_hdbscan.py\

HDBSCAN Clustering:\
python smiles_to_umap_hdbscan_centers.py smiles.tsv umap_hdbscan_centers_plot.png
After obtaining the UMAP embeddings, we apply HDBSCAN clustering with hdbscan.HDBSCAN(min_cluster_size=5).\
You can adjust the min_cluster_size parameter based on your specific dataset and clustering needs.\
Plotting with Clusters:\

The plot_umap_with_clusters function generates a scatter plot of the UMAP-reduced data,\
coloring the points based on their cluster assignments from HDBSCAN.\
Each unique cluster will be represented by a different color.\
Input and Output:\

The script takes a tab-separated file with SMILES strings in the first column and outputs a PNG file of the UMAP plot with HDBSCAN clustering.\
Running the Script:\
Prepare your SMILES file:\
Create a tab-separated text file (e.g., smiles.tsv) where the first column contains the SMILES strings.\
python smiles_to_umap_hdbscan.py smiles.tsv umap_hdbscan_plot.png\

smiles_to_umap_hdbscan_centers.py\
Explanations:\
Store Molecule Objects:\
The smiles_to_fingerprints function now returns both the fingerprints and the corresponding RDKit molecule objects.\
Identify and Plot Cluster Centers:\
In plot_umap_with_clusters, after performing UMAP and HDBSCAN, the script identifies cluster centers by averaging the points in each cluster.\
For each cluster, the center coordinates and a representative molecule are identified.\
The molecule images are overlaid on the UMAP plot at the cluster center coordinates.\
Overlay Molecule Images:\
RDKit's Draw.MolToImage function is used to generate images of the molecules, which are then overlaid on the UMAP plot using plt.imshow.\
Running the Script\
Prepare your SMILES file:\
Create a tab-separated text file (e.g., smiles.tsv) where the first column contains the SMILES strings.\
Run the script:\
python smiles_to_umap_hdbscan_centers.py smiles.tsv umap_hdbscan_centers_plot.png\

Script: umap_kmeans.py\
Replace your_input_file.tsv with the path to your tab-separated file containing SMILES strings, and 5 with the number of clusters you want for the K-means algorithm.

Explanation
Reading SMILES: The script reads SMILES strings from a tab-separated input file.\
Converting to Molecules: The SMILES strings are converted to RDKit molecule objects.\
Computing UMAP Embedding: UMAP is used to reduce the dimensionality of the molecular fingerprints.\
Applying K-means: The UMAP embeddings are clustered using K-means.\
python umap_kmeans.py your_input_file.tsv [# of clusters]\


