#!/bin/bash
#SBATCH --job-name=add_hydrogens
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4  # Adjust based on the number of cores you want to use
#SBATCH --time=01:00:00
#SBATCH --mem=2G

# Load the conda environment
source /path/to/conda.sh  # Adjust the path to where your conda.sh is located
conda activate pymol_env  # Replace with the name of your conda environment

# Define directories
INPUT_DIR="/path/to/your/pdb_files"
OUTPUT_DIR="/path/to/save/processed_files"
NUM_WORKERS=$SLURM_CPUS_PER_TASK

# Run the master Python script
python process_pdb_files.py $INPUT_DIR $OUTPUT_DIR $NUM_WORKERS

# Deactivate the conda environment
conda deactivate

