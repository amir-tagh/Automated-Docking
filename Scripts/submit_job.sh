#!/bin/bash

#SBATCH --job-name=add_hydrogens
#SBATCH -n 1
##SBATCH --partition=gpu
##SBATCH --nodes=1 
##SBATCH --ntasks=1
##SBATCH --gres=gpu:1
#SBATCH --output=RNA-error.out
##SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
##SBATCH --nodelist=c1007a-s11
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept-b
#SBATCH --mem=20G

# Load necessary modules
module purge
python/3.10
module load openbabel/3.1.1



# Define directories
INPUT_DIR="/path/to/your/pdb_files"
OUTPUT_DIR="/path/to/save/processed_files"

# Run the Python script
python add_hydrogens_and_convert.py $INPUT_DIR $OUTPUT_DIR

