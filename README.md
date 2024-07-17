Explanations:\
a series of scripts to prepare the small molecule and RNA for docking
Python Script: The Python script add_hydrogens_and_convert.py takes two command-line arguments: the input directory containing PDB files and the output directory where the processed files will be saved.
It processes each PDB file, adds hydrogens, and converts it to the PDBQT format using OpenBabel.
SLURM Batch Script: The batch script submit_job.sh sets up the job parameters (such as job name, output/error logs, resources required) and runs the Python script.
It loads the necessary modules (pymol and openbabel) and specifies the input and output directories.
Submission: The job is submitted to the HPC scheduler with sbatch submit_job.sh, which schedules and runs the job on the available compute nodes.
