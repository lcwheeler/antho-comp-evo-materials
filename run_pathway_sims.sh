#!/bin/bash
#SBATCH --nodes=1                    		# Number of requested nodes
#SBATCH --time=23:30:00               		# Max wall time
#SBATCH --partition=shas             		# Specify Summit haswell nodes
#SBATCH --ntasks=24          	 	        # Number of tasks per job
#SBATCH --job-name=simulations	         	# Job submission name
#SBATCH --output=simulations.%j.out	     	# Output file name with Job ID


# Written by:	Lucas Wheeler
# Date:		21 December 2018
# Updated:	21 December 2018
# Purpose: 	This script runs my custom Python3 pathway simulation script

# purge all existing modules
module purge

# The directory where you want the job to run
cd /scratch/summit/luwh7529/pathway_simulations/all_mutations/reformulated_sims/del_90_percent

# Run simulation script. Log start and end times.
date
python3 pathway_sims_final.py
date
