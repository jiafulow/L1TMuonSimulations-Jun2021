#!/bin/bash
#SBATCH --partition=hpg2-compute     # Partition
#SBATCH --job-name=slurm             # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Number of processes
#SBATCH --mem=3gb                    # Total memory limit
#SBATCH --time=08:00:00              # Time limit hrs:min:sec
#SBATCH --output=logs/slurm-%j.out   # Standard output log
#SBATCH --error=logs/slurm-%j.err    # Standard error log


# Remember to get a valid proxy!
#voms-proxy-init -voms cms -valid 192:00 -out ${HOME}/.x509up_u`id -u`
#export X509_USER_PROXY=${HOME}/.x509up_u`id -u`

# Remember to create the subdirectory 'logs'
#mkdir -p logs

date; echo "Stand back I'm going to try Science!"

cmsRun pset_SingleNeutrino_PU200_Phase2HLTTDRSummer20.py

date; echo "DONE"
