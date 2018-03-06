#!/bin/bash
#SBATCH -J run_dpmm
#SBATCH -o sbatch_stdout
#SBATCH -e sbatch_stderr
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=bigmem

# First time
#python dpmm.py --vmatrix_fname ../datruf_run/dmel/backup/peak/peak.3.unit.1.seqs.V --out_cluster_pickle clustering.large.sm.pkl --sampling_type sm --n_iteration 1000

# Additional sampling
python dpmm.py --in_cluster_pickle clustering.large.sm.pkl --out_cluster_pickle clustering.large.sm.pkl --sampling_type sm --n_iteration 99000
#python dpmm.py --in_cluster_pickle clustering.large.sm.pkl --out_cluster_pickle clustering.large.gibbs.pkl --sampling_type gibbs --n_iteration 1000000
