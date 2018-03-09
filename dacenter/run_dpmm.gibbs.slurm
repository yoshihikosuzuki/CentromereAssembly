#!/bin/bash
#SBATCH -J run_dpmm_gibbs
#SBATCH -o sbatch_stdout.gibbs
#SBATCH -e sbatch_stderr.gibbs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=bigmem

# First time
#python dpmm.py --vmatrix_fname ../datruf_run/dmel/backup/peak/peak.3.unit.1.seqs.V --out_cluster_pickle clustering.large.sm.pkl --sampling_type sm --n_iteration 1000

# Additional sampling
#python dpmm.py --in_cluster_pickle clustering.large.sm.pkl --out_cluster_pickle clustering.large.sm.pkl --sampling_type sm --n_iteration 5000000

# First Gibbs sampling
python dpmm.py --in_cluster_pickle clustering.large.sm.gibbs.pkl --out_cluster_pickle clustering.large.sm.gibbs.pkl --sampling_type gibbs --start_from_max_s --n_iteration 20000000

# Additional Gibbs sapmling
#python dpmm.py --in_cluster_pickle clustering.large.sm.gibbs.pkl --out_cluster_pickle clustering.large.sm.gibbs.pkl --sampling_type gibbs --n_iteration 10000000
