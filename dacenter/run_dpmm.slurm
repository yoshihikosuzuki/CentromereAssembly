#!/bin/bash
#SBATCH -J run_dpmm
#SBATCH -o sbatch_stdout
#SBATCH -e sbatch_stderr
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=bigmem

python dpmm.py --vmatrix_fname /home/ysuzuki/work/RepeatAssembly/datruf_run/dmel/backup/peak/peak.3.unit.1.seqs.V --n_iteration 100000 --n_parallel 48
