#!/bin/bash
#SBATCH -J run_dpmm
#SBATCH -o sbatch_stdout
#SBATCH -e sbatch_stderr
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=batch

python dpmm.py
