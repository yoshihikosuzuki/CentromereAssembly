#!/bin/bash
#SBATCH -J run_datruf
#SBATCH -o sbatch_stdout
#SBATCH -e sbatch_stderr
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=batch

EXEC="/home/ysuzuki/work/RepeatAssembly/datruf/datruf_run.py"
DB_FILE="DMEL.db"
LAS_FILE="TAN.DMEL.las"
END_ID=100000

python ${EXEC} ${DB_FILE} ${LAS_FILE} --start_dbid 1 --end_dbid ${END_ID} > out
