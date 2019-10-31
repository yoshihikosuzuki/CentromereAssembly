#!/bin/bash

N_ARGS=2
if [ $# -ne ${N_ARGS} ]; then
  echo "Usage: bash fasta_to_db.sh <in_fasta> <out_db_prefix>"
  exit 1
fi

IN_FASTA=$1
OUT_DB_PREFIX=$2

fasta2DB ${OUT_DB_PREFIX} ${IN_FASTA}
DBsplit -s500 ${OUT_DB_PREFIX}
