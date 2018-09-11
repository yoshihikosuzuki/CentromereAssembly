#!/bin/bash

#MASTER_UNITS=./master_units.fasta
#MASTER_UNITS=./master_orig.fasta
MASTER_UNITS=./homo_units.fasta
#READS=./reads.fasta
READS=./reads.peak3.fasta
OUT=./homo.blastn

source ~/.bash_profile
#makeblastdb -dbtype nucl -in ${READS} -parse_seqids
#blastn -num_threads 12 -db ${READS} -query ${MASTER_UNITS} -out master_mapping.blastn -word_size 7 -qcov_hsp_perc 60 -outfmt "17 SQ" -parse_deflines
blastn -num_threads 12 -db ${READS} -query ${MASTER_UNITS} -out ${OUT} -word_size 7 -qcov_hsp_perc 60 -outfmt "6 qseqid sseqid pident length qstart qend sstart send" -parse_deflines
