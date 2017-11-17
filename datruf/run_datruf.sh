#!/bin/bash

#UNIT_LEN=100
#COPY_NUM=50
#FLANKING=0
#MM=0
#INS=0
#DEL=0

UNIT_LEN=$1
COPY_NUM=$2
FLANKING=$3
MM=$4
INS=$5
DEL=$6


./count_perfect/rand_seq tr.fasta unit ${UNIT_LEN} ${COPY_NUM} ${MM} ${INS} ${DEL} ${FLANKING} ${FLANKING} 1000

rm READS* .READS* TAN* datander* tanmask.sh dbid_header dotplots/* all_reads.fasta

fasta2DB READS tr.fasta
DBsplit -s500 READS
HPC.TANmask READS > tanmask.sh
sed -i 's/rm TAN.READS.las/LAmerge TAN.READS.las TAN.READS.*.las\nrm TAN.READS.*.las/' tanmask.sh
bash tanmask.sh
LAdump -c READS.db TAN.READS.las > datander_ladump
DBdump -r -h -mtan READS.db > datander_dbdump
DBshow READS > all_reads.fasta
cat all_reads.fasta | awk -F'>' 'BEGIN {counter = 1} NF > 1 {print counter "\t" $2; counter++}' > dbid_header
python /home/ysuzuki/work/RepeatAssembly/datruf/datruf_v2.py
mv datander_result datruf_${UNIT_LEN}_${COPY_NUM}_${MM}_${INS}_${DEL}_${FLANKING}_${FLANKING}_1000
