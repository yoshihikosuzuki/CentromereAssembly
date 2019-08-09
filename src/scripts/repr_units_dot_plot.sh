#!/bin/bash
# NOTE: Use python2 for this script

REF_UNITS=$HOME/work2/project/CentromereAssembly/assets/dmel_reference_units.fasta
REPR_UNITS=$HOME/work2/project/CentromereAssembly/develop/dmel_small/repr_units
FLEXI_DOT=$HOME/work2/software/flexidot/code/flexidot_v1.05.py

awk -F'\t' 'NR > 1 {print ">" $1 "/" $2 "/0_" $5 "\n" $6}' ${REPR_UNITS} > ${REPR_UNITS}.fasta

python2 ${FLEXI_DOT} -i ${REF_UNITS},${REPR_UNITS}.fasta -p 2 -x 1 -X -f 1
