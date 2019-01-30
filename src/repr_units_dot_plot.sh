#!/bin/bash

REF_UNITS=/home/yoshihiko_s/work2/project/CentromereAssembly/assets/dmel_reference_units.fasta
REPR_UNITS=/home/yoshihiko_s/work2/project/CentromereAssembly/develop/dmel_small/repr_units
FLEXI_DOT=/home/yoshihiko_s/work2/software/flexidot/code/flexidot_v1.05.py

awk -F'\t' 'NR > 1 {print ">" $2 "/" $3 "/0_" $6 "\n" $7}' ${REPR_UNITS} > ${REPR_UNITS}.fasta

source activate py27
python ${FLEXI_DOT} -i ${REF_UNITS},${REPR_UNITS}.fasta -p 2 -x 1 -X -f 1
conda deactivate
