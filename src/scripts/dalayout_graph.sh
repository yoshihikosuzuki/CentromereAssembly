#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

# Load config file
. $1

mkdir -p dalayout

# For now only sge
dalayout_run_graph.py -p ${N_DISTRIBUTE_DALAYOUT} \
                      -n ${N_CORE_DALAYOUT} \
                      -j ${JOB_SCHEDULER} \
                      -c ${SUBMIT_JOB} \
                      $([ -n "${VARIANT_FRAC}" ] && echo "-v var_vec_global${VARIANT_FRAC}$([ -n "${VCALL_HC}" ] && echo "_hc")") \
                      $([ -n "${QUEUE_OR_PARTITION}" ] && echo "-q ${QUEUE_OR_PARTITION}")
