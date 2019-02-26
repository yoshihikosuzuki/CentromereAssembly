#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

# Load config file
. $1

mkdir -p dalayout
echo "dalayout_run_encode.py -n ${N_CORE_ENCODE} $([ -n "${N_INDEX}" ] && echo "-i ${N_INDEX}") $([ -n "${VARIANT_FRAC}" ] && echo "-t ${VARIANT_FRAC}") $([ -n "${VCALL_HC}" ] && echo "-H") $([ -n "${DEBUG_MODE}" ] && echo "-D")" > dalayout/run_dalayout_encode.sh

if ${USE_JOB_SCHEDULER}; then
    echo "Starting encoding"
    python -m BITS.submit_job dalayout/run_dalayout_encode.sh ${JOB_SCHEDULER} ${SUBMIT_JOB} \
           job_name="run_dalayout" \
           out_log="dalayout/dalayout_log.stdout" \
           err_log="dalayout/dalayout_log.stderr" \
           n_core=${N_CORE_ENCODE} \
           $([ -n "${QUEUE_OR_PARTITION}" ] && echo "queue_or_partition=${QUEUE_OR_PARTITION}") \
           $([ -n "${TIME_LIMIT}" ] && echo "time_limit=${TIME_LIMIT}") \
           $([ -n "${MEM_LIMIT}" ] && echo "mem_limit=${MEM_LIMIT}") \
           wait=True
    echo "Finished encoding"
else
    bash dalayout/run_dalayout_encode.sh > dalayout/dalayout.log 2>&1
fi
