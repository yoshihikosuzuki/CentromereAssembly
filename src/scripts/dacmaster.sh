#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

# Load config file
. $1

DBshow -w10000000 ${DB_PREFIX}.db > reads.fasta

echo "dacmaster_run.py -m ${MIN_N_UNITS} -n ${N_CORE_DACMASTER} $([ -n "${DEBUG_MODE}" ] && echo "-D")" > run_dacmaster.sh

if ${USE_JOB_SCHEDULER}; then
    echo "Starting dacmaster"
    python -m BITS.submit_job run_dacmaster.sh ${JOB_SCHEDULER} ${SUBMIT_JOB} \
           job_name="run_dacmaster" \
           out_log="dacmaster_log.stdout" \
           err_log="dacmaster_log.stderr" \
           n_core=${N_CORE_DACMASTER} \
           $([ -n "${QUEUE_OR_PARTITION}" ] && echo "queue_or_partition=${QUEUE_OR_PARTITION}") \
           $([ -n "${TIME_LIMIT}" ] && echo "time_limit=${TIME_LIMIT}") \
           $([ -n "${MEM_LIMIT}" ] && echo "mem_limit=${MEM_LIMIT}") \
           wait=True
    echo "Finished dacmaster"
else
    bash run_dacmaster.sh > dacmaster.log 2>&1
fi
