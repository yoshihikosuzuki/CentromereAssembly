#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

# Load config file
. $1

rm -f .${DB_PREFIX}.*.tan.* .${DB_PREFIX}.tan.* TAN.*
mkdir -p datander
HPC.TANmask -T${N_CORE_DATANDER} ${DB_PREFIX}.db > datander/run_datander.sh

if ${USE_JOB_SCHEDULER}; then
    echo "Starting datander"
    python -m BITS.submit_job datander/run_datander.sh ${JOB_SCHEDULER} ${SUBMIT_JOB} \
           job_name="run_datander" \
           out_log="datander/datander_log.stdout" \
           err_log="datander/datander_log.stderr" \
           n_core=${N_CORE_DATANDER} \
           $([ -n "${QUEUE_OR_PARTITION}" ] && echo "queue_or_partition=${QUEUE_OR_PARTITION}") \
           $([ -n "${TIME_LIMIT}" ] && echo "time_limit=${TIME_LIMIT}") \
           $([ -n "${MEM_LIMIT}" ] && echo "mem_limit=${MEM_LIMIT}") \
           wait=True
    echo "Finished datander"
else
    bash datander/run_datander.sh > datander/datander.log 2>&1
fi

Catrack -v ${DB_PREFIX} tan
rm .${DB_PREFIX}.*.tan.*
