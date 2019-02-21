#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

# Load config file
. $1

echo "dalayout_run.py -n ${N_CORE_DACMASTER} $([ -n "${DEBUG_MODE}" ] && echo "-D")" > run_dalayout.sh

if ${USE_JOB_SCHEDULER}; then
    echo "Starting dalayout"
    python -m BITS.submit_job run_dalayout.sh ${JOB_SCHEDULER} ${SUBMIT_JOB} \
           job_name="run_dalayout" \
           out_log="dalayout_log.stdout" \
           err_log="dalayout_log.stderr" \
           n_core=${N_CORE_DATANDER} \
           $([ -n "${QUEUE_OR_PARTITION}" ] && echo "queue_or_partition=${QUEUE_OR_PARTITION}") \
           $([ -n "${TIME_LIMIT}" ] && echo "time_limit=${TIME_LIMIT}") \
           $([ -n "${MEM_LIMIT}" ] && echo "mem_limit=${MEM_LIMIT}") \
           wait=True
    echo "Finished dalayout"
else
    bash run_dalayout.sh > dalayout.log 2>&1
fi

# For now only sge
dalayout_run_graph.py -p ${N_DISTRIBUTE_DALAYOUT} -n ${N_CORE_DALAYOUT}
# After finishing all jobs
#dalayout_gather_ava_sub.py
