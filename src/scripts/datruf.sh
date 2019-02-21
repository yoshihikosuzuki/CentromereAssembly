#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

# Load config file
. $1

if ${USE_JOB_SCHEDULER}; then
    datruf_run_distribute.py -n ${N_CORE_DATRUF} \
                             -p ${N_DISTRIBUTE} \
                             -j ${JOB_SCHEDULER} \
                             -c ${SUBMIT_JOB} \
                             $([ -n "${QUEUE_OR_PARTITION}" ] && echo "-q ${QUEUE_OR_PARTITION}") \
                             $([ -n "${DEBUG_MODE}" ] && echo "-D") \
                             ${DB_PREFIX}.db \
                             TAN.${DB_PREFIX}.las
else
    datruf_run.py -n ${N_CORE_DATRUF} \
                  $([ -n "${DEBUG_MODE}" ] && echo "-D") \
                  ${DB_PREFIX}.db \
                  TAN.${DB_PREFIX}.las \
                  > datruf.log 2>&1
fi
