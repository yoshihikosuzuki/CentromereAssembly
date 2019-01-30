###############################################################################
#################### User/environment-specific parameters #####################
###############################################################################

## JOB SCHEDULER (COMMENT-IN SINGLE LINE FOR EACH VARIABLE)

#USE_JOB_SCHEDULER=true
USE_JOB_SCHEDULER=false
JOB_SCHEDULER=sge
#JOB_SCHEDULER=slurm
SUBMIT_JOB=qsub
#SUBMIT_JOB=sbatch

## DAZZLER DATABASE (${DB_PREFIX}.db and .${DB_PREFIX}[.bps|.idx] must exist)

DB_PREFIX=DMEL

## DATANDER

N_CORE_DATANDER=8
MEM_DATANDER=5000

## DATRUF

N_DISTRIBUTE=8
N_CORE_DATRUF=4

## DACMASTER

MIN_N_UNITS=10
N_CORE_DACMASTER=8

################################################################################
################################################################################


## Run datander

HPC.TANmask -T${N_CORE_DATANDER} ${DB_PREFIX}.db > run_datander.sh

if ${USE_JOB_SCHEDULER}; then
	python -m BITS.${JOB_SCHEDULER}_nize \
		   run_datander.sh \
		   job_name="run_datander" \
		   n_core=${N_CORE_DATANDER} \
		   wait=False
	${SUBMIT_JOB} run_datander.sh.${JOB_SCHEDULER}
else
	bash run_datander.sh > datander.log 2>&1
fi

Catrack -v ${DB_PREFIX} tan
rm .${DB_PREFIX}.*.tan.*


## Run datruf

if ${USE_JOB_SCHEDULER};then
	datruf_run_distribute.py -n ${N_CORE_DATRUF} \
							 -p ${N_DISTRIBUTE} \
							 -j ${JOB_SCHEDULER} \
							 -c ${SUBMIT_JOB} \
							 -D \
							 ${DB_PREFIX}.db \
							 TAN.${DB_PREFIX}.las
	# TODO: automated execution of finalize.sh!; otherwise, successive dacmaster will fail
else
	datruf_run.py -n ${N_CORE_DATRUF} \
				  -D \
				  ${DB_PREFIX}.db \
				  TAN.${DB_PREFIX}.las \
				  > datruf.log 2>&1
fi


## Run dacmaster

DBshow -w10000000 ${DB_PREFIX}.db \
	| awk -F'>' 'BEGIN {print "dbid\theader\tlength\tsequence"; count = 1}
                 count % 2 == 1 {header = $2}
                 count % 2 == 0 {print (count / 2) "\t" header "\t" length($1) "\t" $1} {count++}' > reads

echo "dacmaster_run.py -m ${MIN_N_UNITS} -n ${N_CORE_DACMASTER} -D -F ${DB_PREFIX}.db" > run_dacmaster.sh

if ${USE_JOB_SCHEDULER}; then
	python -m BITS.${JOB_SCHEDULER}_nize \
		   run_dacmaster.sh \
		   job_name="run_dacmaster" \
		   n_core=${N_CORE_DACMASTER} \
		   wait=False
	${SUBMIT_JOB} run_dacmaster.sh.${JOB_SCHEDULER}
else
	bash run_dacmaster.sh > dacmaster.log 2>&1
fi
