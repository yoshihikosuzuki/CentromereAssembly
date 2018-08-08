## User parameters ##

# sge or slurm
JOB_SCHEDULER = slurm

# qsub or sbatch
SUBMIT_JOB = sbatch

ORIGINAL_DB_DIR = /projects/dazzler/runs/dmel/
DB_PREFIX = DMEL

# for datander
N_CORE_DATANDER = 16
MEM_DATANDER = 5000

# for datruf
N_DISTRIBUTE = 32
N_CORE_DATRUF = 4

#####################



## Copy dazzler db files

DB_OBJ = $(DB_PREFIX).db .$(DB_PREFIX).bps .$(DB_PREFIX).idx .$(DB_PREFIX).qvs

# Macro for copying a file
define COPY
$(1): $(2)
	cp $(2) $(1)
endef

# Define rules for each file to be copied
$(foreach FILE, $(DB_OBJ), $(eval $(call COPY, $(FILE), $(addprefix $(ORIGINAL_DB_DIR), $(FILE)))))


## Run datander

TAN.$(DB_PREFIX).las: $(DB_PREFIX).db
	HPC.TANmask -T$(N_CORE_DATANDER) $^ > run_datander.sh
	python -m BITS.$(JOB_SCHEDULER)_nize run_datander.sh job_name="run_datander" n_core=$(N_CORE_DATANDER) mem_per_cpu=$(MEM_DATANDER)
	$(SUBMIT_JOB) run_datander.sh.$(JOB_SCHEDULER)
	Catrack -v $(DB_PREFIX) tan
	rm .$(DB_PREFIX).*.tan.*


## Run datruf

datruf_result datruf_units.fasta: $(DB_PREFIX).db TAN.$(DB_PREFIX).las
	datruf_run_distribute.py -n $(N_CORE_DATRUF) -p $(N_DISTRIBUTE) -j $(JOB_SCHEDULER) -c $(SUBMIT_JOB) $^


## Run dacmaster

peaks.pkl: datruf_units.fasta
	sbatch damaster_run.py $^


## Targets

datruf: $(DB_OBJ) TAN.$(DB_PREFIX).las datruf_units.fasta
damaster: peaks.pkl

clean_datruf:
	rm datruf_* run_datruf.* finalize_datruf.sh sbatch*
