## Copy dazzler db files

ORIGINAL_DMEL_DB_DIR = /projects/dazzler/runs/dmel/
DMEL_DB_OBJ = DMEL.db .DMEL.bps .DMEL.idx .DMEL.qvs

# Macro for copying a file
define COPY
$(1): $(2)
	cp $(2) $(1)
endef

# Define rules for each file to be copied
$(foreach FILE, $(DMEL_DB_OBJ), $(eval $(call COPY, $(FILE), $(addprefix $(ORIGINAL_DMEL_DB_DIR), $(FILE)))))

## Run datander

#TAN.DMEL.las: DMEL.db
#       sbatch hpc_tanmask.slurm


## Targets

.PHONY: all
all: $(DMEL_DB_OBJ)
