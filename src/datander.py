from logzero import logger
from BITS.util.proc import run_command

dir_name = "datander"
script_fname = f"{dir_name}/run_datander.sh"
log_fname = f"{dir_name}/log"


def calc_n_blocks(db_fname):
    """Extract the number of blocks from the db file."""
    with open(db_fname, "r") as f:
        for line in f:
            if line.startswith("blocks"):
                return int(line.split("=")[1].strip())
    logger.error(f"No information on the number of blocks in {db_fname}")


def run_datander(db_prefix, n_core=1, scheduler=None):
    run_command(f"rm -f .{db_prefix}.*.tan.* .{db_prefix}.tan.* TAN.*")
    run_command(f"mkdir -p {dir_name}; rm -f {dir_name}/*")

    script = run_command(f"HPC.TANmask -T{n_core} {db_prefix}.db")
    if calc_n_blocks(f"{db_prefix}.db") > 1:
        script += f"Catrack -v {db_prefix} tan" + "\n" + f"rm .{db_prefix}.*.tan.*"

    if scheduler is None:
        with open(script_fname, "w") as f:
            f.write(f"{script}\n")
        run_command(f"bash {script_fname} > {log_fname} 2>&1")
    else:
        scheduler.submit(script,
                         script_fname,
                         job_name=dir_name,
                         log_fname=log_fname,
                         n_core=n_core,
                         wait=True)
