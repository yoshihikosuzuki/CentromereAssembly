import numpy as np
from BITS.util.scheduler import Scheduler


if __name__ == "__main__":
    n_distribute = 20
    n_core = 10

    s = Scheduler("sge", "qsub", "all.q")

    for i in range(n_distribute):
        index = str(i + 1).zfill(int(np.log10(n_distribute) + 1))
        out_fname = f"ovlps.{index}.pkl"
        script = (f"python $HOME/work2/project/CentromereAssembly/src/ava_unsync_reads_overlap.py "
                  f"{i} {n_distribute} {n_core} {out_fname}")
        script_fname = f"run_ava.sh.{index}"
        
        s.submit(script,
                 script_fname,
                 job_name="ava_ovlp",
                 n_core=n_core)
