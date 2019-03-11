def add_scheduler_args(p):
    p.add_argument("--job_scheduler",
                   type=str,
                   default=None,
                   help="Job scheduler name. ('sge' or 'slurm)' [None]")

    p.add_argument("--submit_command",
                   type=str,
                   default=None,
                   help="Command name to submit a job with the specified scheduler. [None]")

    p.add_argument("--queue_name",
                   type=str,
                   default=None,
                   help="Name of queue (SGE) or partition (SLURM) to which jobs are submitted. [None]")
