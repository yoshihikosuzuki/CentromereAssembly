import argparse

from datruf_utils import run_command


def main():
    args = load_args()

    command = "DBdump %s | awk 'NR == 1 {print $3}'" % args.db_file
    max_dbid = int(run_command(command).strip())

    if args.start_dbid < 1:
        args.start_dbid = 1
    if args.end_dbid < 1 or args.end_dbid > max_dbid:
        args.end_dbid = max_dbid

    n_dbid_part = -(-args.end_dbid // args.n_distribute)

    for i in range(args.n_distribute):
        script_fname = "run_datruf.slurm.%d" % (i + 1)
        with open(script_fname, 'w') as f:
            f.write(
"""#!/bin/bash
#SBATCH -J run_datruf
#SBATCH -o sbatch_stdout
#SBATCH -e sbatch_stderr
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c %d
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=batch

python %s --start_dbid %d --end_dbid %d --n_core %d --out_main_fname %s --out_units_fname %s %s %s
"""
                % (args.n_core,
                   args.datruf_exec,
                   args.start_dbid + i * n_dbid_part,
                   args.start_dbid + (i + 1) * n_dbid_part - 1,
                   args.n_core,
                   args.out_main_fname + "." + str(i + 1),
                   args.out_units_fname + "." + str(i + 1),
                   args.db_file,
                   args.las_file))

        command = ("sbatch %s" % (script_fname))
        run_command(command)

    print("After finishing all jobs, you have to run following commands:\n"
          "$ cat %s.* > %s\n"
          "$ cat %s.* > %s.cat\n"
          "$ awk -F'\\t' 'NR == 1 {print $0} $1 != \"\" {print $0}' %s.cat > %s\n"
          "$ rm %s.*; rm %s.*"
          % (args.out_units_fname,
             args.out_units_fname,
             args.out_main_fname,
             args.out_main_fname,
             args.out_main_fname,
             args.out_main_fname,
             args.out_units_fname,
             args.out_main_fname,))


def load_args():
    parser = argparse.ArgumentParser(description="Run datruf with many reads")

    parser.add_argument("db_file",
                        help="DAZZ_DB file")

    parser.add_argument("las_file",
                        help="output of TANmask")

    parser.add_argument("--datruf_exec",
                        default="/home/ysuzuki/work/RepeatAssembly/datruf/datruf_run.py",
                        help="Path to datruf_run.py")

    parser.add_argument("--dbdump",
                        type=str,
                        default="datander_dbdump",
                        help=("Output of `DBdump -r -h -mtan <db_file>`. "
                              "This will be automatically generated if it "
                              "does not exist. In <on_the_fly> mode, this "
                              "is not used. [datander_dbdump]"))

    parser.add_argument("--ladump",
                        type=str,
                        default="datander_ladump",
                        help=("Output of `LAdump -c <db_file> <las_file>`. "
                              "This will be automatically generated if it "
                              "does not exist. In <on_the_fly> mode, this "
                              "is not used. [datander_ladump]"))

    parser.add_argument("--start_dbid",
                        type=int,
                        default=1,
                        help=("Start read ID, which is used in DAZZ_DB. "
                              "Set <= 1 for starting from the first read. "
                              "[1]"))

    parser.add_argument("--end_dbid",
                        type=int,
                        default=-1,
                        help=("End read ID. Set < 1 for ending at the last "
                              "read. [-1]"))

    parser.add_argument("--n_core",
                        type=int,
                        default=1,
                        help=("Degree of parallelization. [1]"))

    parser.add_argument("--n_distribute",
                        type=int,
                        default=1,
                        help=("Degree of distributing computation. [1]"))

    parser.add_argument("--out_main_fname",
                        type=str,
                        default="datruf_result",
                        help=("Write main results to this file. "
                              "[datruf_result]"))

    parser.add_argument("--out_units_fname",
                        type=str,
                        default="datruf_units.fasta",
                        help=("Write unit sequences to this file. "
                              "[datruf_units.fasta]"))

    parser.add_argument("--only_interval",
                        action="store_true",
                        default=False,
                        help=("Stop calculation just after obtaining TR "
                              "intervals. Also filtering by CV of the unit  "
                              "lengths is not applied. "
                              "[False]"))

    parser.add_argument("--on_the_fly",
                        action="store_true",
                        default=False,
                        help=("Generate dump data for each read on the fly. "
                              "This mode is very slow and used only when "
                              "whole data are huge and one just wants to "
                              "look at results of only several reads. "
                              "[False]"))

    return parser.parse_args()


if __name__ == "__main__":
    main()
