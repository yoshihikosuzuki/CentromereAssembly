# CentromereAssembly

A bundle of modules for centromere assembly only with PacBio long reads.

## Requirements

* Python version must be 3
* [`DAZZ_DB`](https://github.com/thegenemyers/DAZZ_DB) and [`DALIGNER`](https://github.com/thegenemyers/DALIGNER)
* Git submodule [`DAMASKER`](https://github.com/yoshihikosuzuki/DAMASKER) in this repository
* Python library [`cython`](https://cython.readthedocs.io/en/latest/src/quickstart/install.html)
* (Other Python libraries required are automatically installed during the installation described below)

## Installation

After you satisfy all the requirements above, the typical installation by `setuptools` is available:

```bash
$ python setup.py install
```

## Command usage

(If you want to understand dacembler's overall workflow before running codes, please read the next section firstly.)

Briefly, all you have to do is run the installed commands in the order listed below. You can instead use a `Makefile` which performs the same series of commands, after specifying some user-dependent parameters located at the top of the file.

```
usage: datruf_run.py [-h] [-d DBDUMP] [-l LADUMP] [-s START_DBID]
                     [-e END_DBID] [-m OUT_MAIN_FNAME] [-u OUT_UNITS_FNAME]
                     [--only_interval] [--on_the_fly] [-n N_CORE] [-D]
                     db_file las_file

Detect tandem repeat intervals and their unit sequences in PacBio reads.

positional arguments:
  db_file               DAZZ_DB file
  las_file              output of TANmask of the modified DAMASTER package

optional arguments:
  -h, --help            show this help message and exit
  -d DBDUMP, --dbdump DBDUMP
                        Output of `DBdump -r -h -mtan <db_file>`. This will be
                        automatically generated if it does not exist. In
                        <on_the_fly> mode, this is not used. [datander_ladump]
  -s START_DBID, --start_dbid START_DBID
                        Start read ID, which is used in DAZZ_DB. Set <= 1 to
                        start from the first read. [1]
  -e END_DBID, --end_dbid END_DBID
                        End read ID. Set < 1 to end at the last read. [-1]
  -m OUT_MAIN_FNAME, --out_main_fname OUT_MAIN_FNAME
                        Write main results to this file. [datruf_result]
  -u OUT_UNITS_FNAME, --out_units_fname OUT_UNITS_FNAME
                        Write unit sequences to this file. [datruf_units]
  --only_interval       Stop calculation just after obtaining TR intervals.
                        Since filtering of TRs by CV of its unit lengths is
                        not applied, (intervals of) TRs with short (<50 bp)
                        units will be output, unlike the default mode. [False]
  --on_the_fly          Generate dump data for each read on the fly. This mode
                        is very slow and used only when whole data are huge
                        and you just want to look at results of only several
                        reads. [False]
  -n N_CORE, --n_core N_CORE
                        Degree of parallelization. [1]
  -D, --debug_mode      Run in debug mode. [False]
```

Or you can use the distributed version instead of `datruf_run.py` (You need to run a script after all the distributed jobs finish as instructed by a logger message of the script):

```
usage: datruf_run_distribute.py [-h] [-d DBDUMP] [-l LADUMP] [-s START_DBID]
                                [-e END_DBID] [-m OUT_MAIN_FNAME]
                                [-u OUT_UNITS_FNAME] [--only_interval]
                                [--on_the_fly] [-n N_CORE] [-D]
                                [-p N_DISTRIBUTE] [-j JOB_SCHEDULER]
                                [-c SUBMIT_JOB]
                                db_file las_file

Distribute datruf_run.py jobs using a job scheduler.

positional arguments:
  db_file               DAZZ_DB file
  las_file              output of TANmask of the modified DAMASTER package

optional arguments:
  -h, --help            show this help message and exit
  -d DBDUMP, --dbdump DBDUMP
                        Output of `DBdump -r -h -mtan <db_file>`. This will be
                        automatically generated if not exist.
                        [datander_dbdump]
  -l LADUMP, --ladump LADUMP
                        Output of `LAdump -c <db_file> <las_file>`. This will
                        be automatically generated if not exist.
                        [datander_ladump]
  -s START_DBID, --start_dbid START_DBID
                        Start read ID, which is used in DAZZ_DB. Set <= 1 to
                        start from the first read. [1]
  -e END_DBID, --end_dbid END_DBID
                        End read ID. Set < 1 to end at the last read. [-1]
  -m OUT_MAIN_FNAME, --out_main_fname OUT_MAIN_FNAME
                        Write main results to this file. [datruf_result]
  -u OUT_UNITS_FNAME, --out_units_fname OUT_UNITS_FNAME
                        Write unit sequences to this file. [datruf_units]
  --only_interval       Stop calculation just after obtaining TR intervals.
                        Since filtering of TRs by CV of its unit lengths is
                        not applied, (intervals of) TRs with short (<50 bp)
                        units will be output, unlike the default mode. [False]
  --on_the_fly          Generate dump data for each read on the fly. This mode
                        is very slow and used only when whole data are huge
                        and you just want to look at results of only several
                        reads. [False]
  -n N_CORE, --n_core N_CORE
                        Degree of parallelization in each distributed job. [1]
  -D, --debug_mode      Run datruf in debug mode. [False]
  -p N_DISTRIBUTE, --n_distribute N_DISTRIBUTE
                        Degree of parallelization in each distributed job. [1]
  -j JOB_SCHEDULER, --job_scheduler JOB_SCHEDULER
                        Job scheduler name. ('sge' or 'slurm)' [sge]
  -c SUBMIT_JOB, --submit_job SUBMIT_JOB
                        Command name to submit a job with the specified
                        scheduler. [qsub]
```

As for the input files of `datuf_run.py` or `datruf_run_distribute.py`:

* `db_file` is a `.db` file generated from the PacBio reads by `DAZZ_DB`. Be careful that all the other related files like `.*.bps` and `.*.idx` must also exist in the same directory.
* `las_file` is a `TAN.*.las` file generated by `HPC.TANmask` of modified `DAMASKER`.

```
usage: dacmaster_run.py [-h] [-u UNITS_FNAME] [-m MIN_N_UNITS] [-n N_CORE]
                        [-F] [-D]
                        db_file

Construct master units and representative units.

positional arguments:
  db_file               DAZZ_DB file. Used for generating reads fasta and
                        dbid_header.

optional arguments:
  -h, --help            show this help message and exit
  -u UNITS_FNAME, --units_fname UNITS_FNAME
                        Input file of the unit sequences reported by datruf.
                        [datruf_units]
  -m MIN_N_UNITS, --min_n_units MIN_N_UNITS
                        Minimum number of units in a single TR whose consensus
                        will be taken. [10]
  -n N_CORE, --n_core N_CORE
                        Degree of parallelization. [1]
  -F, --from_scratch    Compute from scratch without loading existing
                        peaks.pkl file. [False]
  -D, --debug_mode      Run in debug mode. [False]
```

The input file `db_file` of `dacmaster_run.py` is the same as that of `datruf_run.py`.

```
usage: dalayout_run.py [-h] [-r TR_READS_FNAME] [-u REPR_UNITS_FNAME]
                       [-p PEAKS_FINDER_FNAME] [-t VARIANT_FRACTION]
                       [-o OUT_PKL_FNAME] [-n N_CORE] [-D]

Perform layout of reads based on the representative units given.

optional arguments:
  -h, --help            show this help message and exit
  -r TR_READS_FNAME, --tr_reads_fname TR_READS_FNAME
                        A file of the TR reads. [tr_reads]
  -u REPR_UNITS_FNAME, --repr_units_fname REPR_UNITS_FNAME
                        A file of the representative units generated by
                        dacmaster. [repr_units]
  -p PEAKS_FINDER_FNAME, --peaks_finder_fname PEAKS_FINDER_FNAME
                        PeaksFinder pickle file which is dacmaster's output.
                        [peaks_finder.pkl]
  -t VARIANT_FRACTION, --variant_fraction VARIANT_FRACTION
                        Value of Consed's -t option. 0.0 means all variant
                        sites will be used. [0.0]
  -o OUT_PKL_FNAME, --out_pkl_fname OUT_PKL_FNAME
                        Output pickle file for encodings. [encodings.pkl]
  -n N_CORE, --n_core N_CORE
                        Degree of parallelization. [1]
  -D, --debug_mode      Run in debug mode. [False]
```


## Overview of the modules

(TBA: a picture of the entire workflow of dacembler)

Every module offers some visualization functions (examples will be shown here).

## datander

Datander has been developed by Gene Myers as a program in his module named [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER).

## datruf

detects tandem repeat regions (with datander by Gene Myers) and tandem repeat units

Datruf detects tandem repeats and unit sequences within them only from PacBio raw reads and without any other external information. Before using datruf, you must run the `HPC.TANmask` program of [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER) originally developed by Gene Myers and slightly modified for datruf, and must have `.db` file and `.las` file.  The minimum length of tandem repeat units accurately infered by datruf is relatively long (approximately >50 bp) due to the resolution of long-read alignment on which datander relies. However, datruf is developed mainly for (core) centromere assembly, whose unit size is known to be long in general (e.g. 171 bp alpha-satellite sequence of human).

Datruf consists of three parts:

* **Runner** applies core algorithm of datruf to arbitrary number of reads without visualization and output the results (read ID, start/end position of tandem repeat, estimated unit length, borders of units, consensus unit sequence) into a file.

* **Viewer** receives one read ID and draw figures that are not told by Runner but help to understand the sequence structure of the read.

* **Plotter** requires results of datruf and/or other software for tandem repeat analysis (TRF, mTR, etc.) and draws figures useful for comparison of these results.

See [Jupyter notebook]() for how to use.

## dacmaster

* infers centromeric satellite repeat units (monomers) from all the tandem repeat units
* performs clustering of the monomers at multiple, stage-specific resolusions
* assigns each monomer in the reads to a representative (consensus) of the clusters

Dacmaster recieves the fasta file of the units reported by datruf, and determine a set of representative centromeric monomers from it. This task is done by

1. Detecting peaks in tandem repeat unit length (which are considered as those of centromeric monomers),
2. Collecting units around each peak
3. Clustering units which are the first units of the first tandem repeats not starting from the end point of each read
4. Adjusting start positions of all monomers by the most similar representative monomer

See [Jupyter notebook]().

## dalayout

layouts the reads encoded with the representative monomers
