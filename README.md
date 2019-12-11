# CentromereAssembly

A bundle of modules for centromere assembly only with PacBio long reads.

Jupyter Notebooks used for Drosophila CCS reads are available [here](https://mlab.cb.k.u-tokyo.ac.jp/~yoshihiko_s/jupyter_nbs_for_submission_1210.zip). Note that our implementation currently depends on Consed, an unpublished work by Dr. Gene Myers, to compute consensus sequences, and thus one cannot directly run these codes. We will prepare alternative codes for it soon.

## Requirements

* Python3
* [`DAZZ_DB`](https://github.com/thegenemyers/DAZZ_DB) and [`DALIGNER`](https://github.com/thegenemyers/DALIGNER)
* [`DAMASKER`](https://github.com/yoshihikosuzuki/DAMASKER) (submodule of this Git repository)
* Python package [`cython`](https://cython.readthedocs.io/en/latest/src/quickstart/install.html)
   * Other Python packages required will be automatically installed

## Installation

After you satisfy all the requirements above, the typical installation by `setuptools` is available:

```bash
$ python setup.py install
```

VCA offers two ways of executions: 1) via command-line or 2) as python modules (I recommend using Jupyter Notebook).

## Command-line execution

The command to run VCA is:

```bash
usage: vca [-h] [-c CONFIG_FNAME] [task_name]

VCA: Vertebrate Centromere Assembler.

positional arguments:
  task_name             Task name. This must be one of {'all', 'datruf',
                        'dacmaster', 'dalayout'}. [all]

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG_FNAME, --config_fname CONFIG_FNAME
                        Config file name. [config]
```

### DAZZ_DB file

As input data, a `.db` file of your sequence data is required. See [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) for details.

### Config file

To prepare the config file required by VCA, first copy `config.template` in the root directory of this repository while renaming as `config`. The template is written in TOML format as follows:

```ini
# Config file must be TOML-formatted [https://github.com/toml-lang/toml].

db_prefix = "DMEL"   # <db_prefix>.db and its related files must exist in the execution firectory

[job_scheduler]

enabled = false   # Use a job scheduler if true

  [job_scheduler.params]

  scheduler_name = "sge"   # Only "sge" or "slurm"; Name of job scheduler
  submit_command = "qsub"  # e.g. "qsub" for SGE and "sbatch" for SLURM
  #queue_name     =        # You can specify the queue name (for SGE) or partition name (for SLURM)

...
```

And then edit `config` according to your data and environment.

### What is `task_name`?

The `all` mode (default) will execute all the tasks of VCA, and the other task names are mainly for development and debug.

## Jupyter Notebook execution

Instead of command-line execution, you can import functions/classes of VCA submodules and run them inside REPL or Jupyter Notebook. Every module of VCA can be used separatedly in a more exploratory and customizable manner through some visualizations. Details of each VCA submodule along with Jupyter Notebook examples are described below.

### datander ([Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/notebooks/usage/1.%20datander.ipynb))

Datander has been developed by Dr. Gene Myers as a program in his module named [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER).

### datruf ([Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/notebooks/usage/2.%20%20datruf.ipynb))

Datruf detects tandem repeat units with the output of datander. Minimum required unit length to be accurately inferred by datruf is approximately >50 bp due to the resolution of long-read alignment on which datander relies. The size of centromeric tandem repeat units is, however, known to be generally longer than it (e.g. ~120 bp and ~360 bp in *Drosophila*).

### ReadViewer ([Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/notebooks/usage/3.%20ReadViewer.ipynb))

This submodule offers a class for visualizing reads with tandem repeats found by datander and units by datruf.
