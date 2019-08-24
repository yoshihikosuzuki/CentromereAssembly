# CentromereAssembly

A bundle of modules for centromere assembly only with PacBio long reads.



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



## How to run

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

The `all` mode (default) will execute all tasks, and the other task names are mainly for development and debug. You can also run VCA with Jupyter Notebook. For details, see the section just below.



## Modules and their usage

Every module of VCA can be used separatedly in a more exploratory and customizable manner. Here are links to Jupyter Notebooks about usage and results with *Drosophila* data for each module offered in VCA:

- [1. datander and datruf](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/ipynbs/1.%20datander%20and%20datruf.ipynb)
- [2. TRReadFilter](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/ipynbs/2.%20TRReadFilter.ipynb)
- [3. ReadViewer](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/ipynbs/3.%20ReadViewer.ipynb)

* [TODO: Links to Jupyter NBs here]





### datander

Datander has been developed by Gene Myers as a program in his module named [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER).



### datruf

Datruf detects tandem repeat units with the output of datander. Minimum required unit length to be accurately inferred by datruf is approximately >50 bp due to the resolution of long-read alignment on which datander relies. The size of centromeric tandem repeat units is, however, known to be generally longer than it (e.g. ~120 bp and ~360 bp in *Drosophila*).

