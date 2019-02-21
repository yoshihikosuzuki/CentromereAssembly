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

## How to run

The simplest way is to use `dacembler.sh` with a config file which can be prepared by modifying `config.template` in this repository root. You need to specify your environment-specific parameters, the number of cores, etc., in the config file. Once you complete it (here we call it as `config_file`), you can run dacembler by:

```bash
$ dacembler.sh config_file
```

Actually `dacembler.sh` just calls `datander.sh`, `datruf.sh`, `dacmaster.sh`, and `dalayout.sh` in this order with `config_file`. And they further execute python scripts.


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
