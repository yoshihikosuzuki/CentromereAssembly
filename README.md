# ECA: Experimental Centromere Assembler

A bundle of modules for centromere assembly with long reads. It currently supports only PacBio CCS reads, although we are planning to adapt to more noisy reads.

Jupyter Notebooks used for assembly with _Drosophila_ CCS reads are available [here](https://mlab.cb.k.u-tokyo.ac.jp/~yoshihiko_s/jupyter_nbs_for_submission_1210.zip). Note that our implementation currently depends on software named Consed, an unpublished work by Dr. Gene Myers (be careful there is a distinct program with the same name), to compute consensus sequences, and thus one cannot directly run these codes. We will prepare alternative codes for it soon.

## Requirements

- Python3
- [`cython`](https://cython.readthedocs.io/en/latest/src/quickstart/install.html)
  - Usually by `$ pip install Cython`
  - Other Python packages required will be automatically installed below

## Installation

First clone this repository including submodules:

```bash
$ git clone --recursive https://github.com/yoshihikosuzuki/CentromereAssembly
```

Change the cloned directory to an individual environment uing `virtualenv`:

```bash
$ virtualenv -p python3 CentromereAssembly
$ source CentromereAssembly/bin/activate   # enter the virtual enviroment
```

Install [`DAZZ_DB`](https://github.com/thegenemyers/DAZZ_DB) and [`DALIGNER`](https://github.com/thegenemyers/DALIGNER) by Dr. Gene Myers:

```bash
$ cd CentromereAssembly
$ git clone https://github.com/thegenemyers/DAZZ_DB; cd DAZZ_DB
$ sed 's|~/bin|../bin|' Makefile > .Makefile; mv .Makefile Makefile   # change install dir
$ make -j; make install -j; cd ..
$ git clone https://github.com/thegenemyers/DALIGNER; cd DALIGNER
$ sed 's|~/bin|../bin|' Makefile > .Makefile; mv .Makefile Makefile   # change install dir
$ make -j; make install -j; cd ..
```

Install modified [`DAMASKER`](https://github.com/yoshihikosuzuki/DAMASKER) (original repository by Dr. Gene Myers is [here](https://github.com/thegenemyers/DAMASKER)):

```bash
$ cd DAMASKER
$ make -j; make install -j; cd ..
```

You need to reload the environment (or re-enter it) to use the executables installed above:

```bash
$ source bin/activate
```

Then install required python packages:

```bash
$ python setup.py install
```

(You will not be able to install `consed_python` because it is for now a private repository as described above.)

## How to use

ECA currently supports only exploratory execution via Jupyter Notebook. That is, ECA offers a bundle of modules and functions for centromere assembly, and one will assemble their focal tandem repeat sequences interactively while looking at some plots. We are planning to implement an integrated single command executable from Terminal in the future.

### Input file

The required input is a [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) file (`.db` file) of PacBio CCS reads. One can convert to it from different file formats like FASTA.

### Workflow overview

1. Detect tandem repeats and tandem repeat units from the reads.
1. Pick up tandem repeats one wishs to assemble, and extract reads that contain such tandem repeats.
1. Compute overlaps between the filtered reads naively by original sequence identity.
1. Polish the overlaps via inference of a repeat model.
1. Construct a string graph and generate contigs from it.

### Modules

#### datander ([Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/notebooks/usage/1.%20datander.ipynb))

Datander has been developed by Dr. Gene Myers as a program in his module named [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER).

#### datruf ([Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/notebooks/usage/2.%20%20datruf.ipynb))

Datruf detects tandem repeat units with the output of datander. Minimum required unit length to be accurately inferred by datruf is approximately >50 bp due to the resolution of long-read alignment on which datander relies. The size of centromeric tandem repeat units is, however, known to be generally longer than it (e.g. ~120 bp and ~360 bp in _Drosophila_).

#### ReadViewer ([Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/notebooks/usage/3.%20ReadViewer.ipynb))

This submodule offers a class for visualizing reads with tandem repeats found by datander and units by datruf.
