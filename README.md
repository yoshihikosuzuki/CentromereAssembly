# ECA: Experimental Centromere Assembler

A bundle of modules for centromere assembly with long reads. It currently supports only PacBio CCS reads, although we are planning to adapt to more noisy reads.

Jupyter Notebooks used for assembly with _Drosophila_ CCS reads are available [here](https://mlab.cb.k.u-tokyo.ac.jp/~yoshihiko_s/jupyter_nbs_for_submission_1210.zip) (commit: `c2f199e`). Note that our implementation currently depends on software named Consed, an unpublished work by Dr. Gene Myers (be careful there is a distinct program with the same name), to compute consensus sequences, and thus one cannot directly run these codes. We will prepare alternative codes for it soon.

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

### Step-by-step running example (headings are links to Jupyter Notebooks)

#### 1. [Detect tandem repeats from the reads](https://nbviewer.jupyter.org/github/yoshihikosuzuki/ECA_docs/blob/master/1.%20datander.ipynb)

This is performed with datander, a program in the modified [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER) module that is originally developed by Dr. Gene Myers.

#### 2. [Detect tandem repeat units from the tandem repeats](https://nbviewer.jupyter.org/github/yoshihikosuzuki/ECA_docs/blob/master/2.%20datruf.ipynb)

This module recieves the output of datander, and split each tandem repeat into a set of unit sequences. Several filterings, e.g. removing noisy units, are applied during it.

#### [Visualize how tandem repeats are detected from reads](https://nbviewer.jupyter.org/github/yoshihikosuzuki/ECA_docs/blob/master/3.%20ReadViewer.ipynb)

This module offers a visualization of reads with annotations of tandem repeats detected above.

#### 3. [Extract reads having tandem repeat units to be assembled](https://nbviewer.jupyter.org/github/yoshihikosuzuki/ECA_docs/blob/master/4.%20TRReadFilter.ipynb)

Through looking at distributions of unit length, unit count, and/or cooccurrence, determine which tandem repeats you assemble. And by specifying some parameters, this module extracts a subset of reads that contain the focal tandem repeat units.

#### 4. Compute overlaps between the filtered reads

#### 5. Polish the overlaps via repeat model inference

#### 6. Generate contigs
