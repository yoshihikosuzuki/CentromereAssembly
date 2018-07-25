# DAzzler Tandem Repeat Unit Finding module

Datruf detects tandem repeats and unit sequences within them only from PacBio raw reads and without any other external information. Before using datruf, you must run the `HPC.TANmask` program of [DAMASKER](https://github.com/yoshihikosuzuki/DAMASKER) originally developed by Gene Myers and slightly modified for datruf, and must have `.db` file and `.las` file.  The minimum length of tandem repeat units accurately infered by datruf is relatively long (approximately >50 bp) due to the resolution of long-read alignment on which datander relies. However, datruf is developed mainly for (core) centromere assembly, whose unit size is known to be long in general (e.g. 171 bp alpha-satellite sequence of human).

Datruf consists of three parts:

* **Runner** applies core algorithm of datruf to arbitrary number of reads without visualization and output the results (read ID, start/end position of tandem repeat, estimated unit length, borders of units, consensus unit sequence) into a file.

* **Viewer** receives one read ID and draw figures that are not told by Runner but help to understand the sequence structure of the read.

* **Plotter** requires results of datruf and/or other software for tandem repeat analysis (TRF, mTR, etc.) and draws figures useful for comparison of these results.

## How to use

See [Jupyter notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/datruf/docs/Usage.ipynb).

<iframe width="900" height="800" frameborder="0" scrolling="no" src="//plot.ly/~yoshihikosuzuki/61.embed"></iframe>
