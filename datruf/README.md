# DAzzler Tandem Repeat Unit Finding module

Datruf detects tandem repeat units from PacBio raw reads after datander module of DAMASKER by Gene Myers identifies the tandem repeat regions. The size of the tandem repeat units accurately infered by datruf is relatively long (basically >100 bp) due to the resolution of long-read alignment used. However, datruf is developed for centromere assembly, whose unit size is known to be generally long (e.g. 171 bp for human).

Datruf consists of three parts:

* **Runner** applies core algorithm of datruf to arbitrary number of reads without visualization and output the results (read ID, start/end position of tandem repeat, estimated unit length, borders of units, consensus unit sequence) into a file.

* **Viewer** receives one read ID and draw figures that are not told by Runner but help to understand the sequence structure of the read.

* **Plotter** requires results of datruf and/or other software for tandem repeat analysis (TRF, mTR, etc.) and draws figures useful for comparison of these results.

## How to use

See [Jupyter notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/datruf/docs/Usage.ipynb).
