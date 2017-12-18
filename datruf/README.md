# datruf

Datruf analyzes tandem repeats with relatively long (>100 bp unit size) using output of datander in DAMASKER software developed by Dr. Gene Myers.

Datruf consists of three parts:

* **Runner** applies core algorithm of datruf to arbitrary number of reads without visualization and output the results (read ID, start/end position of tandem repeat, estimated unit length, borders of units, consensus unit sequence) into a file.

* **Viewer** receives one read ID and draw figures that are not told by Runner but help to understand the sequence structure of the read.

* **Plotter** requires results of datruf and/or other software for tandem repeat analysis (TRF, mTR, ...) and draws figures useful for comparison of these results.

## How to use

Read [Jupyter notebook](http://nbviewer.jupyter.org/github/yoshihikosuzuki/RepeatAssembly/blob/master/datruf/Usage%20and%20examples%20of%20datruf%20viewer%20and%20runner.ipynb).
