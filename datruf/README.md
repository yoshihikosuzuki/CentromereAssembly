# datruf

Datruf analyzes relatively long (>100 bp unit size) tandem repeats using datander in DAMASKER software developed by Dr. Gene Myers.

Datruf mainly consists of three parts: Runner, Viewer, and Plotter. Read the Jupyter notebook (in nbviewer because its file size is large) for basic usage.

* Runner executes datruf with many number of reads and output the results.

* Viewer receives one read ID and draw figures that help understand about the sequence structure of the read.

* Plotter requires results of datruf or other software for tandem repeat analysis (TRF, mTR, ...) and draws figures for comparison of these results.