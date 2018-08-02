# DAzzler Centromeric Monomer AdjuSTER

Dacmaster recieves the fasta file of the units reported by datruf, and determine a set of representative centromeric monomers from it. This task is done by

1. Detecting peaks in tandem repeat unit length (which are considered as those of centromeric monomers),
2. Collecting units around each peak
3. Clustering units which are the first units of the first tandem repeats not starting from the end point of each read
4. Adjusting start positions of all monomers by the most similar representative monomer

## How to use

See [Jupyter notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/dacmaster/docs/Usage.ipynb).
