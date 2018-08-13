# DAzzler Centromeric Monomer AdjuSTER

Dacmaster recieves the fasta file of the units reported by datruf, and determine a set of representative centromeric monomers from it. This task is done by

1. Detecting peaks in tandem repeat unit length (which are considered as those of centromeric monomers),
2. Collecting units around each peak
3. Clustering units which are the first units of the first tandem repeats not starting from the end point of each read
4. Adjusting start positions of all monomers by the most similar representative monomer

## How to use

See [Jupyter notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/dacmaster/docs/Usage.ipynb).

# DAzzler CENTromere clustERing module (combined to dacmaster)

Dacenter recieves a set of tandem repeat unit sequences as input from datruf, and performs clustering of them with several resolutions required in each stage of the assembly.

Dacenter performs following steps:

* Align start positions of the units
* Peak unit lengths detection
* Detect representative unit sequences

In addition, dacenter provides following visualizations:

* Unit length histogram
* Clustering details

## Usage

See [Jupyter notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/CentromereAssembly/blob/master/dacenter/docs/Usage.ipynb) (maybe a little bit slow to load due to its large file size ~40MB).