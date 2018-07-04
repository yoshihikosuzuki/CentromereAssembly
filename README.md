# CentromereAssembly

A budle of python modules for centromere assembly with PacBio long reads:

* datruf
  * detects tandem repeat regions (with datander by Gene Myers) and tandem repeat units
* dacenter
  * infers centromeric satellite repeat units (monomers) from all the tandem repeat units
  * performs clustering of the monomers at multiple, stage-specific resolusions
  * assigns monomers in the reads into the representative (consensus) of the clusters
* dacembler
  * layouts the reads encoded with the representative monomers

Every module offers some visualization functions (examples will be shown here).

## Requirements

Python must be version 3. You need to install modified DAMASKER and BITS submodules specified. BITS needs [pbcore](https://github.com/yoshihikosuzuki/pbcore) modified for python3.