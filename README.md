# fastqc

This is a python script, trying to replicate some and add specific functionality to the FastQC tool[[1]](#1).

## Python version

This script is developed for Python version 3.11.

## Installation

1. Clone repository: `git clone https://github.com/dgruending/fastqc.git`
1. Install requirements: `pip install -r requirements.txt`
1. Call script: `python [-O] fastqc.py -f FILE [args]` Consider using `python -O` for optimization.

## Command line arguments

* [-f/--file] (required): path to an un-/compressed fastq file; multiple files are not supported at this moment
* [-o/--output][default=working directory]: directory for output data (csv files and plots as html files); makes it with parent directories, if they not exist
* [-m/--motifs]: motifs to match, ex. adapter sequences; takes 1+ arguments; arguments need to be separated with at least 1 space
* [-a/--polya]: minimal length of 'A' to be considered a poly A tail; integer > 0
* [-p/--plot][default=False]: toogle plotting
* [-n/--njobs][default=1]: number of processes to spawn for aggregation and evaluation; integer, if less than < 1 the default is used; processes are spawned in addition to the main process
* [-b/--bar]: Sequence coding a barcode at the front of all sequences; Format: `NNXXXNN`, where `N` is any character and is ignored by the program (leading characters are necessary) and `X` marks the Barcode characters (case sensitive)

## Example
The repository contains in the data directory a compressed fastq file with 100,000 entries taken from the first reads of a bigger sequencing run[[2]](#2).

For more see the [wiki](https://github.com/dgruending/fastqc/wiki/Example).

# References
<a id="1">[1]</a>
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
<a id="2">[2]</a>
[SRR24391474](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR24391474&display=metadata)
