Probhap: probabilistic single-individual haplotyping
====================================================

Probhap is a software package for reconstructing parental haplotypes
from long reads. It implements the algorithm proposed in:

```
V. Kuleshov, Probabilistic single-individual haplotyping, Bioinformatics (2014) 30 (17): 379-385.
```

Probhap works best with very long reads at a relatively shallow coverage (<= 12X). Its advantages
are improved accuracy as well as confidence scores at phased positions that let the user tradeoff
accuracy for haplotype completeness.

## Installation

To install Probhap, clone the git repo and make sure the folder is in your `PYTHONPATH`:

```
git clone git@github.com:kuleshov/ProbHap.git;
cd probhap;
export PYTHONPATH="`pwd`:$PYTHONPATH"
```

### Sample data

To download the data used for evaluating Probhap in the 2014
ECCB/Bioinformatics paper, type:

```
wget http://www.stanford.edu/~kuleshov/probhap/phasing-matrices.tar.gz
wget http://www.stanford.edu/~kuleshov/probhap/validation.tar.gz

tar -zxvf phasing-matrices.tar.gz
tar -zxvf validation.tar.gz
```

The first file contains phasing matrices derived from fosmid reads in the paper
of Duitama et al. (2010).
The second file contains the true phase of NA12878 determined by trio phasing.

### Requirements

Probhap requires numpy and scipy-weave to run. During the first run, it will compile some inline
C code using Weave. If these package are not on your system, install them using `pip` or `easy-install`.

## Running Probhap

Probhap takes as input a set of long reads in the same format as Probhap. See the 

## Reproducing experiments
