Probhap: probabilistic single-individual haplotyping
====================================================

Probhap is a software package for reconstructing parental haplotypes
from long reads. It implements the algorithm proposed in:

```
V. Kuleshov, Probabilistic single-individual haplotyping, Bioinformatics (2014) 30 (17): 379-385.
```

Probhap works best with very long reads at a relatively shallow coverage (<= 12X). Its advantages
are improved accuracy and an assessment of confidence through probability scores.

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
C code using Weave. If these packages are not on your system, install them using `pip` or `easy-install`.

## Running Probhap

Probhap takes as input a set of long reads in the same format as RefHap.
See the paper by Duitama et al. (2010) and the documentation of the algorithm
RefHap for more details.

To run Probhap, type:

```
python probhap.py \
  --reads data/input/chr22.matrix.SORTED \
  --parsed-reads chr22.probhap.reads \
  --phase chr22.probhap.out \
  --assignments chr22.probhap.assignments
```

This will take the reads in `chr22.matrix.SORTED` and produce phased blocks `chr22.probhap.out`
in the same format as that of RefHap as well as `chr22.probhap.assignments`, which assigns
reads to blocks. Before performing phasing, `ProbHap` does a merging
and filtering step on the reads; the filtered reads are saved in `chr22.probhap.reads`.

Kuleshov (2014) also describes a post-processing heurstic on the phased blocks that typically
results in a noticeable improvement in performance. To apply this post-processing step, type:

```
python probhap.py \
  --filtered-reads chr22.probhap.reads \
  --assignments chr22.probhap.assignments
  --blocks chr22.probhap.out \
  --corrected-blocks chr22.probhap.corrected.out \
```

## Reproducing experiments

To reproduce the experiments of Kuleshov (2014), use the script in the `/experiments` folder.
This requires first downloading the RefHap and MixSIH packages and setting their paths
in the `Makefile`. Also, make sure that you downloaded the sample data to `/data`.

To run the experiments, type:
```
make run
make other
```

## Feedback

Please send comments and bug reports to `kuleshov@stanford.edu`.
