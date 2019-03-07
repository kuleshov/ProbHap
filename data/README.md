Probhap evaluation data
=======================

To unzip the data used for evaluating Probhap in the 2014
ECCB/Bioinformatics paper, type:

```
unzip phasing-matrices.zip
unzip validation.zip
```

The first file contains phasing matrices derived from fosmid reads in the paper
of Duitama et al. (2010).
The second file contains the true phase of NA12878 determined by trio phasing.
We moved these true phasing calls into the same coordinate system as the phasing
matrices by lifting the matrix coordinates to b37 and by selecting those positions
that were heterozygous in the vcf file of NA12878.
