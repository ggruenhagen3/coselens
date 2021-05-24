# coselenss
COnditional SELection on the Excess of NonSynonymous Substitutions (coselenss) is a R package to detect conditional selection in genes between two groups of individuals. It makes extensive use of the ```dndscv``` R package. In short, the dndscv method estimates the number of nonsynonymous mutations that would be expected in the absence of selection by combining a nucleotide substitution model (196 rates encompassing all possible substitutions in each possible trinucleotide context, estimated from synonymous mutations in the dataset) and a set of genomic covariates which greatly improve the performance of the method at low mutation loads. Then, dn/ds is calculated as the ratio between the observed number of nonsynonymous mutations (n_obs) and its neutral expectation (n_exp). Finally, the significance of dn/ds is computed through a likelihood ratio test, where the null hypothesis corresponds to dn/ds=1.

A full length tutorial on how to use coselenss can be found here.

# Installation
coselenss makes heavy use of ```dndscv```, but a slightly customized version of this package is used under-the-hood and does **not** require installation. However, coselenss does require the following dependencies: ```BiocManager```, ```devtools```, ```seqinr```, ```MASS```, ```GenomicRanges```, ```Biostrings```, ```IRanges```, and ```MASS```. These can be installed by either ```install.packages()``` or ```BiocManager::install()```. To install coselenss use:
```
devtools::install_github("ggruenhagen3/coselenss")
```

# Input
A table of mutations for two groups of individuals is required. Each table of mutations should have 5 columns: sampleID, chr, pos, ref, alt. Only list independent events as mutations. An example of the format of the table:
.

|sampleID | chr | pos | ref | alt|
|---------|-----|-----|-----|----|
|sample1  | 1   | 123 | A   | T  |
|sample1  | 5   | 456 | C   | G  |
|sample2  | 2   | 321 | T   | A  |
|sample3  | 3   | 789 | A   | G  |
|sample3  | 11  | 987 | G   | C  |

# Output
text2

# How to Use
text3
