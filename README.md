# coselenss
COnditional SELection on the Excess of NonSynonymous Substitutions (coselenss) is an R package to detect gene-level differential selection between two groups of samples. If the samples are grouped based on the value of a binary variable (e.g., the presence/absence of some environmental stress or phenotypic trait), coselenss identifies genes that are differentially selected depending on the grouping variable and provides maximum likelihood estimates of the effect sizes. Coselenss makes extensive use of the ```dndscv``` R package. In short, the dndscv method estimates the number of nonsynonymous mutations that would be expected in the absence of selection by combining a nucleotide substitution model (196 rates encompassing all possible substitutions in each possible trinucleotide context, estimated from synonymous mutations in the dataset) and a set of genomic covariates which greatly improve the performance of the method at low mutation loads (read more about dndscv [here](https://github.com/im3sanger/dndscv). Then, dn/ds is calculated as the ratio between the observed number of nonsynonymous substitutions (n_obs) and its neutral expectation (n_exp). Finally, the significance of dn/ds is computed through a likelihood ratio test, where the null hypothesis corresponds to dn/ds=1.

Coselenss expands the dndscv method in two substantial ways. First, it quantifies the difference (rather than the ratio) between the observed number of nonsynonymous substitutions and its neutral expectation (Δn = n_obs - n_exp). This change of perspective is most relevant for applications in which the variable of interest is the number (rather than the fraction) of mutations subject to positive (or negative) selection, such as for the inference of driver mutations in cancer. Second, it uses a modified likelihood ratio test that allows for comparison of Δn between two sets of samples, whose mutation rates and profiles are independently estimated. 

A full length tutorial on how to use coselenss can be found here.

# Installation
Coselenss makes heavy use of ```dndscv```, but a slightly customized version of this package is used under-the-hood and does **not** require installation. However, coselenss does require the following dependencies: ```BiocManager```, ```devtools```, ```seqinr```, ```MASS```, ```GenomicRanges```, ```Biostrings```, ```IRanges```, and ```MASS```. These can be installed by either ```install.packages()``` or ```BiocManager::install()```. To install coselenss use:
```
devtools::install_github("ggruenhagen3/coselenss")
```

# Input
* group1: mutation file for the first group of samples (for example, those that possess the trait of interest)
* group2: mutation file for the second group of samples (for example, those that do NOT possess the trait of interest)
* subset.genes.by (optional): genes to subset results by
* refdb (optional): reference database of coding sequences as an .rda file, the default is GRCh37/hg19

The input parameters group1 and group2 are two dataframes of mutations, one for each group of samples. Each dataframe of mutations should have 5 columns: sampleID, chr (chromosome), pos (position within the chromosome), ref (reference base), alt (mutated base). Only list independent events as mutations. An example of the format of the table:

|sampleID | chr | pos | ref | alt|
|---------|-----|-----|-----|----|
|sample1  | 1   | 123 | A   | T  |
|sample1  | 5   | 456 | C   | G  |
|sample2  | 2   | 321 | T   | A  |
|sample3  | 3   | 789 | A   | G  |
|sample3  | 11  | 987 | G   | C  |

By default, coselenss assumes that the mutation data is mapped to the GRCh37/hg19 assembly of the human reference genome. To use coselenss with different species or assemblies, an alternative reference database (RefCDS object) must be provided with the option refdb. The generation of alternative reference databases can be done using the dndscv package and is explained in [this tutorial](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html):

# Output
Coselenss returns a dataframe with effect sizes and p-values for differential selection in of the reference genome. If a list of genes is provided through the subset.genes.by option, the results and Benjamini-Hochberg corrections are restricted to those genes. The columns returned are described as follows:
* gene_name: name of gene in which differential selection was studied.
* num.drivers.group1: estimate of the excess of non-synonymous mutations (Δn) in group 1, with respect to the neutral expectation. In the absence of negative selection, this number corresponds to the average number of driver (i.e., positively selected) mutations per sample  in that gene.
* num.drivers.group2: idem, for group 2.
* pmis: p-value for differential selection in missense mutations.
* ptrunc: p-value for differential selection in truncating mutations.
* pall: p-value for differential selection in all substitutions (pmis and ptrunc).
* pind: p-value for differential selection in small indels.
* pglobal: Fisher's combined p-value for pall and pind.
* qglobal: q-value of pglobal using Benjamini-Hochberg correction.

Note that the Fisher’s combined value of pglobal/qglobal may be too conservative if the sensitivity of the pall or pind test is low. In our experience with cancer somatic mutation data, the indel test (pind) only reaches acceptable sensitivity levels if the sample size is large (>100) and indels are frequent. Otherwise, the low sensitivity of the indel test results in non-significant values of pglobal, despite the presence of substantial differential selection on substitutions. Thus, we recommend using pall/qall to assess significance of differential selection on substitutions, while restricting pglobal/qglobal to datasets with large sample sizes and genes in which indels are the main subject of selection.

# Example
Coselenss was developed to discover epistatic interactions between cancer genes in specific cancer types. To do that, we searched for differential selection in cancer genes when mutations in another cancer gene were present/absent. As an example, we took a somatic mutation dataset built from biopsies from patients with colorectal cancer (COAD) and split them into two groups, those with mutations in APC and those without mutations in APC. Let's begin, by loading these two mutation datasets.

```
library("coselenss")
data("group1", package = "coselenss")  # mutations from patients with    mutations in APC
data("group2", package = "coselenss")  # mutations from patients without mutations in APC
data("cancer_genes", package = "coselenss")  # cancer genes
```

Next, let's detect differential selection in genes with APC in COAD (runtime ~5.5 minutes on a laptop).

```
coselenss_res = coselenss(group1, group2, subset.genes.by = cancer_genes)
```

Use ```head(coselenss_res)```, the output should look like this:

| gene_name | num.drivers.group1 | num.drivers.group2 | pmis | ptrunc | pall | pind| pglobal | qglobal |
|-----------|--------------------|--------------------|------|--------|------|-----|---------|---------|
|ABL1|-8.257161e-05|0.021727512 |0.20264758| 1.0000000| 0.4441491|0.5359755| 0.5797215|       1|
|ACO1|-1.254983e-03|-0.004122997 |0.94708034| 0.6725580| 0.9125477|0.2687322| 0.5899164|       1|
|ACVR1|3.590744e-04|0.008633316 |0.40628237| 1.0000000| 0.7083432|1.0000000| 0.9525987|       1|
|ACVR1B|3.921533e-02|0.039206586 |0.95634987| 0.8988483| 0.9904685|0.3391245| 0.7023388|       1|
|ACVR2A|1.293825e-02|0.018095464|0.43847259| 0.6886162| 0.6835661|1.0000000| 0.9436165|      1|
|ACVR2B|8.398367e-03|0.028399806|0.08498634| 0.6456863| 0.2041043|1.0000000| 0.5284514|       1|

Let's find the genes subject to significant significant differential selection by doing the following:

```
coselenss_res[which(coselenss_res$qglobal < 0.05),]
```

The output should look like this:
| gene_name | num.drivers.group1 | num.drivers.group2 | pmis | ptrunc | pall | pind| pglobal | qglobal |
|-----------|--------------------|--------------------|------|--------|------|-----|---------|---------|
|APC|0.98535779|-0.04327314| 4.599838e-02| 1.669909e-08|1.661006e-08| 0.2386762| 8.065981e-08| 3.532899e-05|
|BRAF|0.03621197|0.22591651| 2.276400e-08| 1.821872e-01|6.747822e-08| 0.4244034| 5.260378e-07| 1.152023e-04|

We detected two genes, but APC was the gene used to separate individuals in the beginning, so it's expected that it should be significant. We can remove APC because it is the trivial solution, but we have just discovered that there may be conditional selection between APC and BRAF in COAD. Feel free to think outside the box and make discoveries of your own using our tool!
