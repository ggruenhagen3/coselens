# coselenss
COnditional SELection on the Excess of NonSynonymous Substitutions (coselenss) is a R package to detect conditional selection in genes between two groups of individuals. It makes extensive use of the ```dndscv``` R package. In short, the dndscv method estimates the number of nonsynonymous mutations that would be expected in the absence of selection by combining a nucleotide substitution model (196 rates encompassing all possible substitutions in each possible trinucleotide context, estimated from synonymous mutations in the dataset) and a set of genomic covariates which greatly improve the performance of the method at low mutation loads. Then, dn/ds is calculated as the ratio between the observed number of nonsynonymous mutations (n_obs) and its neutral expectation (n_exp). Finally, the significance of dn/ds is computed through a likelihood ratio test, where the null hypothesis corresponds to dn/ds=1.

A full length tutorial on how to use coselenss can be found here.

# Installation
Coselenss makes heavy use of ```dndscv```, but a slightly customized version of this package is used under-the-hood and does **not** require installation. However, coselenss does require the following dependencies: ```BiocManager```, ```devtools```, ```seqinr```, ```MASS```, ```GenomicRanges```, ```Biostrings```, ```IRanges```, and ```MASS```. These can be installed by either ```install.packages()``` or ```BiocManager::install()```. To install coselenss use:
```
devtools::install_github("ggruenhagen3/coselenss")
```

# Input
* group1: group of individuals (for example those that contain a mutation in a split_gene)
* group2: another group of individuals that do NOT contain a mutation in a split_gene
* subset.genes.by: genes to subset results by
* refdb: Reference database (path to .rda file)

The input parameters group1 and groups are a dataframe of mutations for the two groups of individuals is required. Each dataframe of mutations should have 5 columns: sampleID, chr, pos, ref, alt. Only list independent events as mutations. An example of the format of the table:

|sampleID | chr | pos | ref | alt|
|---------|-----|-----|-----|----|
|sample1  | 1   | 123 | A   | T  |
|sample1  | 5   | 456 | C   | G  |
|sample2  | 2   | 321 | T   | A  |
|sample3  | 3   | 789 | A   | G  |
|sample3  | 11  | 987 | G   | C  |

# Output
A dataframe is returned with p-values for conditional selection in each gene in the human reference genome (hg19). The columns returned are described as follows:
* gene_name: name of gene that conditional selection was calculated in
* num.drivers.group1: estimate of the number of drivers in group 1 based excess of non-synonymous mutations
* num.drivers.group2: estimate of the number of drivers in group 2 based excess of non-synonymous mutations
* pmis: p-value for conditional selection in missense mutations
* ptrunc: p-value for conditional selection in truncating mutations
* pall: p-value for conditional selection in all substitutions (pmis and ptrunc)
* pind: p-value for conditional selection in small indels
* pglobal: Fisher's combined p-value for pall and pind
* qglobal: q-value of pglobal using Benjamini-Hochberg correction

Note that pglobal/qglobal may be too conservative if the sensitivity for the pall or pind test is low due to low sample sizes.

# Example
Coselenss was created in order to discover epistatic interactions between cancers genes in specific cancer types. In other words, we detected conditional selection in cancer genes when mutations in another cancer gene were present/absent. As an example, let's take patients with  COAD (colon cancer) and split them into two groups, those with mutations in APC and those without mutations in APC.

```
library("coselenss")
data("group1", package = "coselenss")  # mutations from patients with    mutations in APC
data("group2", package = "coselenss")  # mutations from patients without mutations in APC
```

Next, let's detect conditional selection in genes with APC in COAD (runtime ~5.5 minutes on a laptop).

```
coselenss_res = coselenss(group1, group2)
```

Use ```head(coselenss_res)```, the output should look like this:

| gene_name | num.drivers.group1 | num.drivers.group2 | pmis | ptrunc | pall | pind| pglobal | qglobal |
|-----------|--------------------|--------------------|------|--------|------|-----|---------|---------|
|A1BG       |-0.0002757828|0.003490537|0.6774850|1.0000000|0.9171490|0.28669884|0.6141904|1|
|A1CF       |-0.0033004353|-0.008295631|0.6058041|1.0000000|0.8753205|1.00000000|0.9918827|1|
|A2M        |0.0037921837|0.029329399|0.3043810|0.5185936|0.4791242|NaN|NaN|NaN|
|A2ML1      |0.0295412248|0.040366524|0.7277502|0.7012937|0.8744515|0.09369533|0.2869149|1|
|A3GALT2    |-0.0047941596|-0.006902210|0.8620675|1.0000000|0.9850200|1.00000000|0.9998872|1|
|A4GALT     |-0.0012767999|0.017041030|0.1117167|1.0000000|0.2822721|0.23622026|0.2472351|1|

Let's genes with significant by doing the following

```
coselenss_res[which(coselenss_res$qglobal < 0.05),]
```

The output should look like this:
| gene_name | num.drivers.group1 | num.drivers.group2 | pmis | ptrunc | pall | pind| pglobal | qglobal |
|-----------|--------------------|--------------------|------|--------|------|-----|---------|---------|
|APC        |0.98535779|-0.04327314|4.599838e-02|1.669909e-08|1.661006e-08|0.2386762|8.065981e-08|0.001610938|
|BRAF       |0.03621197|0.22591651|2.276400e-08|1.821872e-01|6.747822e-08|0.4244034|5.260378e-07|0.005253014|

We detected two genes, but APC was the gene used to separate individuals in the beginning, so it's expected that it should be significant. We can remove APC because it is the trivial solution, but we have just discovered that there may be conditional selection between APC and BRAF in COAD. Feel free to think outside the box and make discories of your own using our tool!
