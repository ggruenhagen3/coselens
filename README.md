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
Coselenss was created in order to discover epistatic interactions between cancers genes in specific cancer types. In other words, we detected conditional selection in cancer genes when mutations in another cancer gene were present/absent. As an example, let's take patients with  COAD (colon cancer) and split them into two groups, those with mutations in BRAF and those without mutations in APC.

```
library("coselenss")
group1 = coselenss::group1  # mutations from patients with    mutations in APC
group2 = coselenss::group2  # mutations from patients without mutations in APC
```

Next, let's detect conditional selection in genes with APC in COAD (runtime ~5.5 minutes on a laptop).

```
coselenss_res = coselenss(group1, group2)
```

Use ```head(coselenss_res)```, the output should look like this:

| gene_name | num.drivers.group1 | num.drivers.group2 | pmis | ptrunc | pall | pind| pglobal | qglobal |
|-----------|--------------------|--------------------|------|--------|------|-----|---------|---------|
|A1BG       |0.0139176502 | -0.001275437 |0.4260251|1.0000000|0.7284637|0.4681041|0.7078692|1|
|A1CF       |0.0130513697 |-0.005499302|0.2791135|1.0000000|0.5567157|1.0000000|0.8827844|1|
|A2M        |0.0640460226 |0.003994582|0.0855121|0.8044341|0.2211083|NaN|NaN|NaN|
|A2ML1      |0.0448769648 |0.031765610|0.5995085|0.6697172|0.7954594|0.1085713|0.2978855|1|
|A3GALT2    |-0.0007904211|-0.003739818|0.8294326|1.0000000|0.9770623|1.0000000|0.9997349|1|
|A4GALT     |0.0197142366 |0.002747199|0.4765606|1.0000000|0.7761870|0.1472146|0.3621349|1|

Let's genes with significant by doing the following

```
coselenss_res[which(coselenss_res$qglobal < 0.05),]
```

The output should look like this:
| gene_name | num.drivers.group1 | num.drivers.group2 | pmis | ptrunc | pall | pind| pglobal | qglobal |
|-----------|--------------------|--------------------|------|--------|------|-----|---------|---------|
|BMPR2      |0.2976423|0.004456581|5.019896e-04|4.728058e-09|8.374397e-11|0.0007870415|2.066347e-12|2.050953e-08|
|BRAF       |0.9742985|-0.005574817|7.599226e-29|7.343197e-03|2.942258e-29|0.3863245615|0.000000e+00|0.000000e+00|

We detected two genes, but BRAF was the gene used to separate individuals in the beginning, so it's expected that it should be significant. We can remove BRAF because it is the trivial solution, but we have just discovered that there may be conditional selection between BRAF and BMPR2 in COAD. Feel free to think outside the box and make discories of your own using our tool!
