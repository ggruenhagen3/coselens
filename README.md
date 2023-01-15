# ![alt text](https://github.com/ggruenhagen3/coselens/blob/master/icon.png?raw=true)  Coselens 
COnditional SELection on the Excess of NonSynonymous Substitutions (coselens) is an R package to detect gene-level differential selection between two groups of samples. If the samples are grouped based on the value of a binary variable (e.g., the presence/absence of some environmental stress or phenotypic trait), coselens identifies genes that are differentially selected depending on the grouping variable and provides maximum likelihood estimates of the effect sizes. Coselens makes extensive use of the ```dndscv``` R package. In short, the dndscv method estimates the number of nonsynonymous mutations that would be expected in the absence of selection by combining a nucleotide substitution model (196 rates encompassing all possible substitutions in each possible trinucleotide context, estimated from synonymous mutations in the dataset) and a set of genomic covariates which greatly improve the performance of the method at low mutation loads (read more about dndscv [here](https://github.com/im3sanger/dndscv)). Then, dn/ds is calculated as the ratio between the observed number of nonsynonymous substitutions (n_obs) and its neutral expectation (n_exp). Finally, the significance of dn/ds is computed through a likelihood ratio test, where the null hypothesis corresponds to dn/ds=1.

Coselens expands the dndscv method in two substantial ways. First, it quantifies the difference (rather than the ratio) between the observed number of nonsynonymous substitutions and its neutral expectation (Δn = n_obs - n_exp). This change of perspective is most relevant for applications in which the variable of interest is the number (rather than the fraction) of mutations subject to positive (or negative) selection, such as for the inference of driver mutations in cancer. Second, it uses a modified likelihood ratio test that allows for comparison of Δn between two sets of samples, whose mutation rates and profiles are independently estimated. 

A full length tutorial on how to use coselens can be found here.

# Installation
Coselens makes heavy use of ```dndscv```, but a slightly customized version of this package is used under-the-hood and does **not** require installation. However, coselens does require the following dependencies: ```BiocManager```, ```devtools```, ```geometry```, ```seqinr```, ```MASS```, ```GenomicRanges```, ```Biostrings```, and ```IRanges```. These can be installed by either ```install.packages()``` or ```BiocManager::install()```. To install coselens use:
```
devtools::install_github("ggruenhagen3/coselens")
```

# Input
* group1: mutation file for the first group of samples (for example, those that possess the trait of interest)
* group2: mutation file for the second group of samples (for example, those that do NOT possess the trait of interest)
* subset.genes.by (optional): genes to subset results by
* sequenced.genes (optional): the gene_list paramater from dndscv, which is a list of genes to restrict the analysis (use for targeted sequencing studies)
* ... other paramters passed to dncdscv, all defaults from ```dndscv``` are used except max_muts_per_gene_per_sample is set to Infinity

The input parameters group1 and group2 are two dataframes of mutations, one for each group of samples. Each dataframe of mutations should have 5 columns: sampleID, chr (chromosome), pos (position within the chromosome), ref (reference base), alt (mutated base). Only list independent events as mutations. An example of the format of the table:

|sampleID | chr | pos | ref | alt|
|---------|-----|-----|-----|----|
|sample1  | 1   | 123 | A   | T  |
|sample1  | 5   | 456 | C   | G  |
|sample2  | 2   | 321 | T   | A  |
|sample3  | 3   | 789 | A   | G  |
|sample3  | 11  | 987 | G   | C  |

By default, coselens assumes that the mutation data is mapped to the GRCh37/hg19 assembly of the human reference genome. To use coselens with different species or assemblies, an alternative reference database (RefCDS object) must be provided with the option refdb. The generation of alternative reference databases can be done using the dndscv package and is explained in [this tutorial](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html):

# Output
Coselens returns several tables with effect sizes and p-values for differential selection in each gene of the reference genome. Coselens returns a list of 6 dataframes named substitutions, indels, missense_sub, truncating_sub, overall_mut and dndscv. The dataframes contain information on differential selection from different types of mutations. Coselens also returns the results of a joint analysis of single-nucleotide substitutions and indels in the table ```coselens_out$overall_mut```. The last element of the list is another list containing two dataframes which are the output of dndscv for group1 and group2.

Column Description of dataframes 1-4 in the output of coselens: 
* Output Column Descriptions
  * gene_name: name of gene that conditional selection was calculated in
  * num.driver.group1: estimate of the number of drivers per sample per gene in group 1
  * num.driver.group2: estimate of the number of drivers per sample per gene in group 2
  * Delta.Nd: absolute difference in the average number of driver mutations per sample (group 1 minus group 2) <details><summary>Details</summary>A major feature of Coselens is that it provides the user with biologically meaningful effect sizes for the magnitude of conditional selection. The most straightforward way to quantify effect sizes is by calculating the difference in the average number of driver mutations in the presence and absence of the condition of interest. We call that measure of effect size ΔNd (column 4, Delta.Nd). The value of ΔNd indicates, in absolute terms, how the grouping variable modifies the average number of driver mutations in a gene.</details>
  * classification: qualitative classification of the association between the grouping variable and the magnitude and sign of selection <details><summary>Details</summary>)Provides a qualitative classification of the association between the grouping variable (the condition of interest) and the magnitude and sign of selection. Independence implies that the grouping variable does not affect selection for drivers; strict dependence implies that drivers are only positively selected in the first group of samples, in which the condition of interest is met; strict inhibition implies that positive selection only acts in the second group, in which the condition is not met. Strict dependence and strict inhibition are the two extremes of conditionality, that is, cases of full conditionality. Instances of partial conditionality are labeled as facilitation and inhibition, respectively. Together with independence, these four classes cover the whole spectrum of dependencies for positively selected driver mutations. If negative (purifying) selection is dominant, other classes of dependency become possible: strict dependence with sign change, if the sign of selection changes from negative to positive when the condition is met; strict inhibition with sign change, if the sign of selection changes from positive to negative when the condition is met; aggravation, if purifying selection against mutations becomes stronger when the condition is met; and relaxation, if selection against mutations becomes weaker when the condition is met. </details>
  * dependency: the association between the grouping variable and the average number of drivers observed in a gene<details><summary>Details</summary>A quantitative measure of these associations (read details of classification), which takes values between 0/NA (no conditionality) and 1 (full conditionality).</details>
  * pval: p-value for conditional selection
  * qval: q-value for conditional selection using Benjamini-Hochberg correction of false discovery rate
# Example
Coselens was developed to discover epistatic interactions between cancer genes in specific cancer types. To do that, we searched for differential selection in cancer genes when mutations in another cancer gene were present/absent. As an example, we took a somatic mutation dataset built from biopsies from patients with colorectal cancer (COAD) and split them into two groups, those with mutations in APC and those without mutations in APC. Let's begin, by loading these two mutation datasets.

```
library("coselens")
data("group1", package = "coselens")  # mutations from patients with    mutations in APC
data("group2", package = "coselens")  # mutations from patients without mutations in APC
data("cancer_genes", package = "coselens")  # cancer genes
```

Next, let's detect differential selection in genes with APC in COAD (runtime ~5.5 minutes on a laptop).

```
coselens_res = coselens(group1, group2, subset.genes.by = cancer_genes)
```

Use ```head(coselens_res[["overall_mut"]])```, the output should look like this:

| gene_name | num.driver.group1 | num.driver.group2 | Delta.Nd | classification | dependency | pval | qval |
|-----------|-------------------|-------------------|----------|----------------|------------|------|------|
|ABL1|-0.0019200888|0.022497687|-0.024417775|independence|NA|0.5135217|1|
|ACO1|0.0049344727|-0.009434191|0.014368663|independence|NA|0.2589758|1|
|ACVR1|-0.0004685018|0.005431772|-0.005900274|independence|NA|0.9525934|1|
|ACVR1B|0.0383831658|0.043194501|-0.004811335|independence|NA|0.6044401|1|
|ACVR2A|0.0248638286|0.028849193|-0.003985365|independence|NA|0.9260953|1|
|ACVR2B|0.0077000832|0.025621900|-0.017921817|independence|NA|0.3617251|1|

Let's find the genes subject to significant significant differential selection by doing the following:

```
head(coselens_res$summary[which(coselens_res$summary$qglobal < 0.05),])
```

The output should look like this:
| gene_name | num.driver.group1 | num.driver.group2 | Delta.Nd | classification | dependency | pval | qval |
|-----------|-------------------|-------------------|----------|----------------|------------|------|------|
|APC|1.40733813|-0.05877917|1.4661173|strict|dependence|1.0000000|0.000000e+00|0.000000e+00|
|BRAF|0.03510272|0.22862023|-0.1935175|strict|inhibition|0.8060196|4.175434e-07|9.206832e-05|

We detected 2 significant genes, but APC was the gene used to separate individuals in the beginning, so it's expected that it should be significant. We can remove APC because it is the trivial solution. The q-value of BRAF is really low this provides strong evidence of conditional selection between APC and BRAF in COAD. Now that you see how the tool works, feel free to think outside the box and make discoveries of your own!
