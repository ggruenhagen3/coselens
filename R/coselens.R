#' COnditional SELection on the Excess of NonSynonymous Substitutions (coselens)
#'
#' coselens identifies genes that are differentially selected depending on the grouping variable and provides maximum likelihood estimates of the effect sizes.
#' Given an input of mutations from two groups of patients, dN/dS is calculated in each of the groups,
#' and the excess of non-synonymous substitutions between them is found as an estimate for the difference
#' in the number of driver mutations. For example, one group can be mutations from individuals in a certain
#' cancer type with mutations in a specific gene vs those without mutations in the specific gene.
#' A p value is returned for differential selection in each gene in the reference genome (default: hg19).
#'
#' @param group1 group of individuals (for example those that contain a mutation in a split_gene)
#' @param group2 another group of individuals that do NOT contain a mutation in a split_gene
#' @param subset.genes.by genes to subset results by
#' @param sequenced.genes the gene_list paramater from dndscv, which is a list of genes to restrict the analysis (use for targeted sequencing studies)
#' @param ... other parameters passed to dncdscv, all defaults from dndscv are used except max_muts_per_gene_per_sample is set to Infinity
#'
#' @return coselens returns a list containing four dataframes: 1) "summary" a summary of coselens output that should be sufficient for most users, 2) "full" which includes more detailed output, 3) fitted substitution models for group 1 from dndscv, and 4) fitted substitution models for group 2 from dndscv
#' @return - gene_name: name of gene that conditional selection was calculated in
#' @return - num.driver.sub.group1: estimate of the number of drivers in group 1 based excess of non-synonymous mutations
#' @return - num.driver.sub.group2: estimate of the number of drivers in group 2 based excess of non-synonymous mutations
#' @return - num.driver.ind.group1: estimate of the number of drivers in group 1 based excess of indels
#' @return - num.driver.ind.group2: estimate of the number of drivers in group 2 based excess of indels
#' @return - psub: p-value for conditional selection in non-synonymous substitutions
#' @return - pind: p-value for conditional selection in indels
#' @return - pglobal: Fisher's combined p-value for psub and pind
#' @return - qsub: q-value of psub using Benjamini-Hochberg correction
#' @return - qind: q-value of pind using Benjamini-Hochberg correction
#' @return - qglobal: q-value of pglobal using Benjamini-Hochberg correction
#'
#' @return - num.driver.mis.group1: estimate of the number of drivers in group 1 based excess of missense mutations
#' @return - num.driver.mis.group2: estimate of the number of drivers in group 2 based excess of missense mutations
#' @return - num.driver.trunc.group1: estimate of the number of drivers in group 1 based excess of truncating mutations
#' @return - num.driver.trunc.group2: estimate of the number of drivers in group 2 based excess of truncating mutations
#' @return - pmis: p-value for conditional selection in missense mutations
#' @return - ptrunc: p-value for conditional selection in truncating mutations
#' @return - psub.group1: p-value for selection in group 1 for non-synonymous substitutions
#' @return - psub.group2: p-value for selection in group 2 for non-synonymous substitutions
#' @return - pmis.group1: p-value for selection in group 1 for missense substitutions
#' @return - pmis.group2: p-value for selection in group 2 for missense substitutions
#' @return - ptrunc.group1: p-value for selection in group 1 for truncating substitutions
#' @return - ptrunc.group2: p-value for selection in group 2 for truncating substitutions
#' @return - pind.group1: p-value for selection in group 1 for indels
#' @return - pind.group2: p-value for selection in group 2 for indels
#' @return - qsub.group1: q-value of psub.group1 using Benjamini-Hochberg correction
#' @return - qsub.group2: q-value of psub.group2 using Benjamini-Hochberg correction
#' @return - qmis.group1: q-value of pmis.group1 using Benjamini-Hochberg correction
#' @return - qmis.group2: q-value of pmis.group2 using Benjamini-Hochberg correction
#' @return - qtrunc.group1: q-value of ptrunc.group1 using Benjamini-Hochberg correction
#' @return - qtrunc.group2: q-value of ptrunc.group2 using Benjamini-Hochberg correction
#' @return - qind.group1: q-value of pind.group1 using Benjamini-Hochberg correction
#' @return - qind.group2: q-value of pind.group2 using Benjamini-Hochberg correction
#'
#' @export

coselens = function(group1, group2, subset.genes.by = NULL, sequenced.genes = NULL, refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = 3000, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, outp = 3, numcode = 1) {
  # Find dN/dS and CIs
  tryCatch({
    # Calculate dN/dS and the confidence intervals
    group1_dndsout <- coselens::dndscv(group1, gene_list = sequenced.genes, refdb = refdb, sm = sm, kc = kc, cv = "hg19", max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, outp = outp, numcode = numcode, outmats=T, outmutrates = T, wg = F, ex = T)
    group2_dndsout <- coselens::dndscv(group2, gene_list = sequenced.genes, refdb = refdb, sm = sm, kc = kc, cv = "hg19", max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, outp = outp, numcode = numcode, outmats=T, outmutrates = T, wg = F, ex = T)

    # Calculate nonsynonymous mutation excess
    group1_ex <- calc_ex(group1_dndsout)
    group2_ex <- calc_ex(group2_dndsout)
    group1_ex$ex_tot <- group1_ex[,2] + group1_ex[,3]
    group2_ex$ex_tot <- group2_ex[,2] + group2_ex[,3]
    ex_values <- merge(x=group1_ex[,c("gene_name","ex_mis","ex_non","ex_tot")], y=group2_ex[,c("gene_name", "ex_mis", "ex_non", "ex_tot")], by ="gene_name", suffixes=c(".group1", ".group2"))

    # Calculate indel mutation excess
    group1_ex_ind <- calc_ex_ind(group1_dndsout)
    group2_ex_ind <- calc_ex_ind(group2_dndsout)
    group1_ex_ind$ex_ind = as.numeric(as.vector(group1_ex_ind$ex_ind))
    group2_ex_ind$ex_ind = as.numeric(as.vector(group2_ex_ind$ex_ind))
    ex_ind_values <- merge(x=group1_ex_ind[,c("gene_name","ex_ind")], y=group2_ex_ind[,c("gene_name","ex_ind")], by="gene_name", suffixes=c(".group1",".group2"))
    ex_values <- merge(x=ex_values, y=ex_ind_values, by="gene_name")

    # Put the data in the same order
    group1_sel_cv = group1_dndsout$sel_cv[order(as.numeric(rownames(group1_dndsout$sel_cv))),,drop=FALSE]
    group2_sel_cv = group2_dndsout$sel_cv[order(as.numeric(rownames(group2_dndsout$sel_cv))),,drop=FALSE]

    # Null Hypothesis Log-Likelihoods
    h0_sel_cv <<- lik_func(group1_dndsout, group2_dndsout)

    # Merge the data into one structure
    lldf <- merge(x=h0_sel_cv[,c("gene_name", "llall", "llmis", "lltrunc", "llmis_check_A", "lltrunc_check_B", "llind0", "llind1")], y=group1_sel_cv[,c("gene_name", "wmis_cv", "wnon_cv", "llall")], by="gene_name", suffixes=c(".0", ".group1"))
    lldf <- merge(x=lldf, y=group2_sel_cv[,c("gene_name", "wmis_cv", "wnon_cv", "llall")], by ="gene_name", suffixes=c(".group1", ".group2"))
    colnames(lldf)[c(ncol(lldf)-2, ncol(lldf)-1, ncol(lldf))] <- c("wmis_cv.group2", "wnon_cv.group2", "llall.group2")

    # Alternative Hypothesis Log-Likelihoods for SNV
    lldf$llall.1   <- lldf$llall.group1 + lldf$llall.group2
    lldf <- lldf[,c(1, ncol(lldf), 3:ncol(lldf)-1)] # reorder llall.1 to the front

    # Optionally subset results by list of genes before p value calculation
    if (length(subset.genes.by) > 0) {
      lldf      <- lldf[     which(lldf$gene_name      %in% subset.genes.by),]
      ex_values <- ex_values[which(ex_values$gene_name %in% subset.genes.by),]
    }

    # Calculate p-values
    lldf$pall    <- pchisq(- 2 * (lldf$llall.0 - lldf$llall.1), df=2, lower.tail = FALSE)
    lldf$pmis    <- pchisq(- 2 * (lldf$llmis   - lldf$llall.1), df=1, lower.tail = FALSE)
    lldf$ptrunc  <- pchisq(- 2 * (lldf$lltrunc - lldf$llall.1), df=1, lower.tail = FALSE)
    lldf$pind    <- 1 - pchisq(2*(lldf$llind1 - lldf$llind0),df=1)
    lldf$pglobal <- 1 - pchisq(-2 * (log(lldf$pall) + log(lldf$pind)), df = 4)
    lldf <- lldf[order(lldf$pglobal),] # order data by p-values

    # Adjust p-values
    lldf$qind = p.adjust(lldf$pind, method = "BH")
    lldf$qall = p.adjust(lldf$pall, method = "BH")
    lldf$qglobal = p.adjust(lldf$pglobal, method = "BH")

    # Add mutation excess data
    lldf <- merge(x = lldf, y = ex_values[,c("gene_name","ex_tot.group1","ex_tot.group2", "ex_ind.group1", "ex_ind.group2", "ex_mis.group1", "ex_mis.group2", "ex_non.group1", "ex_non.group2")], by = "gene_name")

    # Single Group Tests
    single.test.names = c("psub.group", "pind.group", "pmis.group", "ptrunc.group", "qsub.group", "qind.group", "qmis.group", "qtrunc.group")
    lldf = merge(x=lldf, y=group1_sel_cv[,c("gene_name", "pallsubs_cv", "pind_cv", "pmis_cv", "ptrunc_cv", "qallsubs_cv", "qind_cv", "qmis_cv", "qtrunc_cv")], by="gene_name")
    colnames(lldf)[(ncol(lldf)-length(single.test.names)+1):ncol(lldf)] = paste0(single.test.names, "1")
    lldf = merge(x=lldf, y=group2_sel_cv[,c("gene_name", "pallsubs_cv", "pind_cv", "pmis_cv", "ptrunc_cv", "qallsubs_cv", "qind_cv", "qmis_cv", "qtrunc_cv")], by="gene_name")
    colnames(lldf)[(ncol(lldf)-length(single.test.names)+1):ncol(lldf)] = paste0(single.test.names, "2")

    # Organize output
    out.list = list()
    cols.for.summary.from.lldf = c("gene_name", "ex_tot.group1", "ex_tot.group2", "ex_ind.group1", "ex_ind.group2", "pall", "pind", "pglobal", "qall", "qind", "qglobal")
    cols.for.summary.from.lldf.names = c("gene_name", "num.driver.sub.group1", "num.driver.sub.group2", "num.driver.ind.group1", "num.driver.ind.group2", "psub", "pind", "pglobal", "qsub", "qind", "qglobal")
    out.list[["summary"]] = lldf[,cols.for.summary.from.lldf]
    colnames(out.list[["summary"]]) = cols.for.summary.from.lldf.names

    cols.for.full.from.lldf = c(cols.for.summary.from.lldf, c("ex_mis.group1", "ex_mis.group2", "ex_non.group1", "ex_non.group2", "pmis", "ptrunc", paste0(single.test.names, "1"), paste0(single.test.names, "2")))
    cols.for.full.from.lldf.names = c(cols.for.summary.from.lldf.names, c("num.driver.mis.group1", "num.driver.mis.group2", "num.driver.trunc.group1", "num.driver.trunc.group2", "pmis", "ptrunc", paste0(single.test.names, "1"), paste0(single.test.names, "2")))
    out.list[["full"]] = lldf[,cols.for.full.from.lldf]
    colnames(out.list[["full"]]) = cols.for.full.from.lldf.names

    out.list[["mle_submodel_group1"]] = group1_dndsout$mle_submodel
    out.list[["mle_submodel_group2"]] = group2_dndsout$mle_submodel

    # Return the output
    print("Done.")
    return(out.list)

  }, error=function(err) {

    print("An error has occured.")
    print(err)
    return(NULL)

  }) # end tryCatch
}
