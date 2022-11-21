#' COnditional SELection on the Excess of NonSynonymous Substitutions (coselens)
#'
#' Coselens identifies genes that are differentially selected depending on the grouping variable and provides maximum likelihood estimates of the effect sizes.
#' Given an input of mutations from two groups of patients, dN/dS is calculated in each of the groups,
#' and the excess of non-synonymous substitutions between them is found as an estimate for the difference
#' in the number of driver mutations. For example, one group can be mutations from individuals in a certain
#' cancer type with mutations in a specific gene vs those without mutations in the specific gene.
#' A p value is returned for differential selection in each gene in the reference genome (default: hg19).
#'
#' @author George Gruehnagen (Georgia Institute of Technology) and Jaime Iranzo (Centro de Biotecnologia y Genomica de Plantas - Universidad Politcnica de Madrid)
#' @details Iranzo J, Gruenhagen G, Calle-Espinosa J, Koonin EV (2022) Pervasive conditional selection of driver mutations and modular epistasis networks in cancer. Cell Reports. 40(8):111272.
#'
#' @param group1 group of samples (for example those that are subject to a given condition)
#' @param group2 another group of samples that are NOT subject to that condition (the control group)
#' @param subset.genes.by genes to subset results by
#' @param sequenced.genes the gene_list paramater from dndscv, which is a list of genes to restrict the analysis (use for targeted sequencing studies)
#' @param ... other parameters passed to dncdscv, all defaults from dndscv are used except max_muts_per_gene_per_sample is set to Infinity
#'
#' @return coselens returns a list containing six objects: 1) "substitutions", a summary table of gene-level conditional selection for single-nucleotide substitutions (including missense, nonsense, and essential splice site mutations); 2) "indels", same for small indels; (3) "missense_sub", same for missense substitutions only; (4) "truncating_sub", same for truncating substitutions only (nonsense and essential splice site mutations); (5) "overall_mut", a summary table with the combined analysis of single-nucleotide substitutions and indels; and (6) "dndscv", a list of objects with the complete output of (non-conditional) selection analyses separately run on each group of samples, as provided by the dndscv package. The first table should be sufficient for most users. If interested in indels, please note that the indel analysis uses a different null model that makes the test for conditional selection notably less sensitive than in the case of substitutions. Such lower sensitivity also extends to the "overall_mut" table. The dataframes (1-5) contain the following:
#' @return - gene_name: name of the gene that was tested for conditional selection
#' @return - num.driver.group1: estimate of the number of drivers per sample per gene in group 1
#' @return - num.driver.group2: estimate of the number of drivers per sample per gene in group 2
#' @return - Delta.Nd: absolute difference in the average number of driver mutations per sample (group 1 minus group 2)
#' @return - classification: classification of conditional selection. The most frequent classes are strict dependence (drivers only in group 1), facilitation (drivers more frequent in group 1), independence, inhibition (drivers less frequent in group 1), and strict inhibition (drivers absent from group 1). If negative selection is present, other possibilities are strict dependence with sign change (drivers positively selected in group 1 but negatively selected in group 2), strict inhibition with sign change (drivers positively selected in group 2 but negatively selected in group 1), aggravation (purifying selection against mutations becomes stronger in group 1), and relaxation (purifying selection against mutations becomes weaker in group 1).
#' @return - dependency: dependency index, measuring the association between the grouping variable (group 1 or 2) and the average number of drivers observed in a gene. It serves as a quantitative measure of the qualitative effect described in "classification". In the most common cases, a value of 1 indicates strict dependence or inhibition (drivers only observed in one group) and a value of 0 (or NA) indicates independence.
#' @return - pval: p-value for conditional selection
#' @return - qval: q-value for conditional selection using Benjamini-Hochberg correction of false discovery rate.
#'
#' @return The "dndscv" list contains two objects. Please, read the documentation of the dndscvpackage for further information about th
#' @return - dndscv_group1: output of dndscv for group 1
#' @return - dndscv_group2: output of dndscv for group 2
#'
#' @export

coselens = function(group1, group2, subset.genes.by = NULL, sequenced.genes = NULL, refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = 3000, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, outp = 3, numcode = 1, mingenecovs = 500) {
  # Find dN/dS and CIs
  tryCatch({
  
    # Calculate dN/dS and the confidence intervals        
    message("[0] Running dndscv for group 1")
    group1_dndsout <- coselens::dndscv(group1, gene_list = sequenced.genes, refdb = refdb, sm = sm, kc = kc, cv = cv, max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, outp = outp, numcode = numcode, outmats=T, outmutrates = T, mingenecovs = mingenecovs)
    message("[0] Running dndscv for group 2")
    group2_dndsout <- coselens::dndscv(group2, gene_list = sequenced.genes, refdb = refdb, sm = sm, kc = kc, cv = cv, max_muts_per_gene_per_sample = max_muts_per_gene_per_sample, max_coding_muts_per_sample = max_coding_muts_per_sample, use_indel_sites = use_indel_sites, min_indels = min_indels, maxcovs = maxcovs, constrain_wnon_wspl = constrain_wnon_wspl, outp = outp, numcode = numcode, outmats=T, outmutrates = T, mingenecovs = mingenecovs)
    
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
    lldf$qmis = p.adjust(lldf$pmis, method = "BH")
    lldf$qtrunc = p.adjust(lldf$ptrunc, method = "BH")
    lldf$qglobal = p.adjust(lldf$pglobal, method = "BH")

    # Add mutation excess data
    lldf <- merge(x = lldf, y = ex_values[,c("gene_name","ex_tot.group1","ex_tot.group2", "ex_ind.group1", "ex_ind.group2", "ex_mis.group1", "ex_mis.group2", "ex_non.group1", "ex_non.group2")], by = "gene_name")

    # Single Group Tests
    single.test.names = c("psub.group", "pind.group", "pmis.group", "ptrunc.group", "pglobal.group", "qsub.group", "qind.group", "qmis.group", "qtrunc.group", "qglobal.group")
    lldf = merge(x=lldf, y=group1_sel_cv[,c("gene_name", "pallsubs_cv", "pind_cv", "pmis_cv", "ptrunc_cv", "pglobal_cv", "qallsubs_cv", "qind_cv", "qmis_cv", "qtrunc_cv", "qglobal_cv")], by="gene_name")
    colnames(lldf)[(ncol(lldf)-length(single.test.names)+1):ncol(lldf)] = paste0(single.test.names, "1")
    lldf = merge(x=lldf, y=group2_sel_cv[,c("gene_name", "pallsubs_cv", "pind_cv", "pmis_cv", "ptrunc_cv", "pglobal_cv", "qallsubs_cv", "qind_cv", "qmis_cv", "qtrunc_cv", "qglobal_cv")], by="gene_name")
    colnames(lldf)[(ncol(lldf)-length(single.test.names)+1):ncol(lldf)] = paste0(single.test.names, "2")

    # Organize output
    full.data.from.lldf = c("gene_name", "ex_tot.group1", "ex_tot.group2", "ex_ind.group1", "ex_ind.group2", "pall", "pind","pglobal","qall", "qind", "qglobal","ex_mis.group1", "ex_mis.group2", "ex_non.group1", "ex_non.group2", "pmis", "ptrunc", "qmis", "qtrunc", paste0(single.test.names, "1"), paste0(single.test.names, "2"))
    full.data.from.lldf.names = c("gene_name", "num.driver.sub.group1", "num.driver.sub.group2", "num.driver.ind.group1", "num.driver.ind.group2", "psub", "pind", "pglobal", "qsub", "qind", "qglobal","num.driver.mis.group1", "num.driver.mis.group2", "num.driver.trunc.group1", "num.driver.trunc.group2", "pmis", "ptrunc", "qmis", "qtrunc", paste0(single.test.names, "1"), paste0(single.test.names, "2"))
    full.data = lldf[,full.data.from.lldf]
    colnames(full.data) = full.data.from.lldf.names
    
    out.list = list()
    out.list[["substitutions"]] = get_effect_size(full.data, mutation.class = "sub")
    out.list[["indels"]] = get_effect_size(full.data, mutation.class = "ind")
    out.list[["missense_sub"]] = get_effect_size(full.data, mutation.class = "mis")
    out.list[["truncating_sub"]] = get_effect_size(full.data, mutation.class = "trunc")
    out.list[["overall_mut"]] = get_effect_size(full.data, mutation.class = "global")
    #out.list[["full"]] = full.data
    out.list[["dndscv"]] = list()
    out.list$dndscv[["dndscv_group1"]] = group1_dndsout
    out.list$dndscv[["dndscv_group2"]] = group2_dndsout

    # Return the output
    message("Done.")
    return(out.list)

  }, error=function(err) {

    print("An error has occurred.")
    print(err)
    return(NULL)

  }) # end tryCatch
}
