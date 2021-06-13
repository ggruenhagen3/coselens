#' COnditional SELection on the Excess of NonSynonymous Substitutions (coselenss)
#' 
#' Given an input of mutations from two groups of patients, dN/dS is calculated in each of the groups, 
#' and the excess of non-synonymous substitutions between them is found as an estimate for the difference 
#' in the number of driver mutations. For example, one group can be mutations from individuals in a certain
#' cancer type with mutations in a specific gene vs those without mutations in the specific gene.
#' A p value is returned for positive conditional selection in each gene in the human reference genome (hg19).
#' 
#' @param group1 group of individuals (for example those that contain a mutation in a split_gene)
#' @param group2 another group of individuals that do NOT contain a mutation in a split_gene
#' @param subset.genes.by genes to subset results by
#' @param refdb Reference database (path to .rda file)
#' 
#' @return coselenss returns a dataframe with rows representing and the following columns
#' @return - gene_name: name of gene that conditional selection was calculated in
#' @return - num.drivers.group1: estimate of the number of drivers in group 1 based excess of non-synonymous mutations
#' @return - num.drivers.group2: estimate of the number of drivers in group 2 based excess of non-synonymous mutations
#' @return - pmis: p-value for conditional selection in missense mutations
#' @return - ptrunc: p-value for conditional selection in truncating mutations
#' @return - pall: p-value for conditional selection in all substitutions (pmis and ptrunc)
#' @return - pind: p-value for conditional selection in small indels
#' @return - pglobal: Fisher's combined p-value for pall and pind
#' @return - qglobal: q-value of pglobal using Benjamini-Hochberg correction
#' 
#' @export

coselenss = function(group1, group2, subset.genes.by = NULL, refdb = "hg19") {
  # Find dN/dS and CIs
  tryCatch({
    # Calculate dN/dS and the confidence intervals
    group1_dndsout <- dndscv(group1, gene_list = NULL, refdb = refdb, sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = Inf, outmats=T, outmutrates = T, wg = F, ex = T)
    group2_dndsout <- dndscv(group2, gene_list = NULL, refdb = refdb, sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = Inf, outmats=T, outmutrates = T, wg = F, ex = T)
    
    # Calculate nonsynonymous mutation excess
    group1_ex <- calc_ex(group1_dndsout)
    group2_ex <- calc_ex(group2_dndsout)
    group1_ex$ex_tot <- group1_ex[,2] + group1_ex[,3]
    group2_ex$ex_tot <- group2_ex[,2] + group2_ex[,3]
    ex_values <- merge(x=group1_ex[,c("gene_name","ex_mis","ex_non","ex_tot")], y=group2_ex[,c("gene_name", "ex_mis", "ex_non", "ex_tot")], by ="gene_name", suffixes=c(".group1", ".group2"))
    
    # Put the data in the same order
    group1_sel_cv = group1_dndsout$sel_cv[order(as.numeric(rownames(group1_dndsout$sel_cv))),,drop=FALSE]
    group2_sel_cv = group2_dndsout$sel_cv[order(as.numeric(rownames(group2_dndsout$sel_cv))),,drop=FALSE]
    
    # If necessary correct for the difference in mutational loads between groups by normalizing dN/dS by
    # a global dN/dS in all cancer genes. 
    h0_sel_cv <<- lik_func(group1_dndsout, group2_dndsout)
    
    # Merge the data into one structure
    lldf <- merge(x=h0_sel_cv[,c("gene_name", "llall", "llmis", "lltrunc", "llmis_check_A", "lltrunc_check_B", "llind0", "llind1")], y=group1_sel_cv[,c("gene_name", "wmis_cv", "wnon_cv", "llall")], by="gene_name", suffixes=c(".0", ".group1"))
    lldf <- merge(x=lldf, y=group2_sel_cv[,c("gene_name", "wmis_cv", "wnon_cv", "llall")], by ="gene_name", suffixes=c(".group1", ".group2"))
    colnames(lldf)[c(ncol(lldf)-2, ncol(lldf)-1, ncol(lldf))] <- c("wmis_cv.group2", "wnon_cv.group2", "llall.group2")
    
    # Alternative Hypothesis Log-Likelihoods for SNV
    lldf$llall.1   <- lldf$llall.group1   + lldf$llall.group2
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
    lldf$qglobal = p.adjust(lldf$pglobal, method = "BH")
    
    # Add mutation excess data
    lldf <- merge(x = lldf, y = ex_values[,c("gene_name","ex_tot.group1","ex_tot.group2")], by = "gene_name")
    
    # Subset the columns returned and reorder them
    lldf = lldf[,c("gene_name", "ex_tot.group1", "ex_tot.group2", "pmis", "ptrunc", "pall", "pind", "pglobal", "qglobal")]
    colnames(lldf)[2:3] = c("num.drivers.group1", "num.drivers.group2")
    
    # Return the output
    print("Done.")
    return(lldf)
    
  }, error=function(err) {
    
    print("An error has occured.")
    print(err)
    return(NULL)
    
  }) # end tryCatch
}