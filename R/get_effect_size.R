#' Get effect sizes and classify instances of conditional selection
#'
#' @param coselens_full: full results table produced by coselens (coselens_out$full)
#' @param mutation.class: class of mutations for which effect sizes will be calculated. Options: "sub" (all coding substitutions, default), "ind" (indels), "mis" (missense substitutions), "trunc" (truncating substitutions, including nonsense and essential splice site substitutions), "global" (combination of coding substitutions and indels, using Fisher's combined test for p-values)
#'
#' @return dataframe with rows representing genes and the following columns
#' @return - gene_name: name of the gene
#' @return - num.driver.group1: estimate of the number of drivers per sample per gene in group 1
#' @return - num.driver.group2: estimate of the number of drivers per sample per gene in group 2
#' @return - Delta.Nd: absolute difference in the average number of driver mutations per sample (group 1 minus group 2)
#' @return - classification: classification of conditional selection. The most frequent classes are strict dependence (drivers only in group 1), facilitation (drivers more frequent in group 1), independence, inhibition (drivers less frequent in group 1), and strict inhibition (drivers absent from group 1). If negative selection is present, other possibilities are strict dependence with sign change (drivers positively selected in group 1 but negatively selected in group 2), strict inhibition with sign change (drivers positively selected in group 2 but negatively selected in group 1), aggravation (purifying selection against mutations becomes stronger in group 1), and relaxation (purifying selection against mutations becomes weaker in group 1).
#' @return - dependency: dependency index, measuring the association between the grouping variable (group 1 or 2) and the average number of drivers observed in a gene. It serves as a quantitative measure of the qualitative effect described in "classification". In the most common cases, a value of 1 indicates strict dependence or inhibition (drivers only observed in one group) and a value of 0 (or NA) indicates independence.
#' @return - pval: p-value for conditional selection
#' @return - qval: q-value for conditional selection using Benjamini-Hochberg correction of false discovery rate.
#'
#' @export

get_effect_size = function(coselens_full, mutation.class = "sub") {

  # Get relevant data for the desired class of mutations (written only for mutation.class = "sub")
  if (mutation.class == "sub") {
    ndr1 <- coselens_full$num.driver.sub.group1
    ndr2 <- coselens_full$num.driver.sub.group2
    pval <- coselens_full$psub
    qval <- coselens_full$qsub
    q1 <- coselens_full$qsub.group1
    q2 <- coselens_full$qsub.group2
  } else if (mutation.class == "ind") {
    ndr1 <- coselens_full$num.driver.ind.group1
    ndr2 <- coselens_full$num.driver.ind.group2
    pval <- coselens_full$pind
    qval <- coselens_full$qind
    q1 <- coselens_full$qind.group1
    q2 <- coselens_full$qind.group2
  } else if (mutation.class == "mis") {
    ndr1 <- coselens_full$num.driver.mis.group1
    ndr2 <- coselens_full$num.driver.mis.group2
    pval <- coselens_full$pmis
    qval <- coselens_full$qmis
    q1 <- coselens_full$qmis.group1
    q2 <- coselens_full$qmis.group2
  } else if (mutation.class == "trunc") {
    ndr1 <- coselens_full$num.driver.trunc.group1
    ndr2 <- coselens_full$num.driver.trunc.group2
    pval <- coselens_full$ptrunc
    qval <- coselens_full$qtrunc
    q1 <- coselens_full$qtrunc.group1
    q2 <- coselens_full$qtrunc.group2
  } else if (mutation.class == "global") {
    ndr1 <- coselens_full$num.driver.sub.group1 + coselens_full$num.driver.ind.group1
    ndr2 <- coselens_full$num.driver.sub.group2 + coselens_full$num.driver.ind.group2
    pval <- coselens_full$pglobal
    qval <- coselens_full$qglobal
    q1 <- coselens_full$qglobal.group1
    q2 <- coselens_full$qglobal.group2
  } else {
    stop("Invalid argument for mutation.class. Valid options are 'sub', 'ind', 'mis', and 'trunc'")
  }

  
  effectSz = data.frame(gene_name = coselens_full$gene_name)
  effectSz$num.driver.group1 = ndr1
  effectSz$num.driver.group2 = ndr2
  
  # Calculate Delta Nd
  Delta.Nd = ndr1 - ndr2
  effectSz$Delta.Nd = Delta.Nd

  # Classify
  f_independence      <- qval > 0.05  # independence
  f_strict_dependence <- qval < 0.05 & ndr1 > ndr2 & ndr1 > 0 & q1 < 0.05 & q2 > 0.05 # strict dependence
  f_facilitation      <- qval < 0.05 & ndr1 > ndr2 & ndr2 > 0 & q1 < 0.05 & q2 < 0.05 # facilitation
  f_inhibition        <- qval < 0.05 & ndr1 < ndr2 & ndr1 > 0 & q1 < 0.05 & q2 < 0.05 # inhibition
  f_strict_inhibition <- qval < 0.05 & ndr1 < ndr2 & ndr2 > 0 & q1 > 0.05 & q2 < 0.05 # strict inhibition
  f_strict_dep_signch <- qval < 0.05 & ndr1 > 0 & ndr2 < 0 & q1 < 0.05 & q2 < 0.05    # strict dependence with sign change
  f_strict_inh_signch <- qval < 0.05 & ndr1 < 0 & ndr2 > 0 & q1 < 0.05 & q2 < 0.05    # strict inhibition with sign change
  f_aggravation       <- qval < 0.05 & ndr1 < ndr2 & ndr2 < 0 & q1 < 0.05 & q2 < 0.05 # aggravation
  f_relaxation        <- qval < 0.05 & ndr1 > ndr2 & ndr1 < 0 & q1 < 0.05 & q2 < 0.05 # relaxation
  f_strict_dependence_alt <- qval < 0.05 & ndr1 > ndr2 & ndr1 > 0 & q1 > 0.05 & q2 > 0.05 # strict dependence (even if dndscv does not detect significant selection in group1)
  f_strict_inhibition_alt <- qval < 0.05 & ndr1 < ndr2 & ndr2 > 0 & q1 > 0.05 & q2 > 0.05 # strict inhibition (even if dndscv does not detect significant selection in group2)
  # f_unclass <- !(f_independence | f_strict_dependence | f_facilitation | f_inhibition | f_strict_inhibition | f_strict_dep_signch | f_strict_inh_signch | f_aggravation | f_relaxation) # unclassified

  sel_class <- rep("unclassified",length(qval))
  sel_class[which(f_independence)] = "independence"
  sel_class[which(f_strict_dependence)] = "strict dependence"
  sel_class[which(f_facilitation)] = "facilitation"
  sel_class[which(f_inhibition)] = "inhibition"
  sel_class[which(f_strict_inhibition)] = "strict inhibition"
  sel_class[which(f_strict_dep_signch)] = "strict dependence with sign change"
  sel_class[which(f_strict_inh_signch)] = "strict inhibition with sign change"
  sel_class[which(f_aggravation)] = "aggravation"
  sel_class[which(f_relaxation)] = "relaxation"
  #sel_class[which(f_strict_dependence_alt)] = "strict dependence"
  sel_class[which(f_strict_inhibition_alt)] = "strict inhibition"

  effectSz$classification = sel_class
  
  # Calculate Dependency: convert to polar, take the angular coordinate, and normalize
  dep_idx = geometry::cart2pol(ndr2,ndr1)[,1]*4/pi - 1
  dep_idx[which(f_independence)] = NA
  dep_idx[which(f_strict_dependence& dep_idx > 1)] = 1                              # Correct cases of strict dependence > 1
  dep_idx[which(f_strict_dependence_alt & dep_idx > 1)] = 1                         # Correct cases of strict dependence > 1
  dep_idx[which(f_inhibition)] = - dep_idx[which(f_inhibition)]                     # Change sign of inhibition
  dep_idx[which(f_strict_inhibition)] = - dep_idx[which(f_strict_inhibition)]       # Change sign of strict inhibition
  dep_idx[which(f_strict_inhibition_alt)] = - dep_idx[which(f_strict_inhibition_alt)]   # Change sign of strict inhibition (cases without dndscv significance)
  dep_idx[which(f_strict_inhibition & dep_idx > 1)] = 1                             # Correct cases of strict inhibition > 1
  dep_idx[which(f_strict_inhibition_alt & dep_idx > 1)] = 1                         # Correct cases of strict inhibition > 1
  dep_idx[which(f_strict_dep_signch)] = (3 - dep_idx[which(f_strict_dep_signch)])/2 # Rescale and reverse strict dependence with sign change
  dep_idx[which(f_strict_inh_signch)] = (dep_idx[which(f_strict_inh_signch)] + 3)/2 # Rescale and reverse strict inhibition with sign change
  dep_idx[which(f_aggravation)] = dep_idx[which(f_aggravation)] + 4                 # Rescale aggravation
  dep_idx[which(f_relaxation)] = - dep_idx[which(f_relaxation)] - 4                 # Rescale relaxation

  effectSz$dependency = dep_idx

  # Add p-val and q-val
  effectSz$pval = pval
  effectSz$qval = qval

  return(effectSz)
}
