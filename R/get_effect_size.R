#' Get effect sizes and classify instances of conditional selection
#'
#' @param coselens_full: full results table produced by coselens (coselens_out$full)
#' @param mutation.class: class of mutations for which effect sizes will be calculated. Options: "sub" (all coding substitutions, default), "ind" (indels), "mis" (missense substitutions), "trunc" (truncating substitutions)
#'
#' @return dataframe with rows representing genes and the following columns
#' @return - gene_name: name of the gene
#' @return - Delta.Nd: absolute difference in the average number of driver mutations per sample (group 1 minus group 2)
#' @return - dependency: dependency index, measuring the association between the grouping variable (group 1 or 2) and the average number of drivers observed in a gene. A value of 1 indicates strict dependence (driver mutations only observed in group 1); a value of 0 indicates independence; a value of -1 indicates strict inhibition (drivers observed in group 2 but never in group 1).
#' @return - class: classification of conditional selection. The most frequent classes are strict dependence (drivers only in group 1), facilitation (drivers more frequent in group 1), independence, partial inhibition (drivers less frequent in group 1), and strict inhibition (drivers absent from group 1). If negative selection is present, other possibilities are strict dependence with sign change (drivers positively selected in group 1 but negatively selected in group 2), strict inhibition with sign change (drivers positively selected in group 2 but negatively selected in group 1), aggravation (purifying selection against mutations becomes stronger in group 1), and relaxation (purifying selection against mutations becomes weaker in group 1).
#' @return - pval: p-value for conditional selection
#' @return - qval: q-value for conditional selection
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
  } else {
    stop("Invalid argument for mutation.class. Valid options are 'sub', 'ind', 'mis', and 'trunc'")
  }

  # Calculate Delta Nd
  Delta.Nd = ndr1 - ndr2
  effectSz = data.frame(Delta.Nd = Delta.Nd)

  # Calculate Dependency: convert to polar, take the angular coordinate, and normalize
  effectSz$dependency = geometry::cart2pol(ndr2,ndr1)[,1]*4/pi - 1

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

  effectSz$classification = sel_class

  # Add p-val and q-val
  effectSz$pval = pval
  effectSz$qval = qval

  return(effectSz)
}
