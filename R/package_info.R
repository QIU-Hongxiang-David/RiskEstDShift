#' @importFrom assertthat assert_that is.string has_name noNA is.count is.flag is.number
#' @importFrom DescTools BinomCI
#' @import SuperLearner
NULL


#' @title RiskEstDShift package
#' 
#' @description
#' In this package, we implement efficient and multiply robust cross-fit target population risk estimators under one of the following four dataset shift conditions:
#' * full-data covariate shift: The distribution of \eqn{Y|X}{Y|X} is same in source and target populations and \eqn{(X,Y)}{(X,Y)} is observed in data from both populations
#' * full-data label shift: The distribution of \eqn{X|Y}{X|Y} is same in source and target populations and \eqn{(X,Y)}{(X,Y)} is observed in data from both populations
#' * concept shift in the features: The distribution of \eqn{X}{X} is same in source and target populations
#' * concept shift in the labels: The distribution of \eqn{Y}{Y} is same in source and target populations
#' 
#' @docType package
#' @name RiskEstDShift-package
NULL
