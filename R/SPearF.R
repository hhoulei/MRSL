#' @title SPearF
#'
#' @description Provides Spearman's footrule to evaluate the estimated
#'      topological sorting of a graph
#'
#'
#' @param v1 The estimated topological sorting of a graph
#' @param v2 The true topological sorting of a graph
#'
#' @return A value of Spearman's footrule
#'
#'
#' @examples
#' v1 <- sample(1:10,10)
#' v2 <- 1:10
#' SPearF(v1,v2)
#'
#' @export
#'
SPearF <- function(v1,v2){
  SpearF <- sum(sapply(seq_along(v1), function(i) abs(i - (which(v2 == v1[i])))))
  return(SpearF)
}
