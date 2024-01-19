#' @title kenDTau
#'
#' @description Provides Kendall’s tau to evaluate the estimated
#'      topological sorting of a graph
#'
#'
#' @param v1 The estimated topological sorting of a graph
#' @param v2 The true topological sorting of a graph
#'
#' @return A value of Kendall’s tau
#'
#'
#' @examples
#' v1 <- sample(1:10,10)
#' v2 <- 1:10
#' kenDTau(v1,v2)
#'
#' @export
#'
kenDTau <- function(v1,v2){
  kenDcor <- cor.test(v1,v2, use = "pairwise", method="kendall")
  kenDcor_res <- c(kenDcor$statistic,kenDcor$p.value,kenDcor$estimate)
  return(kenDcor_res)
}
